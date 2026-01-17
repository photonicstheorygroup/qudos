// QUantum Dot Open-source Simulator (QUDOS)
//
// Copyright (c) 2025
// Tommy Murphy and Christopher A. Broderick
//
// This file is part of QUDOS.
//
// QUDOS is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either
// version 3 of the License, or (at your discretion) any later version.
//
// QUDOS is distributed in the hope that it will be useful, but without any warranty and
// without an implied warranty of merchantability or fitness for a particular purpose.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with QUDOS.
// If not, see <https://www.gnu.org/licenses/>.
//
// Alternative commercial licensing terms may be available from the copyright holders upon request.

#include "8_band.h"

void piezo_reciprocal_space_fft( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell )
{

	// Process limit for Dfti FFT library
	MKL_INT process_limit = ( 2*strain_basis.N_Gx_T - 1 <= 2*strain_basis.N_Gy_T - 1 ? 2*strain_basis.N_Gx_T - 1 : 2*strain_basis.N_Gy_T - 1 );

	// Sub-communicator splits communicator into 2 for process IDs under limit and above limit
	MPI_Comm MPI_SUB_COMM;
	if ( nprocs > process_limit )
	{

		MPI_Comm_split( MPI_COMM_WORLD, ( myid < process_limit ? 0 : 1 ), myid, &MPI_SUB_COMM );

	}
	else
	{

		MPI_Comm_dup( MPI_COMM_WORLD, &MPI_SUB_COMM );

	}

	// Variables relating to the extraction of the computationally relevant central portion of the calculated convolution
	MKL_INT local_start_reduced;
	MKL_INT local_rows_reduced;
	MKL_INT local_size_reduced;

	// Local electronic supercell characteristic function
	struct characteristic_functions characteristic_function_G_local;
	// Local strain supercell piezoelectric potential
	struct piezo_potential phi_piezo_G_local;
	// Local piezoelectric potential and electronic supercell characteristic function convolution
	struct convolution_char_piezo convolutions_char_piezo_local;
	// Central (reduced) portion of local piezoelectric potential and electronic supercell characteristic function convolution
	struct convolution_char_piezo convolutions_char_piezo_reduced_local;

	// The following store information about the local reduced convolution array in order to later gather it on a single process
	int *recvcounts = ( int * ) malloc( nprocs*sizeof( int ) );
	int *displs = ( int * ) malloc( nprocs*sizeof( int ) );

	for ( int proc = 0; proc < nprocs; proc++ )
	{

		recvcounts[proc] = 1;

		displs[proc] = ( proc == 0 ? 0 : displs[proc - 1] + recvcounts[proc - 1] );
	
	}

	int *recvcounts_convolutions = ( int * ) malloc( nprocs*sizeof( int ) );
	int *displs_convolutions = ( int * ) malloc( nprocs*sizeof( int ) );	


	if ( myid < process_limit ) // Processes with IDs less than the process limit
	{

		// Total number of elements in each dimension for zero-padded arrays involved in Dfti FFT calculations as part of convolution calculation
		MKL_LONG sizes[] = { 2*strain_basis.N_Gx_T - 1, 2*strain_basis.N_Gy_T - 1, 2*strain_basis.N_Gz_T - 1 };

		// Global size of involved arrays
		MKL_LONG global_size = sizes[0]*sizes[1]*sizes[2];

		// Dfti precision & array type
		MKL_LONG DFTI_PRECISION = DFTI_DOUBLE;
		MKL_LONG DFTI_TYPE = DFTI_COMPLEX;

		// Dfti task descriptor
		DFTI_DESCRIPTOR_DM_HANDLE task;

		// Initialise Dfti task descriptor
		check_dfti( DftiCreateDescriptorDM( MPI_SUB_COMM, &task, DFTI_PRECISION, DFTI_TYPE, 3, sizes ), "DftiCreateDescriptorDM" );

		// Set Dfti backward scale to represent inverse FFT
		check_dfti( DftiSetValueDM( task, DFTI_BACKWARD_SCALE, 1.0 / global_size ), "DftiSetValueDM" );

		// Local memory size for arrays involved in Dfti FFT calculations
		MKL_LONG local_memory_size;
		check_dfti( DftiGetValueDM( task, CDFT_LOCAL_SIZE, &local_memory_size ), "DftiGetValueDM" );
		// Global position of first element of local arrays involved in Dfti FFT calculations
		MKL_LONG local_start;
		check_dfti( DftiGetValueDM( task, CDFT_LOCAL_X_START, &local_start ), "DftiGetValueDM" );
		// Number of global rows in local arrays involved in Dfti FFT calculations
		MKL_LONG local_rows;
		check_dfti( DftiGetValueDM( task, CDFT_LOCAL_NX, &local_rows ), "DftiGetValueDM" );
		// Number of elements in local arrays involved in Dfti FFT calculations
		MKL_LONG local_size = local_rows*sizes[1]*sizes[2];

		// Commit Dfti descriptor
		check_dfti( DftiCommitDescriptorDM( task ), "DftiCommitDescriptorDM" );


		// Will hold local part of electronic structure supercell characteristic function
		characteristic_function_G_local.chi = ( complex double* ) calloc( local_memory_size, sizeof( complex double ) );
		// Will hold local part of strain supercell piezoelectric potential
		phi_piezo_G_local.pz = ( complex double* ) calloc( local_memory_size, sizeof( complex double ) );

		// Extracting local parts of zero-padded electronic structure supercell characteristic function and piezoelectric potential from global
		#pragma omp parallel for collapse( 3 )
		for ( int global_row = local_start; global_row < local_start + local_rows; global_row++ )
		{

			for ( int idx_y = 0; idx_y < sizes[1]; idx_y++ )
			{

				for ( int idx_z = 0; idx_z < sizes[2]; idx_z++ )
				{

					if ( global_row >= 0 && global_row < strain_basis.N_Gx_T && idx_y >= 0 && idx_y <  strain_basis.N_Gy_T && idx_z >= 0 && idx_z < strain_basis.N_Gz_T )
					{

						characteristic_function_G_local.chi[( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = characteristic_function_strain_G[0].chi[ ( global_row )*( strain_basis.N_Gy_T )*( strain_basis.N_Gz_T ) + ( idx_y )*( strain_basis.N_Gz_T ) + ( idx_z ) ];
						
						phi_piezo_G_local.pz[( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = phi_piezo_region_G.pz[ ( global_row )*( strain_basis.N_Gy_T )*( strain_basis.N_Gz_T ) + ( idx_y )*( strain_basis.N_Gz_T ) + ( idx_z )];

					}

				}

			}

		}

		// Compute FFT using Dfti library of electronic structure supercell characteristic function
		check_dfti( DftiComputeForwardDM( task, characteristic_function_G_local.chi ), "DftiComputeForwardDM" );

		// Compute FFT using Dfti library of piezoelectric potential
		check_dfti( DftiComputeForwardDM( task, phi_piezo_G_local.pz ), "DftiComputeForwardDM" );


		// Setting local variables enabling extraction of central portion of convolved array, i.e. number of components of central portion that exist on each process and start point in local data
		if ( local_start >= strain_basis.N_Gx_T - 1 - supercell[0].RL1*( basis.N_Gx_T / 2 ) && local_start < strain_basis.N_Gx_T + supercell[0].RL1*( basis.N_Gx_T / 2 ) )
		{

			local_start_reduced = ( ( local_start - ( strain_basis.N_Gx_T - 1 - supercell[0].RL1*( basis.N_Gx_T / 2 ) ) ) % supercell[0].RL1 == 0 ? 0 : supercell[0].RL1 - ( local_start - ( strain_basis.N_Gx_T - 1 - supercell[0].RL1*( basis.N_Gx_T / 2 ) ) ) % supercell[0].RL1 );
			if ( local_start + local_rows < strain_basis.N_Gx_T + supercell[0].RL1*( basis.N_Gx_T / 2 ) )
			{

				local_rows_reduced = ( local_rows / supercell[0].RL1 ) + ( local_start_reduced < local_rows % supercell[0].RL1 ? 1 : 0 );

			}
			else
			{

				local_rows_reduced = ( ( strain_basis.N_Gx_T + supercell[0].RL1*( basis.N_Gx_T / 2 ) - local_start ) / supercell[0].RL1 ) + ( local_start_reduced < ( strain_basis.N_Gx_T + supercell[0].RL1*( basis.N_Gx_T / 2 ) - local_start ) % supercell[0].RL1 ? 1 : 0 );

			}

		}
		else if ( local_start < strain_basis.N_Gx_T - 1 - supercell[0].RL1*( basis.N_Gx_T / 2 ) )
		{

			if ( local_start + local_rows >= strain_basis.N_Gx_T - 1 - supercell[0].RL1*( basis.N_Gx_T / 2 ) && local_start + local_rows < strain_basis.N_Gx_T - 1 + supercell[0].RL1*( basis.N_Gx_T / 2 ) )
			{

				local_start_reduced = strain_basis.N_Gx_T - 1 - supercell[0].RL1*( basis.N_Gx_T / 2 ) - local_start;
				local_rows_reduced = ( ( local_rows - local_start_reduced ) / supercell[0].RL1 ) + ( 0 < ( local_rows - local_start_reduced ) % supercell[0].RL1 ? 1 : 0 );

			}
			else if ( local_start + local_rows >= strain_basis.N_Gx_T - 1 + supercell[0].RL1*( basis.N_Gx_T / 2 ) )
			{

				local_start_reduced = strain_basis.N_Gx_T - 1 - supercell[0].RL1*( basis.N_Gx_T / 2 ) - local_start;
				local_rows_reduced = basis.N_Gx_T;

			}
			else
			{

				local_start_reduced = 0;
				local_rows_reduced = 0;

			}

		}
		else
		{

			local_start_reduced = 0;
			local_rows_reduced = 0;

		}
		local_size_reduced = local_rows_reduced*basis.N_Gy_T*basis.N_Gz_T;

		// Gathering information regarding local number of elements in reduced array on each process that must be sent to all processes
		MPI_Allgatherv( &local_size_reduced, 1, MPI_INT, recvcounts_convolutions, recvcounts, displs, MPI_INT, MPI_COMM_WORLD );

		// Calulating reduced array global displacements based on gathered data
		for ( int proc = 0; proc < nprocs; proc++ )
		{

			displs_convolutions[proc] = ( proc == 0 ? 0 : displs_convolutions[proc - 1] + recvcounts_convolutions[proc - 1] );
		
		}


		// Local array to store element by element multiplcation of FFT arrays
		convolutions_char_piezo_local.char_piezo = ( complex double* ) calloc( local_memory_size, sizeof( complex double ) );

		// Local array to then store local extracted central portion
		convolutions_char_piezo_reduced_local.char_piezo = ( complex double* ) calloc( local_size_reduced, sizeof( complex double ) );

		// Element by element multiplcation of FFT arrays
		vzMul( local_size, characteristic_function_G_local.chi, phi_piezo_G_local.pz, convolutions_char_piezo_local.char_piezo );


		// Compute inverse FFT using Dfti library of element by element multiplied arrays to complete convolution calculation
		check_dfti( DftiComputeBackwardDM( task, convolutions_char_piezo_local.char_piezo ), "DftiComputeBackwardDM" );
		

		free( characteristic_function_G_local.chi );

		free( phi_piezo_G_local.pz );


		// Extracting central portion of local convolution array
		#pragma omp parallel for collapse( 3 )
		for ( int idx_x = 0; idx_x < local_rows_reduced; idx_x++ )
		{

			for ( int idx_y = 0; idx_y < basis.N_Gy_T; idx_y++ )
			{

				for ( int idx_z = 0; idx_z < basis.N_Gz_T; idx_z++ )
				{

					convolutions_char_piezo_reduced_local.char_piezo[idx_x*basis.N_Gy_T*basis.N_Gz_T + idx_y*basis.N_Gz_T + idx_z] = convolutions_char_piezo_local.char_piezo[ ( supercell[0].RL1*idx_x + local_start_reduced )*( sizes[1] )*( sizes[2] ) + ( supercell[0].RL2*idx_y + strain_basis.N_Gy_T - 1 - supercell[0].RL2*( basis.N_Gy_T / 2 ) )*( sizes[2] ) + ( supercell[0].RL3*idx_z + strain_basis.N_Gz_T - 1 - supercell[0].RL3*( basis.N_Gz_T / 2 ) ) ];

				}

			}

		}


		// Gathering central portion of local convolution arrays
		MPI_Allgatherv( convolutions_char_piezo_reduced_local.char_piezo, local_size_reduced, MPI_DOUBLE_COMPLEX, phi_piezo_G.pz, recvcounts_convolutions, displs_convolutions, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
		

		check_dfti( DftiFreeDescriptorDM( &task ), "DftiFreeDescriptorDM" );

		free( convolutions_char_piezo_local.char_piezo );

		free( convolutions_char_piezo_reduced_local.char_piezo );

	}
	else // Processes with IDs greater than or equal to the process limit
	{

		local_start_reduced = 0;
		local_rows_reduced = 0;
		local_size_reduced = 0;


		// Gathering information on number of local reduced convolution array elements that will be sent by each process
		MPI_Allgatherv( &local_size_reduced, 1, MPI_INT, recvcounts_convolutions, recvcounts, displs, MPI_INT, MPI_COMM_WORLD );

		for ( int proc = 0; proc < nprocs; proc++ )
		{

			displs_convolutions[proc] = ( proc == 0 ? 0 : displs_convolutions[proc - 1] + recvcounts_convolutions[proc - 1] );
		
		}


		// Gathering local reduced convolution array elements
		MPI_Allgatherv( convolutions_char_piezo_reduced_local.char_piezo, local_size_reduced, MPI_DOUBLE_COMPLEX, phi_piezo_G.pz, recvcounts_convolutions, displs_convolutions, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );

	}

	
	free( recvcounts );
	free( recvcounts_convolutions );

	free( displs );
	free( displs_convolutions );

}
