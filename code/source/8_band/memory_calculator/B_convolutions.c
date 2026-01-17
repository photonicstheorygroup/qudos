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

void B_convolutions( struct wave_vectors basis, struct wave_vectors coulomb_basis, struct supercell_geometry *supercell, int *N_states, struct eigenstate *eigenstate_vb_padded, struct eigenstate *eigenstate_cb_padded, struct eigenstate *eigenstate_vb_padded_conj, struct eigenstate *eigenstate_cb_padded_conj, MKL_INT process_num, MKL_INT process_limit, DFTI_DESCRIPTOR_DM_HANDLE task, MKL_LONG local_memory_size, MKL_LONG local_size, struct convolution_state *convolution_states, struct convolution_state *convolution_states_reversed, int *sendcounts, int *displs_send, int coulomb_region_type )
{

	// Variables for gathering B Fourier coefficients after initial calculation using electronic structure supercell
	
	int N_G_conv = ( 2*basis.N_Gx - 1 )*( 2*basis.N_Gy - 1 )*( 2*basis.N_Gz - 1 );
	int N_GC_conv = ( 2*coulomb_basis.N_Gx - 1 )*( 2*coulomb_basis.N_Gy - 1 )*( 2*coulomb_basis.N_Gz - 1 );
	
	int *recvcounts;// = ( int * ) malloc( nprocs*sizeof( int ) );
	int *displs;// = ( int * ) malloc( nprocs*sizeof( int ) );

	for ( int proc = 0; proc < nprocs; proc++ )
	{

		//recvcounts[proc] = 1;

		//displs[proc] = ( proc == 0 ? 0 : displs[proc - 1] + recvcounts[proc - 1] );
	
	}

	int *recvcounts_convolutions;// = ( int * ) malloc( nprocs*sizeof( int ) );
	int *displs_convolutions;// = ( int * ) malloc( nprocs*sizeof( int ) );

	//MPI_Allgatherv( &local_size, 1, MPI_INT, recvcounts_convolutions, recvcounts, displs, MPI_INT, MPI_COMM_WORLD );

	for ( int proc = 0; proc < nprocs; proc++ )
	{

		//displs_convolutions[proc] = ( proc == 0 ? 0 : displs_convolutions[proc - 1] + recvcounts_convolutions[proc - 1] );
	
	}


	// Variables for calculating B Fourier coefficients for Coulomb calculation supercell if using coulomb_region_type == 0
	
	// Local size of arrays
	int local_size_coulomb = coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T / nprocs + ( ( myid < coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T % nprocs ) ? 1 : 0 );
	
	struct convolution_state convolution_states_interim_local;
	int *recvcounts_coulomb;
	int *displs_coulomb; 
	
	if ( coulomb_region_type == 0 )
	{

		//convolution_states_interim_local.conv_state = ( complex double* ) calloc( local_size_coulomb, sizeof( complex double ) );
		local_mem += ( local_size_coulomb )*sizeof( complex double );	// Adding memory to store local parts of downsampled B array
		
		//recvcounts_coulomb = ( int * ) malloc( nprocs*sizeof( int ) );
		//displs_coulomb = ( int * ) malloc( nprocs*sizeof( int ) );

		// Will store displacement from first value of the "global" characteristic function for the first value of the characteristic functions calculated by each process
	
		for ( int proc = 0; proc < nprocs; proc++ )
		{
		
			//recvcounts_coulomb[proc] = coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T / nprocs + ( ( proc < coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T % nprocs ) ? 1 : 0 );
			//displs_coulomb[proc] = ( proc == 0 ? 0 : displs_coulomb[proc - 1] + recvcounts_coulomb[proc - 1] ); 
		
		}

	}


	// Variables for calculating B Fourier coefficients for Coulomb calculation supercell if using coulomb_region_type == 1

	struct characteristic_functions characteristic_function_coulomb_G;
	
	MKL_INT process_limit_B_supercell = ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) <= coulomb_basis.N_Gy_T + 2*supercell[0].RCL2*( basis.N_Gy_T / 2 ) ? coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) : coulomb_basis.N_Gy_T + 2*supercell[0].RCL2*( basis.N_Gy_T / 2 ) );
	MKL_INT process_num_B_supercell;

	MPI_Comm MPI_SUB_COMM_B_supercell;

	struct characteristic_functions characteristic_function_coulomb_G_local;
	struct convolution_char_state convolution_char_state_interim_local;
	struct convolution_char_state convolution_char_state_interim_reduced_local;

	int *recvcounts_convolutions_B_supercell;
	int *displs_convolutions_B_supercell;

	if ( coulomb_region_type == 1 )
	{

		if ( nprocs > process_limit_B_supercell )
		{

			MPI_Comm_split( MPI_COMM_WORLD, ( myid < process_limit_B_supercell ? 0 : 1 ), myid, &MPI_SUB_COMM_B_supercell );
			process_num_B_supercell = process_limit_B_supercell;

		}
		else
		{

			MPI_Comm_dup( MPI_COMM_WORLD, &MPI_SUB_COMM_B_supercell );
			process_num_B_supercell = nprocs;

		}

		if ( myid < process_num_B_supercell )
		{

			//recvcounts_convolutions_B_supercell = ( int * ) malloc( process_num_B_supercell*sizeof( int ) );
			//displs_convolutions_B_supercell = ( int * ) malloc( process_num_B_supercell*sizeof( int ) );

		}

	}

	MKL_LONG status_B_supercell;

	MKL_LONG sizes_B_supercell[] = { coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ), coulomb_basis.N_Gy_T + 2*supercell[0].RCL2*( basis.N_Gy_T / 2 ), coulomb_basis.N_Gz_T + 2*supercell[0].RCL3*( basis.N_Gz_T / 2 ) };

	MKL_LONG global_size_B_supercell = sizes_B_supercell[0]*sizes_B_supercell[1]*sizes_B_supercell[2];

	MKL_LONG DFTI_PRECISION_B_supercell = DFTI_DOUBLE;
	MKL_LONG DFTI_TYPE_B_supercell = DFTI_COMPLEX;

	DFTI_DESCRIPTOR_DM_HANDLE task_B_supercell;

	MKL_LONG local_memory_size_B_supercell;
	
	MKL_LONG local_start_B_supercell;
	
	MKL_LONG local_rows_B_supercell;
	
	MKL_LONG local_size_B_supercell;


	MKL_INT local_start_reduced_B_supercell;
	MKL_INT local_rows_reduced_B_supercell;
	MKL_INT local_size_reduced_B_supercell;


	if ( coulomb_region_type == 1 )
	{

		//characteristic_function_coulomb_G.chi = ( complex double* ) calloc( coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T, sizeof( complex double ) );

		local_mem += ( coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T )*sizeof( complex double ); // Adding memory to store characteristic function of electronic structure supercell in Coulomb supercell

		//characteristic_function_coulomb_supercell( coulomb_basis, supercell, characteristic_function_coulomb_G );

		if ( myid < process_limit_B_supercell )
		{

			status_B_supercell = DftiCreateDescriptorDM( MPI_SUB_COMM_B_supercell, &task_B_supercell, DFTI_PRECISION_B_supercell, DFTI_TYPE_B_supercell, 3, sizes_B_supercell );

			status_B_supercell = DftiSetValueDM( task_B_supercell, DFTI_BACKWARD_SCALE, 1.0 / global_size_B_supercell );

			status_B_supercell = DftiGetValueDM( task_B_supercell, CDFT_LOCAL_SIZE, &local_memory_size_B_supercell );
			
			status_B_supercell = DftiGetValueDM( task_B_supercell, CDFT_LOCAL_X_START, &local_start_B_supercell );
			
			status_B_supercell = DftiGetValueDM( task_B_supercell, CDFT_LOCAL_NX, &local_rows_B_supercell );

			local_size_B_supercell = local_rows_B_supercell*sizes_B_supercell[1]*sizes_B_supercell[2];

			//status_B_supercell = DftiCommitDescriptorDM( task_B_supercell );


			//characteristic_function_coulomb_G_local.chi = ( complex double* ) calloc( local_memory_size_B_supercell, sizeof( complex double ) );
			//convolution_states_interim_local.conv_state = ( complex double* ) calloc( local_memory_size_B_supercell, sizeof( complex double ) );
			//convolution_char_state_interim_local.conv_char_state = ( complex double* ) calloc( local_memory_size_B_supercell, sizeof( complex double ) );

			local_mem += ( 3*local_memory_size_B_supercell )*sizeof( complex double ); // Adding memory to store local electronic supercell characteristic function, interim sum of local components of convolutions between B of electronic structure supercell and electronic structure characteristic function, and local B array and electronic structure supercell characteristic function convolutions

			#pragma omp parallel for collapse( 3 )
			for ( int global_row = local_start_B_supercell; global_row < local_start_B_supercell + local_rows_B_supercell; global_row++ )
			{

				for ( int idx_y = 0; idx_y < sizes_B_supercell[1]; idx_y++ )
				{

					for ( int idx_z = 0; idx_z < sizes_B_supercell[2]; idx_z++ )
					{

						if ( global_row >= 0 && global_row < coulomb_basis.N_Gx_T && idx_y >= 0 && idx_y < coulomb_basis.N_Gy_T && idx_z >= 0 && idx_z < coulomb_basis.N_Gz_T )
						{

							//characteristic_function_coulomb_G_local.chi[( global_row - local_start_B_supercell )*sizes_B_supercell[1]*sizes_B_supercell[2] + idx_y*sizes_B_supercell[2] + idx_z] = characteristic_function_coulomb_G.chi[ ( global_row )*( coulomb_basis.N_Gy_T )*( coulomb_basis.N_Gz_T ) + ( idx_y )*( coulomb_basis.N_Gz_T ) + ( idx_z ) ];

						}

					}

				}

			}

			//status_B_supercell = DftiComputeForwardDM( task_B_supercell, characteristic_function_coulomb_G_local.chi );

			if ( local_start_B_supercell >= ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / 2 - coulomb_basis.N_Gx_T / 2 && local_start_B_supercell < ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / 2 - coulomb_basis.N_Gx_T / 2 + coulomb_basis.N_Gx_T )
			{

				local_start_reduced_B_supercell = 0;
				if ( local_start_B_supercell + local_rows_B_supercell < ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / 2 - coulomb_basis.N_Gx_T / 2 + coulomb_basis.N_Gx_T )
				{

					local_rows_reduced_B_supercell = local_rows_B_supercell;

				}
				else
				{

					local_rows_reduced_B_supercell = ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / 2 - coulomb_basis.N_Gx_T / 2 + coulomb_basis.N_Gx_T - local_start_B_supercell;

				}

			}
			else if ( local_start_B_supercell < ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / 2 - coulomb_basis.N_Gx_T / 2 )
			{

				if ( local_start_B_supercell + local_rows_B_supercell >= ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / 2 - coulomb_basis.N_Gx_T / 2  && local_start_B_supercell + local_rows_B_supercell < ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / 2 - coulomb_basis.N_Gx_T / 2 + coulomb_basis.N_Gx_T )
				{

					local_start_reduced_B_supercell = ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / 2 - coulomb_basis.N_Gx_T / 2 - local_start_B_supercell;
					local_rows_reduced_B_supercell = local_rows_B_supercell - local_start_reduced_B_supercell;

				}
				else if ( local_start_B_supercell + local_rows_B_supercell >= ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / 2 - coulomb_basis.N_Gx_T / 2 + coulomb_basis.N_Gx_T )
				{

					local_start_reduced_B_supercell = ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / 2 - coulomb_basis.N_Gx_T / 2 - local_start_B_supercell;
					local_rows_reduced_B_supercell = coulomb_basis.N_Gx_T;

				}
				else
				{

					local_start_reduced_B_supercell = 0;
					local_rows_reduced_B_supercell = 0;

				}

			}
			else
			{

				local_start_reduced_B_supercell = 0;
				local_rows_reduced_B_supercell = 0;

			}
			local_size_reduced_B_supercell = local_rows_reduced_B_supercell*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T;

			//MPI_Allgatherv( &local_size_reduced_B_supercell, 1, MPI_INT, recvcounts_convolutions_B_supercell, recvcounts, displs, MPI_INT, MPI_SUB_COMM_B_supercell );

			for ( int proc = 0; proc < process_num_B_supercell; proc++ )
			{

				//displs_convolutions_B_supercell[proc] = ( proc == 0 ? 0 : displs_convolutions_B_supercell[proc - 1] + recvcounts_convolutions_B_supercell[proc - 1] );
			
			}

			//convolution_char_state_interim_reduced_local.conv_char_state = ( complex double* ) calloc( local_size_reduced_B_supercell, sizeof( complex double ) );

			local_mem += ( local_size_reduced_B_supercell )*sizeof( complex double ); // Adding memory to store local B array and electronic structure supercell characteristic function convolutions

			if ( local_mem > local_mem_peak )
			{
				
				local_mem_peak = local_mem; // Store if peak memory 

			}

		}

		//free( characteristic_function_coulomb_G.chi );
		local_mem -= ( coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T )*sizeof( complex double );

	}

	
	
	
	if ( myid < process_limit )
	{

		struct convolution_state convolution_states_local_bands;
		struct convolution_state convolution_states_local_interim;
		struct convolution_state convolution_states_local;
		struct convolution_state convolution_states_interim;
		struct convolution_state convolution_states_global;
		struct convolution_state convolution_states_global_reversed;

		//convolution_states_local_bands.conv_state = calloc( N_BAND*local_memory_size, sizeof( complex double ) );
		//convolution_states_local_interim.conv_state = calloc( local_size, sizeof( complex double ) );
		//convolution_states_local.conv_state = calloc( local_size, sizeof( complex double ) );

		//convolution_states_interim.conv_state = calloc( N_G_conv, sizeof( complex double ) );

		//convolution_states_global.conv_state = calloc( N_GC_conv, sizeof( complex double ) );
		//convolution_states_global_reversed.conv_state = calloc( N_GC_conv, sizeof( complex double ) );

		local_mem += ( N_BAND*local_memory_size )*sizeof( complex double ); // Adds memory to store local component of convolution between envelope functions for each individual band, to later calculate B
		local_mem += ( 2*local_size )*sizeof( complex double ); // Adds memory to store interim sum of local components of convolutions between envelope function band components, to later calculate B, and memory to store total sum of local components of convolutions between envelope function band components, to calculate B
		local_mem += ( N_G_conv )*sizeof( complex double ); // Adds memory to store global B for electronic structure supercell
		local_mem += ( 2*N_GC_conv )*sizeof( complex double );	// Adds memory to store global B for Coulomb supercell, and wave vector reversed version of global B for Coulomb supercell


		if ( local_mem > local_mem_peak )
		{
			
			local_mem_peak = local_mem; // Store if peak memory

		}

		
		//free( convolution_states_local_bands.conv_state );
		//free( convolution_states_local_interim.conv_state );
		//free( convolution_states_local.conv_state );
		//free( convolution_states_interim.conv_state );
		//free( convolution_states_global.conv_state );
		//free( convolution_states_global_reversed.conv_state );

		local_mem -= ( N_BAND*local_memory_size )*sizeof( complex double );
		local_mem -= ( 2*local_size )*sizeof( complex double );
		local_mem -= ( 2*N_GC_conv )*sizeof( complex double );
		local_mem -= ( N_G_conv )*sizeof( complex double );

	}
	else
	{

		struct convolution_state convolution_states_interim;
		struct convolution_state convolution_states_global;
		struct convolution_state convolution_states_global_reversed;

		//convolution_states_interim.conv_state = calloc( N_G_conv, sizeof( complex double ) );
		//convolution_states_global.conv_state = calloc( N_GC_conv, sizeof( complex double ) );

		struct convolution_state convolution_states_local;
		

		local_mem += ( N_G_conv )*sizeof( complex double ); // Adds memory to store global B for electronic structure supercell
		local_mem += ( N_GC_conv )*sizeof( complex double ); // Adds memory to store global B for Coulomb supercell

		if ( local_mem > local_mem_peak )
		{
			
			local_mem_peak = local_mem;

		}

		//free( convolution_states_interim.conv_state );
		//free( convolution_states_global.conv_state );

		local_mem -= ( N_GC_conv )*sizeof( complex double );
		local_mem -= ( N_G_conv )*sizeof( complex double );

	}

	//free( recvcounts );
	//free( displs );
	
	//free( recvcounts_convolutions );
	//free( displs_convolutions );

	if ( coulomb_region_type == 0 )
	{

		//free( convolution_states_interim_local.conv_state );
		local_mem -= ( local_size_coulomb )*sizeof( complex double );
		
		//free( recvcounts_coulomb );
		//free( displs_coulomb );

	}
	else if ( coulomb_region_type == 1 )
	{

		if ( myid < process_limit_B_supercell )
		{

			//DftiFreeDescriptorDM( &task_B_supercell );
			//free( convolution_states_interim_local.conv_state );
			//free( characteristic_function_coulomb_G_local.chi );
			//free( convolution_char_state_interim_local.conv_char_state );
			//free( convolution_char_state_interim_reduced_local.conv_char_state );

			//free( recvcounts_convolutions_B_supercell );
			//free( displs_convolutions_B_supercell );

			local_mem -= ( 3*local_memory_size_B_supercell )*sizeof( complex double );
			local_mem -= ( local_size_reduced_B_supercell )*sizeof( complex double );

		}

	}

}
