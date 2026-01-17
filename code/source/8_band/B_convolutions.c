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

	MPI_Allgatherv( &local_size, 1, MPI_INT, recvcounts_convolutions, recvcounts, displs, MPI_INT, MPI_COMM_WORLD );

	for ( int proc = 0; proc < nprocs; proc++ )
	{

		displs_convolutions[proc] = ( proc == 0 ? 0 : displs_convolutions[proc - 1] + recvcounts_convolutions[proc - 1] );
	
	}


	// Set-up for later calculating B Fourier coefficients for Coulomb calculation supercell if using free Coulomb supercell lengths
	
	// Local size of B array for Coulomb supercell
	int local_size_coulomb = coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T / nprocs + ( ( myid < coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T % nprocs ) ? 1 : 0 );
	
	// Will store interim sum of local components of convolutions between B of electronic structure supercell and electronic structure characteristic function
	struct convolution_state convolution_states_interim_local;
	
	// The following store information about the local sizes and displacements of the Coulomb supercell B array in order to later gather it on a single process
	int *recvcounts_coulomb;
	int *displs_coulomb; 
	
	if ( coulomb_region_type == 0 ) // Free Coulomb supercell lengths
	{

		convolution_states_interim_local.conv_state = ( complex double* ) calloc( local_size_coulomb, sizeof( complex double ) );
		recvcounts_coulomb = ( int * ) malloc( nprocs*sizeof( int ) );
		displs_coulomb = ( int * ) malloc( nprocs*sizeof( int ) );
	
		for ( int proc = 0; proc < nprocs; proc++ )
		{
		
			recvcounts_coulomb[proc] = coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T / nprocs + ( ( proc < coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T % nprocs ) ? 1 : 0 );
			displs_coulomb[proc] = ( proc == 0 ? 0 : displs_coulomb[proc - 1] + recvcounts_coulomb[proc - 1] ); 
		
		}

	}


	// Set-up for later calculating B Fourier coefficients for Coulomb calculation supercell if using Coulomb supercell lengths commensurate with those of the electronic structure supercell

	// Will store Fourier transform of electronic structure supercell characteristic function in Coulomb supercell
	struct characteristic_functions characteristic_function_coulomb_G;
	
	// Process limit for Dfti FFT library
	MKL_INT process_limit_B_supercell = ( coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) <= coulomb_basis.N_Gy_T + 2*supercell[0].RCL2*( basis.N_Gy_T / 2 ) ? coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) : coulomb_basis.N_Gy_T + 2*supercell[0].RCL2*( basis.N_Gy_T / 2 ) );
	MKL_INT process_num_B_supercell;

	// Sub-communicator splits communicator into 2 for process IDs under limit and above limit
	MPI_Comm MPI_SUB_COMM_B_supercell;

	// Local electronic supercell characteristic function
	struct characteristic_functions characteristic_function_coulomb_G_local;
	// Local B function and electronic structure supercell characteristic function convolutions
	struct convolution_char_state convolution_char_state_interim_local;
	// Local central portion (reduced) B function and electronic structure supercell characteristic function convolutions
	struct convolution_char_state convolution_char_state_interim_reduced_local;

	// The following store information about the local reduced convolution array in order to later gather it on a single process
	int *recvcounts_convolutions_B_supercell;
	int *displs_convolutions_B_supercell;

	if ( coulomb_region_type == 1 ) // Commensurate electronic structure supercell and Coulomb supercell lengths
	{

		if ( nprocs > process_limit_B_supercell ) // No. processes greater than process limit for Dfti FFT library
		{

			MPI_Comm_split( MPI_COMM_WORLD, ( myid < process_limit_B_supercell ? 0 : 1 ), myid, &MPI_SUB_COMM_B_supercell );
			process_num_B_supercell = process_limit_B_supercell;

		}
		else
		{

			MPI_Comm_dup( MPI_COMM_WORLD, &MPI_SUB_COMM_B_supercell );
			process_num_B_supercell = nprocs;

		}

		if ( myid < process_num_B_supercell ) // Processes with IDs less than the process limit
		{

			recvcounts_convolutions_B_supercell = ( int * ) malloc( process_num_B_supercell*sizeof( int ) );
			displs_convolutions_B_supercell = ( int * ) malloc( process_num_B_supercell*sizeof( int ) );

		}

	}

	// Total number of elements in each dimension for zero-padded arrays involved in Dfti FFT calculations as part of convolution calculation
	MKL_LONG sizes_B_supercell[] = { coulomb_basis.N_Gx_T + 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ), coulomb_basis.N_Gy_T + 2*supercell[0].RCL2*( basis.N_Gy_T / 2 ), coulomb_basis.N_Gz_T + 2*supercell[0].RCL3*( basis.N_Gz_T / 2 ) };

	// Global size of involved arrays
	MKL_LONG global_size_B_supercell = sizes_B_supercell[0]*sizes_B_supercell[1]*sizes_B_supercell[2];

	// Dfti precision & array type
	MKL_LONG DFTI_PRECISION_B_supercell = DFTI_DOUBLE;
	MKL_LONG DFTI_TYPE_B_supercell = DFTI_COMPLEX;

	// Dfti task descriptor
	DFTI_DESCRIPTOR_DM_HANDLE task_B_supercell;

	// Local memory size for arrays involved in Dfti FFT calculations
	MKL_LONG local_memory_size_B_supercell;
	
	// Global position of first element of local arrays involved in Dfti FFT calculations
	MKL_LONG local_start_B_supercell;
	
	// Number of global rows in local arrays involved in Dfti FFT calculations
	MKL_LONG local_rows_B_supercell;
	
	// Number of elements in local arrays involved in Dfti FFT calculations
	MKL_LONG local_size_B_supercell;


	// Variables relating to the extraction of the computationally relevant central portion of the calculated convolution
	MKL_INT local_start_reduced_B_supercell;
	MKL_INT local_rows_reduced_B_supercell;
	MKL_INT local_size_reduced_B_supercell;


	if ( coulomb_region_type == 1 ) // Commensurate electronic structure supercell and Coulomb supercell lengths
	{

		characteristic_function_coulomb_G.chi = ( complex double* ) calloc( coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T, sizeof( complex double ) );

		// Calculating Fourier transform of electronic structure supercell characteristic function in Coulomb supercell
		characteristic_function_coulomb_supercell( coulomb_basis, supercell, characteristic_function_coulomb_G );

		if ( myid < process_limit_B_supercell ) // Processes with IDs less than the process limit
		{

			// Initialising Dfti task descriptor
			check_dfti( DftiCreateDescriptorDM( MPI_SUB_COMM_B_supercell, &task_B_supercell, DFTI_PRECISION_B_supercell, DFTI_TYPE_B_supercell, 3, sizes_B_supercell ), "DftiCreateDescriptorDM" );

			// Setting Dfti backward scale to represent inverse FFT
			check_dfti( DftiSetValueDM( task_B_supercell, DFTI_BACKWARD_SCALE, 1.0 / global_size_B_supercell ), "DftiSetValueDM" );

			check_dfti( DftiGetValueDM( task_B_supercell, CDFT_LOCAL_SIZE, &local_memory_size_B_supercell ), "DftiGetValueDM" );
			
			check_dfti( DftiGetValueDM( task_B_supercell, CDFT_LOCAL_X_START, &local_start_B_supercell ), "DftiGetValueDM" );
			
			check_dfti( DftiGetValueDM( task_B_supercell, CDFT_LOCAL_NX, &local_rows_B_supercell ), "DftiGetValueDM" );

			local_size_B_supercell = local_rows_B_supercell*sizes_B_supercell[1]*sizes_B_supercell[2];

			// Committing Dfti descriptor
			check_dfti( DftiCommitDescriptorDM( task_B_supercell ), "DftiCommitDescriptorDM" );


			characteristic_function_coulomb_G_local.chi = ( complex double* ) calloc( local_memory_size_B_supercell, sizeof( complex double ) );
			convolution_states_interim_local.conv_state = ( complex double* ) calloc( local_memory_size_B_supercell, sizeof( complex double ) );
			convolution_char_state_interim_local.conv_char_state = ( complex double* ) calloc( local_memory_size_B_supercell, sizeof( complex double ) );

			// Extracting local parts of zero-padded electronic structure supercell characteristic function from global
			#pragma omp parallel for collapse( 3 )
			for ( int global_row = local_start_B_supercell; global_row < local_start_B_supercell + local_rows_B_supercell; global_row++ )
			{

				for ( int idx_y = 0; idx_y < sizes_B_supercell[1]; idx_y++ )
				{

					for ( int idx_z = 0; idx_z < sizes_B_supercell[2]; idx_z++ )
					{

						if ( global_row >= 0 && global_row < coulomb_basis.N_Gx_T && idx_y >= 0 && idx_y < coulomb_basis.N_Gy_T && idx_z >= 0 && idx_z < coulomb_basis.N_Gz_T )
						{

							characteristic_function_coulomb_G_local.chi[( global_row - local_start_B_supercell )*sizes_B_supercell[1]*sizes_B_supercell[2] + idx_y*sizes_B_supercell[2] + idx_z] = characteristic_function_coulomb_G.chi[ ( global_row )*( coulomb_basis.N_Gy_T )*( coulomb_basis.N_Gz_T ) + ( idx_y )*( coulomb_basis.N_Gz_T ) + ( idx_z ) ];

						}

					}

				}

			}

			// Compute FFT of zero-padded electronic structure supercell characteristic function
			check_dfti( DftiComputeForwardDM( task_B_supercell, characteristic_function_coulomb_G_local.chi ), "DftiComputeForwardDM" );

			// Setting local variables enabling extraction of central portion of convolved array, i.e. number of components of central portion that exist on each process and start point in local data
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

			// Gathering information regarding local number of elements in reduced array on each process that must be sent to all processes
			MPI_Allgatherv( &local_size_reduced_B_supercell, 1, MPI_INT, recvcounts_convolutions_B_supercell, recvcounts, displs, MPI_INT, MPI_SUB_COMM_B_supercell );

			// Calulating reduced array global displacements based on gathered data
			for ( int proc = 0; proc < process_num_B_supercell; proc++ )
			{

				displs_convolutions_B_supercell[proc] = ( proc == 0 ? 0 : displs_convolutions_B_supercell[proc - 1] + recvcounts_convolutions_B_supercell[proc - 1] );
			
			}

			convolution_char_state_interim_reduced_local.conv_char_state = ( complex double* ) calloc( local_size_reduced_B_supercell, sizeof( complex double ) );

		}

		free( characteristic_function_coulomb_G.chi );

	}

	
	// Calculating B functions for electronic structure supercell and upsampling them for Coulomb supercell
	
	if ( myid < process_limit ) // Processes with IDs less than the process limit
	{

		// Stores local component of convolution between envelope functions for each individual band, to later calculate B
		struct convolution_state convolution_states_local_bands;
		// Stores interim sum of local components of convolutions between envelope function band components, to later calculate B
		struct convolution_state convolution_states_local_interim;
		// Stores total sum of local components of convolutions between envelope function band components, to calculate B
		struct convolution_state convolution_states_local;
		// Stores global B for electronic structure supercell
		struct convolution_state convolution_states_interim;
		// Stores global B for Coulomb supercell
		struct convolution_state convolution_states_global;
		// Stores wave vector reversed version of global B for Coulomb supercell
		struct convolution_state convolution_states_global_reversed;

		convolution_states_local_bands.conv_state = calloc( N_BAND*local_memory_size, sizeof( complex double ) );
		convolution_states_local_interim.conv_state = calloc( local_size, sizeof( complex double ) );
		convolution_states_local.conv_state = calloc( local_size, sizeof( complex double ) );

		convolution_states_interim.conv_state = calloc( N_G_conv, sizeof( complex double ) );

		convolution_states_global.conv_state = calloc( N_GC_conv, sizeof( complex double ) );
		convolution_states_global_reversed.conv_state = calloc( N_GC_conv, sizeof( complex double ) );


		for ( int idx_state_1 = 0; idx_state_1 < 2*N_states[0] + 2*N_states[1]; idx_state_1++ ) // Loop over state 1 for B calculation
		{

			for ( int idx_state_2 = idx_state_1; idx_state_2 < 2*N_states[0] + 2*N_states[1]; idx_state_2++ ) // Loop over state 2 for B calculation
			{

				memset( convolution_states_local_interim.conv_state, 0, local_size*sizeof( complex double ) ); // Reset sum over bands to 0 for each set of states

				for ( int idx_band = 0; idx_band < N_BAND; idx_band++ ) // Loop over bands
				{

					if ( idx_state_1 < 2*N_states[0] && idx_state_2 < 2*N_states[0] ) // State 1 and 2 both VB states
					{

						// Element-by-element multiplication as part of convolution calculation via FFT method
						vzMul( local_size, &eigenstate_vb_padded_conj[idx_state_1].fourier[idx_band*local_memory_size], &eigenstate_vb_padded[idx_state_2].fourier[idx_band*local_memory_size], &convolution_states_local_bands.conv_state[idx_band*local_memory_size] );

					}
					else if ( idx_state_1 < 2*N_states[0] && idx_state_2 >= 2*N_states[0] ) // State 1 VB state and state 2 CB state
					{

						// Element-by-element multiplication as part of convolution calculation via FFT method
						vzMul( local_size, &eigenstate_vb_padded_conj[idx_state_1].fourier[idx_band*local_memory_size], &eigenstate_cb_padded[idx_state_2 - 2*N_states[0]].fourier[idx_band*local_memory_size], &convolution_states_local_bands.conv_state[idx_band*local_memory_size] );

					}
					else if ( idx_state_1 >= 2*N_states[0] && idx_state_2 < 2*N_states[0] ) // State 1 CB state and state 2 VB state
					{

						// Element-by-element multiplication as part of convolution calculation via FFT method
						vzMul( local_size, &eigenstate_cb_padded_conj[idx_state_1 - 2*N_states[0]].fourier[idx_band*local_memory_size], &eigenstate_vb_padded[idx_state_2].fourier[idx_band*local_memory_size], &convolution_states_local_bands.conv_state[idx_band*local_memory_size] );

					}
					else // State 1 and 2 both CB states
					{

						// Element-by-element multiplication as part of convolution calculation via FFT method
						vzMul( local_size, &eigenstate_cb_padded_conj[idx_state_1 - 2*N_states[0]].fourier[idx_band*local_memory_size], &eigenstate_cb_padded[idx_state_2 - 2*N_states[0]].fourier[idx_band*local_memory_size], &convolution_states_local_bands.conv_state[idx_band*local_memory_size] );

					}

					// Complete covolution calculation by calculating inverse FFT
					check_dfti( DftiComputeBackwardDM( task, &convolution_states_local_bands.conv_state[idx_band*local_memory_size] ), "DftiComputeBackwardDM" );

					// Adding contributions from each band
					vzAdd( local_size, &convolution_states_local_bands.conv_state[idx_band*local_memory_size], convolution_states_local_interim.conv_state, convolution_states_local.conv_state );

					// Storing interim sum if not final band
					if ( idx_band != N_BAND - 1 )
					{

						cblas_zcopy( local_size, convolution_states_local.conv_state, 1, convolution_states_local_interim.conv_state, 1 );

					}
					
				}


				// Gathering local parts of B to global
				MPI_Allgatherv( convolution_states_local.conv_state, local_size, MPI_DOUBLE_COMPLEX, convolution_states_interim.conv_state, recvcounts_convolutions, displs_convolutions, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
				
				
				// Upsampling B coefficients for Coulomb supercell
				
				if ( coulomb_region_type == 0 ) // If using free Coulomb supercell lengths
				{

					B_coulomb_supercell( basis, coulomb_basis, supercell, local_size_coulomb, convolution_states_interim, convolution_states_interim_local, recvcounts_coulomb, displs_coulomb, convolution_states_global );

				}
				else if ( coulomb_region_type == 1 ) // If using Coulomb supercell lengths that are commensurate with electronic structure supercell lengths
				{

					if ( myid < process_limit_B_supercell ) // If process ID is less than process limit
					{

						B_coulomb_supercell_fft( basis, coulomb_basis, supercell, process_limit_B_supercell, MPI_SUB_COMM_B_supercell, sizes_B_supercell, task_B_supercell, local_start_B_supercell, local_rows_B_supercell, local_size_B_supercell, local_memory_size_B_supercell, convolution_states_interim, convolution_states_interim_local, characteristic_function_coulomb_G_local, convolution_char_state_interim_local, local_start_reduced_B_supercell, local_rows_reduced_B_supercell, local_size_reduced_B_supercell, recvcounts_convolutions_B_supercell, displs_convolutions_B_supercell, convolution_char_state_interim_reduced_local, convolution_states_global );

					}

				}
				
				// Assigning wave vector reversed version of B 
				if ( myid == 0 )
				{

					for ( int idx_N = 0; idx_N < N_GC_conv; idx_N++ )
					{

						convolution_states_global_reversed.conv_state[idx_N] = convolution_states_global.conv_state[N_GC_conv - 1 - idx_N];

					}

				}

				// Distributing Coulomb supercell B among processes for later Vijkl calculation
				MPI_Scatterv( convolution_states_global.conv_state, sendcounts, displs_send, MPI_DOUBLE_COMPLEX, convolution_states[(idx_state_2 * (idx_state_2 + 1) / 2) + idx_state_1].conv_state, sendcounts[myid], MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD );
				MPI_Scatterv( convolution_states_global_reversed.conv_state, sendcounts, displs_send, MPI_DOUBLE_COMPLEX, convolution_states_reversed[(idx_state_2 * (idx_state_2 + 1) / 2) + idx_state_1].conv_state, sendcounts[myid], MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD );

			}

		}
		
		free( convolution_states_local_bands.conv_state );
		free( convolution_states_local_interim.conv_state );
		free( convolution_states_local.conv_state );
		free( convolution_states_interim.conv_state );
		free( convolution_states_global.conv_state );
		free( convolution_states_global_reversed.conv_state );

	}
	else // Processes with IDs greater than or equal to the process limit
	{

		// Stores global B for electronic structure supercell
		struct convolution_state convolution_states_interim;
		// Stores global B for Coulomb supercell
		struct convolution_state convolution_states_global;
		// Stores wave vector reversed version of global B for Coulomb supercell
		struct convolution_state convolution_states_global_reversed;

		convolution_states_interim.conv_state = calloc( N_G_conv, sizeof( complex double ) );
		convolution_states_global.conv_state = calloc( N_GC_conv, sizeof( complex double ) );

		// Stores total sum of local components of convolutions between envelope function band components, to calculate B
		struct convolution_state convolution_states_local;
		

		for ( int idx_state_1 = 0; idx_state_1 < 2*N_states[0] + 2*N_states[1]; idx_state_1++ ) // Loop over state 1 for B calculation
		{

			for ( int idx_state_2 = idx_state_1; idx_state_2 < 2*N_states[0] + 2*N_states[1]; idx_state_2++ ) // Loop over state 2 for B calculation
			{

				// Gathering local parts of B to global
				MPI_Allgatherv( convolution_states_local.conv_state, local_size, MPI_DOUBLE_COMPLEX, convolution_states_interim.conv_state, recvcounts_convolutions, displs_convolutions, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );

				
				// Upsampling B coefficients for Coulomb supercell
				
				if ( coulomb_region_type == 0 ) // If using free Coulomb supercell lengths
				{

					B_coulomb_supercell( basis, coulomb_basis, supercell, local_size_coulomb, convolution_states_interim, convolution_states_interim_local, recvcounts_coulomb, displs_coulomb, convolution_states_global );

				}
				else if ( coulomb_region_type == 1 ) // If using Coulomb supercell lengths that are commensurate with electronic structure supercell lengths
				{

					if ( myid < process_limit_B_supercell ) // If process ID is less than process limit
					{

						B_coulomb_supercell_fft( basis, coulomb_basis, supercell, process_limit_B_supercell, MPI_SUB_COMM_B_supercell, sizes_B_supercell, task_B_supercell, local_start_B_supercell, local_rows_B_supercell, local_size_B_supercell, local_memory_size_B_supercell, convolution_states_interim, convolution_states_interim_local, characteristic_function_coulomb_G_local, convolution_char_state_interim_local, local_start_reduced_B_supercell, local_rows_reduced_B_supercell, local_size_reduced_B_supercell, recvcounts_convolutions_B_supercell, displs_convolutions_B_supercell, convolution_char_state_interim_reduced_local, convolution_states_global );

					}

				}

				// Distributing Coulomb supercell B among processes for later Vijkl calculation
				MPI_Scatterv( convolution_states_global.conv_state, sendcounts, displs_send, MPI_DOUBLE_COMPLEX, convolution_states[(idx_state_2 * (idx_state_2 + 1) / 2) + idx_state_1].conv_state, sendcounts[myid], MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD );
				MPI_Scatterv( convolution_states_global_reversed.conv_state, sendcounts, displs_send, MPI_DOUBLE_COMPLEX, convolution_states_reversed[(idx_state_2 * (idx_state_2 + 1) / 2) + idx_state_1].conv_state, sendcounts[myid], MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD );

			}

		}

		free( convolution_states_interim.conv_state );
		free( convolution_states_global.conv_state );

	}

	free( recvcounts );
	free( displs );
	
	free( recvcounts_convolutions );
	free( displs_convolutions );

	if ( coulomb_region_type == 0 )
	{

		free( convolution_states_interim_local.conv_state );
		free( recvcounts_coulomb );
		free( displs_coulomb );

	}
	else if ( coulomb_region_type == 1 )
	{

		if ( myid < process_limit_B_supercell )
		{

			check_dfti( DftiFreeDescriptorDM( &task_B_supercell ), "DftiFreeDescriptorDM" );
			free( convolution_states_interim_local.conv_state );
			free( characteristic_function_coulomb_G_local.chi );
			free( convolution_char_state_interim_local.conv_char_state );
			free( convolution_char_state_interim_reduced_local.conv_char_state );

			free( recvcounts_convolutions_B_supercell );
			free( displs_convolutions_B_supercell );

		}

	}

}
