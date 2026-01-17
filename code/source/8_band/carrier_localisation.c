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

complex double carrier_localisation( int idx_shape, struct wave_vectors basis, struct supercell_geometry *supercell, int *N_states, struct characteristic_functions *characteristic_function_G_local, int *shape_indicator, int state_type, struct eigenstate *eigenstate_stored, int idx_state )
{
	
	
	complex double result = 0;

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Calculating convolutions between state Fourier coefficients and shape characteristic function ////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Process limit for Dfti FFT library
	MKL_INT process_limit = ( basis.N_Gx + basis.N_Gx_T - 1 <= basis.N_Gy + basis.N_Gy_T - 1 ? basis.N_Gx + basis.N_Gx_T - 1 : basis.N_Gy + basis.N_Gy_T - 1 );

	// Sub-communicator splits communicator into 2 for process IDs under limit and above limit
	MPI_Comm MPI_SUB_COMM;
	if ( nprocs > process_limit ) // If number of MPI processes greater than process limit
	{

		MPI_Comm_split( MPI_COMM_WORLD, ( myid < process_limit ? 0 : 1 ), myid, &MPI_SUB_COMM ); // Split communicator

	}
	else // If number of MPI processes less than or equal to process limit
	{

		MPI_Comm_dup( MPI_COMM_WORLD, &MPI_SUB_COMM ); // Duplicate communicator

	}


	// Enveleope function plane wave basis set size
	int N_G = basis.N_Gx*basis.N_Gy*basis.N_Gz;
	
	if ( myid < process_limit ) // Processes with IDs less than the process limit
	{

		// Total number of elements in each dimension for zero-padded arrays involved in Dfti FFT calculations as part of convolution calculation
		MKL_LONG sizes[] = { basis.N_Gx + basis.N_Gx_T - 1, basis.N_Gy + basis.N_Gy_T - 1, basis.N_Gz + basis.N_Gz_T - 1 };

		// Global size of involved arrays
		MKL_LONG global_size = sizes[0]*sizes[1]*sizes[2];

		// Dfti precision & array type
		MKL_LONG DFTI_PRECISION = DFTI_DOUBLE;
		MKL_LONG DFTI_TYPE = DFTI_COMPLEX;

		// Dfti task descriptor
		DFTI_DESCRIPTOR_DM_HANDLE task;

		// Initialising Dfti task descriptor
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

		// Committing Dfti descriptor
		check_dfti( DftiCommitDescriptorDM( task ), "DftiCommitDescriptorDM" );


		// Calculate FFT of shape characteristic function if not already done so
		if ( shape_indicator[idx_shape] == 0 )
		{

			// Will hold local part of shape characteristic function
			characteristic_function_G_local[idx_shape].chi = ( complex double* ) calloc( local_memory_size, sizeof( complex double ) );
			shape_indicator[idx_shape] = 1;

			// Extracting local parts of zero-padded shape characteristic functions from global
			#pragma omp parallel for collapse( 3 )
			for ( int global_row = local_start; global_row < local_start + local_rows; global_row++ )
			{

				for ( int idx_y = 0; idx_y < sizes[1]; idx_y++ )
				{

					for ( int idx_z = 0; idx_z < sizes[2]; idx_z++ )
					{

						if ( global_row >= 0 && global_row < basis.N_Gx_T && idx_y >= 0 && idx_y < basis.N_Gy_T && idx_z >= 0 && idx_z < basis.N_Gz_T )
						{

							characteristic_function_G_local[idx_shape].chi[( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = characteristic_function_G[idx_shape].chi[ ( global_row )*( basis.N_Gy_T )*( basis.N_Gz_T ) + ( idx_y )*( basis.N_Gz_T ) + ( idx_z ) ];

						}

					}

				}

			}

			// Compute FFT of shape characteristic function using DFTI library
			check_dfti( DftiComputeForwardDM( task, characteristic_function_G_local[idx_shape].chi ), "DftiComputeForwardDM" );

		}


		// Will hold local parts of padded envelope functions involved in convolution calculation
		struct eigenstate eigenstate_padded;
		
		eigenstate_padded.fourier = calloc( N_BAND*local_memory_size, sizeof( complex double ) );

		// Interim eigenstate to facilitate state redistribution among processes
		struct eigenstate eigenstate_interim;

		eigenstate_interim.fourier = calloc( N_H     , sizeof( complex double ) );

		// Combining the block cyclically distributed eigenstate from all processes and sending a copy to all processes
		redistribute_eigenstate( basis, N_states, state_type, idx_state, eigenstate_stored, eigenstate_interim );

		// Will store computed convolution
		complex double *convolution_results =  calloc( N_BAND*local_memory_size, sizeof( complex double ) );

		
		for ( int idx_band = 0; idx_band < N_BAND; idx_band++ ) // Iterate over bands
		{

			// Setting the local parts of the padded envelope functions to be used in the convolution calculation
			#pragma omp parallel for collapse( 3 )
			for ( int global_row = local_start; global_row < local_start + local_rows; global_row++ )
			{

				for ( int idx_y = 0; idx_y < sizes[1]; idx_y++ )
				{

					for ( int idx_z = 0; idx_z < sizes[2]; idx_z++ )
					{

						if ( global_row >= 0 && global_row < basis.N_Gx && idx_y >= 0 && idx_y < basis.N_Gy && idx_z >= 0 && idx_z < basis.N_Gz )
						{

							eigenstate_padded.fourier[idx_band*local_memory_size + ( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = eigenstate_interim.fourier[idx_band*N_G + ( global_row )*( basis.N_Gy )*( basis.N_Gz ) + ( idx_y )*( basis.N_Gz ) + ( idx_z )];

						}

					}

				}

			}

			// Compute FFT of band envelope function using DFT library
			check_dfti( DftiComputeForwardDM( task, &eigenstate_padded.fourier[idx_band*local_memory_size] ), "DftiComputeForwardDM" );

			// Element-by-element multiplication as part of convolution calculation via FFT method
			vzMul( local_size, &eigenstate_padded.fourier[idx_band*local_memory_size], characteristic_function_G_local[idx_shape].chi, &convolution_results[idx_band*local_memory_size] );

			// Complete covolution calculation by calculating inverse FFT
			check_dfti( DftiComputeBackwardDM( task, &convolution_results[idx_band*local_memory_size] ), "DftiComputeBackwardDM" );


			#pragma omp parallel for reduction(+:result) schedule(guided)
			for ( int global_row = local_start; global_row < local_start + local_rows; global_row++ )
			{
				for ( int idx_y = 0; idx_y < sizes[1]; idx_y++ )
				{

					for ( int idx_z = 0; idx_z < sizes[2]; idx_z++ )
					{

						if ( global_row >= ( basis.N_Gx + basis.N_Gx_T - 1 )/2 - basis.N_Gx/2 && global_row < ( basis.N_Gx + basis.N_Gx_T - 1 )/2 + basis.N_Gx/2 + 1 && idx_y >= ( basis.N_Gy + basis.N_Gy_T - 1 )/2 - basis.N_Gy/2 && idx_y < ( basis.N_Gy + basis.N_Gy_T - 1 )/2 + basis.N_Gy/2 + 1 && idx_z >= ( basis.N_Gz + basis.N_Gz_T - 1 )/2 - basis.N_Gz/2 && idx_z < ( basis.N_Gz + basis.N_Gz_T - 1 )/2 + basis.N_Gz/2 + 1 )
						{

							result += convolution_results[idx_band*local_memory_size + ( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z]*conj( eigenstate_interim.fourier[idx_band*N_G + ( global_row - ( basis.N_Gx + basis.N_Gx_T - 1 )/2 + basis.N_Gx/2 )*( basis.N_Gy )*( basis.N_Gz ) + ( idx_y - ( basis.N_Gy + basis.N_Gy_T - 1 )/2 + basis.N_Gy/2 )*( basis.N_Gz ) + ( idx_z - ( basis.N_Gz + basis.N_Gz_T - 1 )/2 + basis.N_Gz/2 )] );

						}

					}

				}

			}

		}

		free( eigenstate_interim.fourier );
		
		free( convolution_results );
		free( eigenstate_padded.fourier );

	}
	else
	{

		// Interim eigenstate to facilitate redistribution
		struct eigenstate eigenstate_interim;

		eigenstate_interim.fourier   = calloc( N_H     , sizeof( complex double ) );

		// Combining the block cyclically distributed eigenstate from all processes and sending a copy to all processes
		redistribute_eigenstate( basis, N_states, state_type, idx_state, eigenstate_stored, eigenstate_interim );

		free( eigenstate_interim.fourier );

	}

	return result;

} 
