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

void state_redistribution( struct wave_vectors basis, struct supercell_geometry *supercell, int *N_states, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb, struct eigenstate *eigenstate_vb_padded, struct eigenstate *eigenstate_cb_padded, struct eigenstate *eigenstate_vb_padded_conj, struct eigenstate *eigenstate_cb_padded_conj, MKL_INT process_limit, MKL_LONG *sizes, MKL_LONG local_memory_size, MKL_LONG local_start, MKL_LONG local_rows )
{

	// Enveleope function plane wave basis set size
	int N_G = basis.N_Gx*basis.N_Gy*basis.N_Gz;
	
	if ( myid < process_limit ) // Processes with IDs less than the process limit
	{

		// Interim eigenstate to facilitate redistribution
		struct eigenstate eigenstate_interim;

		eigenstate_interim.fourier   = calloc( N_H     , sizeof( complex double ) );

		for ( int idx_state_vb = 0; idx_state_vb < 2*N_states[0]; idx_state_vb++ ) // Iterating over VB states
		{

			memset( eigenstate_interim.fourier, 0, N_H*sizeof( complex double ) ); // Intiallly setting Fourier coefficients of interim eigenstate to 0

			// Combining the block cyclically distributed eigenstate from all processes and sending a copy to all processes
			redistribute_eigenstate( basis, N_states, 0, idx_state_vb, eigenstate_vb, eigenstate_interim );

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

								eigenstate_vb_padded[idx_state_vb].fourier[idx_band*local_memory_size + ( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = eigenstate_interim.fourier[idx_band*N_G + ( global_row )*( basis.N_Gy )*( basis.N_Gz ) + ( idx_y )*( basis.N_Gz ) + ( idx_z )];
								eigenstate_vb_padded_conj[idx_state_vb].fourier[idx_band*local_memory_size + ( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = conj( eigenstate_interim.fourier[idx_band*N_G + ( basis.N_Gx - 1 - global_row )*( basis.N_Gy )*( basis.N_Gz ) + ( basis.N_Gy - 1 - idx_y )*( basis.N_Gz ) + ( basis.N_Gz - 1 - idx_z )] );

							}
							else
							{

								eigenstate_vb_padded[idx_state_vb].fourier[idx_band*local_memory_size + ( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = 0;
								eigenstate_vb_padded_conj[idx_state_vb].fourier[idx_band*local_memory_size + ( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = 0;

							}

						}

					}

				}

			}

		}

		for ( int idx_state_cb = 0; idx_state_cb < 2*N_states[1]; idx_state_cb++ ) // Iterating over CB states
		{

			memset( eigenstate_interim.fourier, 0, N_H*sizeof( complex double ) ); // Initiallly setting Fourier coefficients of interim eigenstate to 0

			// Combining the block cyclically distributed eigenstate from all processes and sending a copy to all processes
			redistribute_eigenstate( basis, N_states, 1, idx_state_cb, eigenstate_cb, eigenstate_interim ); 

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

							if ( global_row >= 0 && global_row < basis.N_Gx && idx_y >= 0 && idx_y <  basis.N_Gy && idx_z >= 0 && idx_z < basis.N_Gz )
							{

								eigenstate_cb_padded[idx_state_cb].fourier[idx_band*local_memory_size + ( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = eigenstate_interim.fourier[idx_band*N_G + ( global_row )*( basis.N_Gy )*( basis.N_Gz ) + ( idx_y )*( basis.N_Gz ) + ( idx_z )];
								eigenstate_cb_padded_conj[idx_state_cb].fourier[idx_band*local_memory_size + ( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = conj( eigenstate_interim.fourier[idx_band*N_G + ( basis.N_Gx - 1 - global_row )*( basis.N_Gy )*( basis.N_Gz ) + ( basis.N_Gy - 1 - idx_y )*( basis.N_Gz ) + ( basis.N_Gz - 1 - idx_z )] );

							}
							else
							{

								eigenstate_cb_padded[idx_state_cb].fourier[idx_band*local_memory_size + ( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = 0;
								eigenstate_cb_padded_conj[idx_state_cb].fourier[idx_band*local_memory_size + ( global_row - local_start )*sizes[1]*sizes[2] + idx_y*sizes[2] + idx_z] = 0;

							}

						}

					}

				}

			}

		}

		free( eigenstate_interim.fourier );

	}
	else // Processes with IDs greater than or equal to the process limit
	{

		struct eigenstate eigenstate_interim;

		eigenstate_interim.fourier   = calloc( N_H     , sizeof( complex double ) );

		// Iterating over VB states
		for ( int idx_state_vb = 0; idx_state_vb < 2*N_states[0]; idx_state_vb++ )
		{

			// Initially setting Fourier coefficients of interim eigenstate to 0
			memset( eigenstate_interim.fourier, 0, N_H*sizeof( complex double ) );

			// Combining the block cyclically distributed eigenstate from all processes and sending a copy to all processes
			redistribute_eigenstate( basis, N_states, 0, idx_state_vb, eigenstate_vb, eigenstate_interim );

		}
		
		// Iterating over CB states
		for ( int idx_state_cb = 0; idx_state_cb < 2*N_states[1]; idx_state_cb++ )
		{

			// Initially setting Fourier coefficients of interim eigenstate to 0
			memset( eigenstate_interim.fourier, 0, N_H*sizeof( complex double ) );

			// Combining the block cyclically distributed eigenstate from all processes and sending a copy to all processes
			redistribute_eigenstate( basis, N_states, 1, idx_state_cb, eigenstate_cb, eigenstate_interim );

		}

		free( eigenstate_interim.fourier );

	}

}
