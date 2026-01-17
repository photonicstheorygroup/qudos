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

void B_coulomb_supercell_fft( struct wave_vectors basis, struct wave_vectors coulomb_basis, struct supercell_geometry *supercell, MKL_INT process_limit_B_supercell, MPI_Comm MPI_SUB_COMM_B_supercell, MKL_LONG *sizes_B_supercell, DFTI_DESCRIPTOR_DM_HANDLE task_B_supercell, MKL_LONG local_start_B_supercell, MKL_LONG local_rows_B_supercell, MKL_LONG local_size_B_supercell, MKL_LONG local_memory_size_B_supercell, struct convolution_state convolution_states_interim, struct convolution_state convolution_states_interim_local, struct characteristic_functions characteristic_function_coulomb_G_local, struct convolution_char_state convolution_char_state_interim_local, MKL_INT local_start_reduced_B_supercell, MKL_INT local_rows_reduced_B_supercell, MKL_INT local_size_reduced_B_supercell, int *recvcounts_convolutions_B_supercell, int *displs_convolutions_B_supercell, struct convolution_char_state convolution_char_state_interim_reduced_local, struct convolution_state convolution_states_global )
{
	
	// Processes with IDs less than the process limit
	if ( myid < process_limit_B_supercell )
	{

		// Setting intial sum to 0 for interim sum of local components of convolutions between B of electronic structure supercell and electronic structure characteristic function
		memset( convolution_states_interim_local.conv_state, 0, local_memory_size_B_supercell*sizeof( complex double ) );

		// Upsampling B for Coulomb supercell basis set
		#pragma omp parallel for collapse( 3 )
		for ( int global_row = local_start_B_supercell; global_row < local_start_B_supercell + local_rows_B_supercell; global_row++ ) // Looping over local wave vectors
		{

			for ( int idx_y = 0; idx_y < sizes_B_supercell[1]; idx_y++ )
			{

				for ( int idx_z = 0; idx_z < sizes_B_supercell[2]; idx_z++ )
				{

					if ( global_row >= 0 && global_row < 2*supercell[0].RCL1*( basis.N_Gx_T / 2 ) + 1 && idx_y >= 0 && idx_y <  2*supercell[0].RCL2*( basis.N_Gy_T / 2 ) + 1 && idx_z >= 0 && idx_z < 2*supercell[0].RCL3*( basis.N_Gz_T / 2 ) + 1 )
					{

						// Setting local components of electronic structure B array from global if Coulomb supercell wave vector is contained in electronic structure supercell basis set, otherwise set to 0
						if ( ( global_row - supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) % supercell[0].RCL1 == 0 && ( idx_y - supercell[0].RCL2*( basis.N_Gy_T / 2 ) ) % supercell[0].RCL2 == 0 && ( idx_z - supercell[0].RCL3*( basis.N_Gz_T / 2 ) ) % supercell[0].RCL3 == 0 )
						{

							convolution_states_interim_local.conv_state[( global_row - local_start_B_supercell )*sizes_B_supercell[1]*sizes_B_supercell[2] + idx_y*sizes_B_supercell[2] + idx_z] = convolution_states_interim.conv_state[ ( basis.N_Gx_T / 2 + ( global_row - supercell[0].RCL1*( basis.N_Gx_T / 2 ) ) / supercell[0].RCL1 )*( basis.N_Gy_T )*( basis.N_Gz_T ) + ( basis.N_Gy_T / 2 + ( idx_y - supercell[0].RCL2*( basis.N_Gy_T / 2 ) ) / supercell[0].RCL2 )*( basis.N_Gz_T ) + ( basis.N_Gz_T / 2 + ( idx_z - supercell[0].RCL3*( basis.N_Gz_T / 2 ) ) / supercell[0].RCL3 ) ];

						}

					}

				}

			}

		}

		// Compute FFT using Dfti library of upsampled electronic structure B
		check_dfti( DftiComputeForwardDM( task_B_supercell, convolution_states_interim_local.conv_state ), "DftiComputeForwardDM" );
		

		// Element by element multiplcation of FFT arrays
		vzMul( local_size_B_supercell, convolution_states_interim_local.conv_state, characteristic_function_coulomb_G_local.chi, convolution_char_state_interim_local.conv_char_state );


		// Compute inverse FFT using Dfti library of element by element multiplied arrays to complete convolution calculation
		check_dfti( DftiComputeBackwardDM( task_B_supercell, convolution_char_state_interim_local.conv_char_state ), "DftiComputeBackwardDM" );


		// Extracting central (reduced portion)
		#pragma omp parallel for collapse( 3 )
		for ( int idx_x = 0; idx_x < local_rows_reduced_B_supercell; idx_x++ )
		{

			for ( int idx_y = 0; idx_y < coulomb_basis.N_Gy_T; idx_y++ )
			{

				for ( int idx_z = 0; idx_z < coulomb_basis.N_Gz_T; idx_z++ )
				{

					convolution_char_state_interim_reduced_local.conv_char_state[idx_x*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T + idx_y*coulomb_basis.N_Gz_T + idx_z] = ( ( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3 ) / ( supercell[0].L1*supercell[0].L2*supercell[0].L3 ) )*convolution_char_state_interim_local.conv_char_state[ ( idx_x + local_start_reduced_B_supercell )*( sizes_B_supercell[1] )*( sizes_B_supercell[2] ) + ( idx_y + ( coulomb_basis.N_Gy_T + 2*supercell[0].RCL2*( basis.N_Gy_T / 2 ) ) / 2 - coulomb_basis.N_Gy_T / 2 )*( sizes_B_supercell[2] ) + ( idx_z + ( coulomb_basis.N_Gz_T + 2*supercell[0].RCL3*( basis.N_Gz_T / 2 ) ) / 2 - coulomb_basis.N_Gz_T / 2 ) ];

				}

			}

		}


		// Gathering B components calculated for Coulomb supercell
		MPI_Allgatherv( convolution_char_state_interim_reduced_local.conv_char_state, local_size_reduced_B_supercell, MPI_DOUBLE_COMPLEX, convolution_states_global.conv_state, recvcounts_convolutions_B_supercell, displs_convolutions_B_supercell, MPI_DOUBLE_COMPLEX, MPI_SUB_COMM_B_supercell );

	}

}
