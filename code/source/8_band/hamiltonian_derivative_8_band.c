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

void hamiltonian_derivative_8_band( int polarisation, int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, complex double* dH_dk )
{

	#pragma omp parallel for collapse(2) schedule(guided)
	for( int local_row_index = 0; local_row_index < na_rows; local_row_index++ ) // Looping over local matrix rows
	{
		
		for( int local_col_index = 0; local_col_index < na_cols; local_col_index++ ) // Looping over local matrix columns
		{
		
			int global_row_index = indxl2g(local_row_index, block_size, my_prow, 0, np_rows); // Global row index
			int global_col_index = indxl2g(local_col_index, block_size, my_pcol, 0, np_cols); // Global column index
			
			if (global_row_index >= global_col_index) // Only setting values of lower triangular portion as Hamiltonian matrix is Hermitian
			{
				
				int idx_row_block = global_row_index / (N_H / 8) + 1; // Identifying band in 8-band row
				int idx_col_block = global_col_index / (N_H / 8) + 1; // Identifying band in 8-band column
				
				int idx_G_1[3], idx_G_2[3]; // Will store global wave vector indices
	
				// Identifying global row wave vector indices
				idx_G_1[0] = (global_row_index % (N_H / 8)) / (basis.N_Gy*basis.N_Gz);
				idx_G_1[1] = ((global_row_index % (N_H / 8)) % (basis.N_Gy*basis.N_Gz)) / basis.N_Gz;
				idx_G_1[2] = ((global_row_index % (N_H / 8)) % (basis.N_Gy*basis.N_Gz)) % basis.N_Gz;
				
				// Identifying global column wave vector indices
				idx_G_2[0] = (global_col_index % (N_H / 8)) / (basis.N_Gy*basis.N_Gz);
				idx_G_2[1] = ((global_col_index % (N_H / 8)) % (basis.N_Gy*basis.N_Gz)) / basis.N_Gz;
				idx_G_2[2] = ((global_col_index % (N_H / 8)) % (basis.N_Gy*basis.N_Gz)) % basis.N_Gz;
				
				MKL_INT idx_row = local_row_index;
				MKL_INT idx_col = local_col_index;
				
				
				MKL_INT index; // Local matrix index
				
				if (PACKED && myid == 0) // If using packed matrix storage
				{
					
					if ( idx_row >= idx_col ) // Only setting values of lower triangular portion as Hamiltonian matrix is Hermitian
					{
					
						index = ( idx_col * ( 2*N_H - idx_col + 1 ) ) / 2 + idx_row - idx_col;
					
					}
				
				}
				else // Otherwise
				{
				
					index = idx_row + na_rows*idx_col;
				
				}

				switch( idx_col_block ) // Setting 8-band matrix element depending on column band
				{

					// Column 1

					case 1:

						switch( idx_row_block ) // Setting 8-band matrix element depending on row band
						{

							// Column 1, row 1

							case 1:
								
								dH_dk[index] = dE_CB_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 1, row 2

							case 2:

								dH_dk[index] = -1.0*sqrt( 3.0 )*dT_minus_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 1, row 3

							case 3:

								dH_dk[index] = sqrt( 2.0 )*dU_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 1, row 4

							case 4:

								dH_dk[index] = -1.0*dU_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 1, row 5

							case 5:

								dH_dk[index] = 0.0 + I*0.0;

								break;


							// Column 1, row 6

							case 6:

								dH_dk[index] = 0.0 + I*0.0;
								
								break;


							// Column 1, row 7

							case 7:

								dH_dk[index] = -1.0*dT_plus_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 1, row 8

							case 8:

								dH_dk[index] = -1.0*sqrt( 2.0 )*dT_plus_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;

						}

						break;


					// Column 2

					case 2:

						switch( idx_row_block )
						{

							// Column 2, row 2

							case 2:

								dH_dk[index] = dE_HH_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 2, row 3

							case 3:

								dH_dk[index] = sqrt( 2.0 )*dS_conj_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 2, row 4

							case 4:

								dH_dk[index] = -1.0*dS_conj_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 2, row 5

							case 5:

								dH_dk[index] = 0.0 + I*0.0;

								break;


							// Column 2, row 6

							case 6:

								dH_dk[index] = 0.0 + I*0.0;

								break;


							// Column 2, row 7

							case 7:

								dH_dk[index] = -1.0*dR_conj_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 2, row 8

							case 8:

								dH_dk[index] = -1.0*sqrt( 2.0 )*dR_conj_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;

						}

						break;


					// Column 3

					case 3:

						switch( idx_row_block )
						{

							// Column 3, row 3

							case 3:

								dH_dk[index] = dE_LH_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 3, row 4

							case 4:

								dH_dk[index] = dQ_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 3, row 5

							case 5:

								dH_dk[index] = dT_plus_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 3, row 6

							case 6:

								dH_dk[index] = dR_conj_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 3, row 7

							case 7:

								dH_dk[index] = 0.0 + I*0.0;

								break;


							// Column 3, row 8

							case 8:

								dH_dk[index] = sqrt( 3.0 )*dS_conj_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;

						}

						break;


					// Column 4

					case 4:

						switch( idx_row_block )
						{

							// Column 4, row 4

							case 4:

								dH_dk[index] = dE_SO_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 4, row 5

							case 5:

								dH_dk[index] = sqrt( 2.0 )*dT_plus_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 4, row 6

							case 6:

								dH_dk[index] = sqrt( 2.0 )*dR_conj_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 4, row 7

							case 7:

								dH_dk[index] = -1.0*sqrt( 3.0 )*dS_conj_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );
								
								break;


							// Column 4, row 8

							case 8:

								dH_dk[index] = 0.0 + I*0.0;

								break;

						}

						break;


					// Column 5

					case 5:

						switch( idx_row_block )
						{

							// Column 5, row 5

							case 5:

								dH_dk[index] = dE_CB_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 5, row 6

							case 6:

								dH_dk[index] = -1.0*sqrt( 3.0 )*dT_plus_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 5, row 7

							case 7:

								dH_dk[index] = sqrt( 2.0 )*dU_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 5, row 8

							case 8:

								dH_dk[index] = -1.0*dU_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;

						}

						break;


					// Column 6

					case 6:

						switch( idx_row_block )
						{

							// Column 6, row 6

							case 6:

								dH_dk[index] = dE_HH_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 6, row 7

							case 7:

								dH_dk[index] = sqrt( 2.0 )*dS_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );
								
								break;


							// Column 6, row 8

							case 8:

								dH_dk[index] = -1.0*dS_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;

						}

						break;


					// Column 7

					case 7:

						switch( idx_row_block )
						{

							// Column 7, row 7

							case 7:

								dH_dk[index] = dE_LH_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;


							// Column 7, row 8

							case 8:

								dH_dk[index] = 1.0*dQ_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;

						}

						break;


					// Column 8

					case 8:

						switch( idx_row_block )
						{

							// Column 8, row 8

							case 8:

								dH_dk[index] = dE_SO_dk( polarisation, N_shape, idx_G_1, idx_G_2, basis, supercell, parameters );

								break;

						}

						break;

				}
				
			}
			
		}

	}

}
