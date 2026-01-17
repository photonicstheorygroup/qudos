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


#include "8_band_linescan.h"

void characteristic_functions_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct characteristic_functions *characteristic_function_G, int strain_region )
{

	int shape_index_start;
	if ( strain_region == 0 || strain_region == 1 ) // Electronic structure supercell or strain supercell with free lengths indicator
	{

		shape_index_start = 1;

	}
	else // Otherwise strain supercell with lengths commensurate with electronic structure supercell lengths
	{

		shape_index_start = 0;

	}

	#pragma omp parallel
	{

		#pragma omp for collapse( 4 )
		for( int idx_shape = shape_index_start; idx_shape < N_shape; idx_shape++ )
		{

			// Loop over wave vectors and evaluate the Fourier transforms of the shape characteristic functions
			for( int idx_Gx = 0; idx_Gx < basis.N_Gx_T; idx_Gx++ )
			{

				for( int idx_Gy = 0; idx_Gy < basis.N_Gy_T; idx_Gy++ )
				{

					for( int idx_Gz = 0; idx_Gz < basis.N_Gz_T; idx_Gz++ )
					{

						int idx_basis = idx_Gz + basis.N_Gz_T*idx_Gy + basis.N_Gy_T*basis.N_Gz_T*idx_Gx;

						characteristic_function_G[idx_shape].chi[idx_basis] = characteristic_function_fourier( basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, idx_shape, strain_region );

					}

				}

			}

		}

	}

}
