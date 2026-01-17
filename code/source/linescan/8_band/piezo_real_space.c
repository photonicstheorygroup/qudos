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

void piezo_real_space( int N_shape, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, complex double *piezo_G, complex double *piezo_r )
{

	int idx_grid, idx_Gx, idx_Gy, idx_Gz, idx_basis;


	// Loop over grid points

	#pragma omp parallel for private( idx_Gx, idx_Gy, idx_Gz, idx_basis )
	for( idx_grid = 0; idx_grid < N_grid; idx_grid++ )
	{

		piezo_r[idx_grid] = 0.0 + I*0.0;


		// Compute piezoelectric potential at current real space grid point

		for( idx_Gx = 0; idx_Gx < basis.N_Gx; idx_Gx++ )
		{

			for( idx_Gy = 0; idx_Gy < basis.N_Gy; idx_Gy++ )
			{

				for( idx_Gz = 0; idx_Gz < basis.N_Gz; idx_Gz++ )
				{

					// Set Fourier coefficient (plane wave) index

					idx_basis = idx_Gz + basis.N_Gz*idx_Gy + basis.N_Gy*basis.N_Gz*idx_Gx;


					int idx_basis_Gx = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_Gx;
					int idx_basis_Gy = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_Gy;
					int idx_basis_Gz = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_Gz;


					// Add Fourier component for the specified wave vector

					piezo_r[idx_grid] += piezo_G[idx_basis]*cexp( I*( basis.Gx[idx_basis_Gx]*x[idx_grid] + basis.Gy[idx_basis_Gy]*y[idx_grid] + basis.Gz[idx_basis_Gz]*z[idx_grid] ) );

				}

			}

		}

	}

}
