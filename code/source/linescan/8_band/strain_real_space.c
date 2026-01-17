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

void strain_real_space( int N_shape, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct strain_tensor strain_G, struct strain_tensor strain_r )
{

	// Loop over grid points
	#pragma omp parallel for
	for( int idx_grid = 0; idx_grid < N_grid; idx_grid++ )
	{
		
		strain_r.xx[idx_grid] = 0.0 + I*0.0;
		strain_r.yy[idx_grid] = 0.0 + I*0.0;
		strain_r.zz[idx_grid] = 0.0 + I*0.0;
		strain_r.yz[idx_grid] = 0.0 + I*0.0;
		strain_r.xz[idx_grid] = 0.0 + I*0.0;
		strain_r.xy[idx_grid] = 0.0 + I*0.0;
		
		// Compute strain at current real space grid point
		
		for( int idx_Gx = 0; idx_Gx < basis.N_Gx; idx_Gx++ )
		{
			
			for( int idx_Gy = 0; idx_Gy < basis.N_Gy; idx_Gy++ )
			{

				for( int idx_Gz = 0; idx_Gz < basis.N_Gz; idx_Gz++ )
				{
					
					// Set Fourier coefficient (plane wave) index
					
					int idx_basis = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_Gz + basis.N_Gz_T*(basis.N_Gy_T/2 - basis.N_Gy/2 + idx_Gy) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 - basis.N_Gx/2 + idx_Gx);

					int idx_basis_Gx = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_Gx;
					int idx_basis_Gy = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_Gy;
					int idx_basis_Gz = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_Gz;
					
					// Calculate exponential part of strain element contribution
					
					complex double exp_part = cexp( I*( basis.Gx[idx_basis_Gx]*x[idx_grid] + basis.Gy[idx_basis_Gy]*y[idx_grid] + basis.Gz[idx_basis_Gz]*z[idx_grid] ) );

					//  Add contribution to each strain tensor component 
					
					strain_r.xx[idx_grid] += strain_G.xx[idx_basis]*exp_part;

					strain_r.yy[idx_grid] += strain_G.yy[idx_basis]*exp_part;

					strain_r.zz[idx_grid] += strain_G.zz[idx_basis]*exp_part;

					strain_r.yz[idx_grid] += strain_G.yz[idx_basis]*exp_part;

					strain_r.xz[idx_grid] += strain_G.xz[idx_basis]*exp_part;

					strain_r.xy[idx_grid] += strain_G.xy[idx_basis]*exp_part;

				}

			}

		}

	}

}
