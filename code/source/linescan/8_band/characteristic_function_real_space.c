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

double characteristic_function_real_space( double x, double y, double z, struct wave_vectors basis, struct supercell_geometry *supercell, struct characteristic_functions characteristic_function_G )
{

	double result = 0.0;

	// Loop over basis set
	for( int idx_Gx = 0; idx_Gx < basis.N_Gx; idx_Gx++ )
	{

		for( int idx_Gy = 0; idx_Gy < basis.N_Gy; idx_Gy++ )
		{

			for( int idx_Gz = 0; idx_Gz < basis.N_Gz; idx_Gz++ )
			{

				// Wave vector index for stored reciprocal space characteristic function array
				int idx_basis = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_Gz + basis.N_Gz_T*(basis.N_Gy_T/2 - basis.N_Gy/2 + idx_Gy) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 - basis.N_Gx/2 + idx_Gx);

				// Wave vector indices for stored basis set
				int idx_basis_Gx = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_Gx;
				int idx_basis_Gy = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_Gy;
				int idx_basis_Gz = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_Gz;

				result += characteristic_function_G.chi[idx_basis]*cexp( I*( basis.Gx[idx_basis_Gx]*x + basis.Gy[idx_basis_Gy]*y + basis.Gz[idx_basis_Gz]*z ) );

			}

		}

	}


	return result;

}
