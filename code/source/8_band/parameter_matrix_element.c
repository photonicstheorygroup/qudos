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

complex double parameter_matrix_element( int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, double *parameter )
{

	complex double result;
	
	// Wave vector index
	int idx_basis = (idx_G_1[2] + basis.N_Gz_T/2 - idx_G_2[2]) + basis.N_Gz_T*(idx_G_1[1] + basis.N_Gy_T/2 - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(idx_G_1[0] + basis.N_Gx_T/2 - idx_G_2[0]);


	// Matrix term
	result = parameter[0]*kronecker_delta( idx_G_1[0], idx_G_2[0] )*kronecker_delta( idx_G_1[1], idx_G_2[1] )*kronecker_delta( idx_G_1[2], idx_G_2[2] );


	// Shape terms
	int idx_shape;

	for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
	{

		result += ( parameter[idx_shape] - parameter[supercell[idx_shape].embedded] )*characteristic_function_G[idx_shape].chi[idx_basis];

	}

	return result;

} // End of function
