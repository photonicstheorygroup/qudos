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

complex double field_fourier_coefficient( struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis )
{

	// Will store result
	complex double result = 0.0 + I*0.0;


	// Compute Fourier transform of electrostatic potential
	if( idx_G_1[0] != idx_G_2[0] ) // x component
	{

		result += ( I*F.x*pow( -1.0, ( double ) idx_G_1[0] - idx_G_2[0] + 1 )*kronecker_delta( idx_G_1[1], idx_G_2[1] )*kronecker_delta( idx_G_1[2], idx_G_2[2] ) )/( basis.Gx[idx_G_2[0]] - basis.Gx[idx_G_1[0]] );

	}

	if( idx_G_1[1] != idx_G_2[1] ) // y component
	{

		result += ( I*F.y*pow( -1.0, ( double ) idx_G_1[1] - idx_G_2[1] + 1 )*kronecker_delta( idx_G_1[0], idx_G_2[0] )*kronecker_delta( idx_G_1[2], idx_G_2[2] ) )/( basis.Gy[idx_G_2[1]] - basis.Gy[idx_G_1[1]] );

	}

	if( idx_G_1[2] != idx_G_2[2] ) // z component
	{

		result += ( I*F.z*pow( -1.0, ( double ) idx_G_1[2] - idx_G_2[2] + 1 )*kronecker_delta( idx_G_1[0], idx_G_2[0] )*kronecker_delta( idx_G_1[1], idx_G_2[1] ) )/( basis.Gz[idx_G_2[2]] - basis.Gz[idx_G_1[2]] );

	}


	return result;

} // End of function
