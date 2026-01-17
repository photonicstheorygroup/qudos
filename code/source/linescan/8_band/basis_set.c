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

void basis_set( struct wave_vectors basis, struct supercell_geometry *geometry, int strain_region )
{

	// Note: geometry[0].L1 = Length of supercell along x
	//       geometry[0].L2 = Length of supercell along y
	//       geometry[0].L3 = Length of supercell along z

	int idx_G;
	
	double L1, L2, L3;

	if ( strain_region == 0 )	// Electronic structure supercell lengths
	{
	
		L1 = geometry[0].L1;
		L2 = geometry[0].L2;
		L3 = geometry[0].L3;
	
	}
	
	else	// Strain supercell lengths
	{
	
		L1 = geometry[0].SL1;
		L2 = geometry[0].SL2;
		L3 = geometry[0].SL3;
	
	}


	// Wave vectors along x direction

	for( idx_G = 0; idx_G < basis.N_Gx_T; idx_G++ )
	{

		basis.Gx[idx_G] = ( ( 2.0*M_PI )/L1 )*( 0.5*( 1.0 - ( double ) basis.N_Gx_T ) + ( double ) idx_G );

	}


	// Wave vectors along y direction

	for( idx_G = 0; idx_G < basis.N_Gy_T; idx_G++ )
	{

		basis.Gy[idx_G] = ( ( 2.0*M_PI )/L2 )*( 0.5*( 1.0 - ( double ) basis.N_Gy_T ) + ( double ) idx_G );

	}


	// Wave vectors along z direction

	for( idx_G = 0; idx_G < basis.N_Gz_T; idx_G++ )
	{

		basis.Gz[idx_G] = ( ( 2.0*M_PI )/L3 )*( 0.5*( 1.0 - ( double ) basis.N_Gz_T ) + ( double ) idx_G );

	}


}
