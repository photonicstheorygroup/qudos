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

complex double piezo_fourier_coefficient( int piezo_order, int idx_shape, int idx_Gx_global, int idx_Gy, int idx_Gz, int idx_Gx_local, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters,  struct strain_tensor strain_region_G, struct convolution_char_strain convolutions_char_strain_region_reduced_local )
{
	
	double Gx = basis.Gx[idx_Gx_global];
	double Gy = basis.Gy[idx_Gy];
	double Gz = basis.Gz[idx_Gz];
	
	double G_mag = sqrt( pow( Gx, 2.0 ) + pow( Gy, 2.0 ) + pow( Gz, 2.0 ) ); // |G|
	
	complex double result = 0.0 + I*0.0;


	// Compute Fourier transform of piezoelectric polarisation

	struct vector_complex P_G;

	P_G.x = 0.0 + I*0.0;
	P_G.y = 0.0 + I*0.0;
	P_G.z = 0.0 + I*0.0;
	
	// Wave vector index for global arrays
	int idx_basis_global = idx_Gx_global*basis.N_Gy_T*basis.N_Gz_T + idx_Gy*basis.N_Gz_T + idx_Gz;

	// Wave vector index for local arrays
	int idx_basis_local = idx_Gx_local*basis.N_Gy_T*basis.N_Gz_T + idx_Gy*basis.N_Gz_T + idx_Gz;


	if( G_mag > EPS ) // |G| != 0
	{

		switch( piezo_order )
		{

			case 1: // First order

				if ( idx_shape == 0) // Matrix material contribution
				{

					P_G.x += 2.0*parameters[0].e14*strain_region_G.yz[idx_basis_global];
					P_G.y += 2.0*parameters[0].e14*strain_region_G.xz[idx_basis_global];
					P_G.z += 2.0*parameters[0].e14*strain_region_G.xy[idx_basis_global];

				}
				else // Shape contribution
				{

					P_G.x += 2.0*( parameters[idx_shape].e14 - parameters[supercell[idx_shape].embedded].e14 )*( convolutions_char_strain_region_reduced_local.char_yz[idx_basis_local] );
					P_G.y += 2.0*( parameters[idx_shape].e14 - parameters[supercell[idx_shape].embedded].e14 )*( convolutions_char_strain_region_reduced_local.char_xz[idx_basis_local] );
					P_G.z += 2.0*( parameters[idx_shape].e14 - parameters[supercell[idx_shape].embedded].e14 )*( convolutions_char_strain_region_reduced_local.char_xy[idx_basis_local] );

				}

				break;

			case 2: // Second order
			
				if ( idx_shape == 0) // Matrix material contribution
				{

					P_G.x += 2.0*parameters[0].B114*convolutions_char_strain_region_reduced_local.xx_yz[idx_basis_local]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region_reduced_local.yy_yz[idx_basis_local]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region_reduced_local.zz_yz[idx_basis_local]
									+ 4.0*parameters[0].B156*convolutions_char_strain_region_reduced_local.xy_xz[idx_basis_local];

					P_G.y += 2.0*parameters[0].B114*convolutions_char_strain_region_reduced_local.yy_xz[idx_basis_local]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region_reduced_local.xx_xz[idx_basis_local]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region_reduced_local.zz_xz[idx_basis_local]
									+ 4.0*parameters[0].B156*convolutions_char_strain_region_reduced_local.xy_yz[idx_basis_local];

					P_G.z += 2.0*parameters[0].B114*convolutions_char_strain_region_reduced_local.zz_xy[idx_basis_local]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region_reduced_local.xx_xy[idx_basis_local]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region_reduced_local.yy_xy[idx_basis_local]
									+ 4.0*parameters[0].B156*convolutions_char_strain_region_reduced_local.xz_yz[idx_basis_local];

				}
				else // Shape contribution
				{

					P_G.x += 2.0*( parameters[idx_shape].B114 - parameters[supercell[idx_shape].embedded].B114 )*convolutions_char_strain_region_reduced_local.char_xx_yz[idx_basis_local]
								+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region_reduced_local.char_yy_yz[idx_basis_local]
								+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region_reduced_local.char_zz_yz[idx_basis_local]
								+ 4.0*( parameters[idx_shape].B156 - parameters[supercell[idx_shape].embedded].B156 )*convolutions_char_strain_region_reduced_local.char_xy_xz[idx_basis_local];

					P_G.y += 2.0*( parameters[idx_shape].B114 - parameters[supercell[idx_shape].embedded].B114 )*convolutions_char_strain_region_reduced_local.char_yy_xz[idx_basis_local]
								+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region_reduced_local.char_xx_xz[idx_basis_local]
								+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region_reduced_local.char_zz_xz[idx_basis_local]
								+ 4.0*( parameters[idx_shape].B156 - parameters[supercell[idx_shape].embedded].B156 )*convolutions_char_strain_region_reduced_local.char_xy_yz[idx_basis_local];

					P_G.z += 2.0*( parameters[idx_shape].B114 - parameters[supercell[idx_shape].embedded].B114 )*convolutions_char_strain_region_reduced_local.char_zz_xy[idx_basis_local]
								+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region_reduced_local.char_xx_xy[idx_basis_local]
								+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region_reduced_local.char_yy_xy[idx_basis_local]
								+ 4.0*( parameters[idx_shape].B156 - parameters[supercell[idx_shape].embedded].B156 )*convolutions_char_strain_region_reduced_local.char_xz_yz[idx_basis_local];

				}

				break;

		}


		// Compute resulting Fourier transform of piezoelectric potential (using refractive index of matrix material to compute relative permittivity)

		result += ( -100.0*I*( Gx*P_G.x + Gy*P_G.y + Gz*P_G.z ) )/( 8.854188*pow( parameters[0].nr, 2.0 )*pow( G_mag, 2.0 ) );

	}

	return result;

} // End of function
