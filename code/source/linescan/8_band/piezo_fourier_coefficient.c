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

complex double piezo_fourier_coefficient( int piezo_order, int idx_shape, int idx_Gx, int idx_Gy, int idx_Gz, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct characteristic_functions *characteristic_function_G, struct strain_tensor strain_G, struct convolution_char_strain convolutions_char_strain_region )
{

	complex double P_Gx = 0.0 + I*0.0, P_Gy = 0.0 + I*0.0, P_Gz = 0.0 + I*0.0, result = 0.0 + I*0.0; // Polarisation vector components and resulting contribution to piezoelectric potential

	// Wave vector indices
	int idx_basis_Gx = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_Gx;
	int idx_basis_Gy = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_Gy;
	int idx_basis_Gz = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_Gz;

	// Wave vector components
	double Gx = basis.Gx[idx_basis_Gx];
	double Gy = basis.Gy[idx_basis_Gy];
	double Gz = basis.Gz[idx_basis_Gz];

	// Wave vector index for stored strain and convolution arrays
	int idx_basis = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_Gz + basis.N_Gz_T*(basis.N_Gy_T/2 - basis.N_Gy/2 + idx_Gy) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 - basis.N_Gx/2 + idx_Gx);


	// Compute Fourier transform of piezoelectric polarisation

	if( fabs( pow( Gx, 2.0 ) + pow( Gy, 2.0 ) + pow( Gz, 2.0 ) ) > EPS ) // |G| != 0
	{

		switch( piezo_order )
		{

			case 1: // First order

				if ( idx_shape == 0 )
				{

					P_Gx += 2.0*parameters[0].e14*strain_G.yz[idx_basis];
					P_Gy += 2.0*parameters[0].e14*strain_G.xz[idx_basis];
					P_Gz += 2.0*parameters[0].e14*strain_G.xy[idx_basis];

				}
				else
				{

					P_Gx += 2.0*( parameters[idx_shape].e14 - parameters[supercell[idx_shape].embedded].e14 )*convolutions_char_strain_region.char_yz[idx_basis];
					P_Gy += 2.0*( parameters[idx_shape].e14 - parameters[supercell[idx_shape].embedded].e14 )*convolutions_char_strain_region.char_xz[idx_basis];
					P_Gz += 2.0*( parameters[idx_shape].e14 - parameters[supercell[idx_shape].embedded].e14 )*convolutions_char_strain_region.char_xy[idx_basis];

				}

				break;

			case 2: // Second order

				if ( idx_shape == 0 )
				{

					P_Gx += 2.0*parameters[0].B114*convolutions_char_strain_region.xx_yz[idx_basis]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region.yy_yz[idx_basis]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region.zz_yz[idx_basis]
									+ 4.0*parameters[0].B156*convolutions_char_strain_region.xy_xz[idx_basis];

					P_Gy += 2.0*parameters[0].B114*convolutions_char_strain_region.yy_xz[idx_basis]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region.xx_xz[idx_basis]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region.zz_xz[idx_basis]
									+ 4.0*parameters[0].B156*convolutions_char_strain_region.xy_yz[idx_basis];

					P_Gz += 2.0*parameters[0].B114*convolutions_char_strain_region.zz_xy[idx_basis]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region.xx_xy[idx_basis]
									+ 2.0*parameters[0].B124*convolutions_char_strain_region.yy_xy[idx_basis]
									+ 4.0*parameters[0].B156*convolutions_char_strain_region.xz_yz[idx_basis];

				}
				else
				{

					P_Gx += 2.0*( parameters[idx_shape].B114 - parameters[supercell[idx_shape].embedded].B114 )*convolutions_char_strain_region.char_xx_yz[idx_basis]
						+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region.char_yy_yz[idx_basis]
						+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region.char_zz_yz[idx_basis]
						+ 4.0*( parameters[idx_shape].B156 - parameters[supercell[idx_shape].embedded].B156 )*convolutions_char_strain_region.char_xy_xz[idx_basis];

					P_Gy += 2.0*( parameters[idx_shape].B114 - parameters[supercell[idx_shape].embedded].B114 )*convolutions_char_strain_region.char_yy_xz[idx_basis]
						+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region.char_xx_xz[idx_basis]
						+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region.char_zz_xz[idx_basis]
						+ 4.0*( parameters[idx_shape].B156 - parameters[supercell[idx_shape].embedded].B156 )*convolutions_char_strain_region.char_xy_yz[idx_basis];

					P_Gz += 2.0*( parameters[idx_shape].B114 - parameters[supercell[idx_shape].embedded].B114 )*convolutions_char_strain_region.char_zz_xy[idx_basis]
						+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region.char_xx_xy[idx_basis]
						+ 2.0*( parameters[idx_shape].B124 - parameters[supercell[idx_shape].embedded].B124 )*convolutions_char_strain_region.char_yy_xy[idx_basis]
						+ 4.0*( parameters[idx_shape].B156 - parameters[supercell[idx_shape].embedded].B156 )*convolutions_char_strain_region.char_xz_yz[idx_basis];

				}

				break;

		}


		// Compute Fourier transform of piezoelectric potential (using refractive index of matrix material to compute relative permittivity)

		result += ( -100.0*I*( Gx*P_Gx + Gy*P_Gy + Gz*P_Gz ) )/( 8.854188*pow( parameters[0].nr, 2.0 )*( pow( Gx, 2.0 ) + pow( Gy, 2.0 ) + pow( Gz, 2.0 ) ) );

	}

	return result;

} // End of function
