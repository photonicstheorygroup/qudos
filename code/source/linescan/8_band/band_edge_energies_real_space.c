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

void band_edge_energies_real_space( int N_grid, int N_shape, double *x, double *y, double *z, double Fx, double Fy, double Fz, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct characteristic_functions *characteristic_function_G, struct strain_tensor strain_r, complex double *piezo_r, struct band_edge_energies band_edges )
{

	int idx_grid, idx_shape = 0;

	mkl_set_num_threads(1);

	#pragma omp parallel for private( idx_shape, idx_grid )
	for( idx_grid = 0; idx_grid < N_grid; idx_grid++ )
	{

		double shape_frac[N_shape]; // A measure of the influence of each shape on the grid point in the case of Gaussian interface softening

		// Determine in which shape the current real space grid point resides

		if( softening == 0 ) // Hard shape interfaces
		{

			idx_shape = determine_idx_shape( N_shape, x[idx_grid], y[idx_grid], z[idx_grid], supercell );

		}
		else // Gaussian interface softening
		{

			for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
			{

				shape_frac[idx_shape] = characteristic_function_real_space( x[idx_grid], y[idx_grid], z[idx_grid], basis, supercell, characteristic_function_G[idx_shape] );

			}

		}


		// Set up and diagonalise the zone-centre bulk k.p Hamiltonian at the current real space grid point

		if( softening == 0 ) // Hard shape interfaces
		{

			eigenstates_real_space( idx_grid, idx_shape, x[idx_grid], y[idx_grid], z[idx_grid], Fx, Fy, Fz, parameters[idx_shape], strain_r, piezo_r, band_edges );

		}
		else // Gaussian interface softening
		{

			struct material_params interface_params;

			// Simple linear interpolation of material parameters based on "shape_frac"
			interpolate_interface( N_shape, shape_frac, supercell, parameters, &interface_params );

			eigenstates_real_space( idx_grid, idx_shape, x[idx_grid], y[idx_grid], z[idx_grid], Fx, Fy, Fz, interface_params, strain_r, piezo_r, band_edges );

		}

	}

}
