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

void real_space_grid_3D( int N_grid, struct supercell_geometry *supercell, double *x, double *y, double *z )
{

	// Set real space grid endpoints and step sizes

	x[0]          = -0.5*supercell[0].L1;
	x[N_grid - 1] =  0.5*supercell[0].L1;

	y[0]          = -0.5*supercell[0].L2;
	y[N_grid - 1] =  0.5*supercell[0].L2;

	z[0]          = -0.5*supercell[0].L3;
	z[N_grid - 1] =  0.5*supercell[0].L3;

	double dx = ( x[N_grid - 1] - x[0] )/( ( double ) N_grid - 1 );
	double dy = ( y[N_grid - 1] - y[0] )/( ( double ) N_grid - 1 );
	double dz = ( z[N_grid - 1] - z[0] )/( ( double ) N_grid - 1 );


	// Populate 3D real space grid

	int idx_grid_1, idx_grid_2, idx_grid_3;

	for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
	{
				
		for( idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++ )
		{

			for( idx_grid_3 = 0; idx_grid_3 < N_grid; idx_grid_3++ )
			{
				
				x[(idx_grid_1*N_grid*N_grid) + idx_grid_2*N_grid + idx_grid_3] = x[0] + idx_grid_1 * dx;
				y[(idx_grid_1*N_grid*N_grid) + idx_grid_2*N_grid + idx_grid_3] = y[0] + idx_grid_2 * dy;	
				z[(idx_grid_1*N_grid*N_grid) + idx_grid_2*N_grid + idx_grid_3] = z[0] + idx_grid_3 * dz;

			}

		}

	}

}
