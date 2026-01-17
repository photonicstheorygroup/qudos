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

void real_space_grid_cross_section( int direction, int N_grid, struct supercell_geometry *supercell, double *x, double *y, double *z )
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


	// Populate real space grid for cross-sections with selected normal direction

	int idx_grid_1, idx_grid_2;

	switch( direction )
	{

		case 0: // Normal along [100]

			for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
			{
				
				for( idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++ )
				{
					y[(idx_grid_1*N_grid) + idx_grid_2] = y[0] + idx_grid_1 * dy;	
					z[(idx_grid_1*N_grid) + idx_grid_2] = z[0] + idx_grid_2 * dz;
					x[(idx_grid_1*N_grid) + idx_grid_2] = 0;

				}

			}

			break;


		case 1: // Normal along [010]

			for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
			{
				
				for( idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++ )
				{
					x[(idx_grid_1*N_grid) + idx_grid_2] = x[0] + idx_grid_1 * dx;	
					z[(idx_grid_1*N_grid) + idx_grid_2] = z[0] + idx_grid_2 * dz;
					y[(idx_grid_1*N_grid) + idx_grid_2] = 0;

				}

			}

			break;


		case 2: // Normal along [001]

			for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
			{
				
				for( idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++ )
				{
					x[(idx_grid_1*N_grid) + idx_grid_2] = x[0] + idx_grid_1 * dx;	
					y[(idx_grid_1*N_grid) + idx_grid_2] = y[0] + idx_grid_2 * dy;
					z[(idx_grid_1*N_grid) + idx_grid_2] = 0;

				}

			}

			break;


		case 3: // Normal along [-110]

			if ( dx < dy ) // Identifying limiting length in case of non-square supercell in (001) plane
			{

				dy = dx;

				y[0] = x[0];

				y[N_grid - 1] = x[N_grid - 1];

			}
			else
			{

				dx = dy;

				x[0] = y[0];

				x[N_grid - 1] = y[N_grid - 1];

			}

			for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
			{
				
				for( idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++ )
				{
					y[(idx_grid_1*N_grid) + idx_grid_2] = y[0] + idx_grid_1 * dy;	
					z[(idx_grid_1*N_grid) + idx_grid_2] = z[0] + idx_grid_2 * dz;
					x[(idx_grid_1*N_grid) + idx_grid_2] = x[0] + idx_grid_1 * dx;

				}

			}

			break;


		case 4: // Normal along [110]

			if ( dx < dy ) // Identifying limiting length in case of non-square supercell in (001) plane
			{

				dy = dx;

				y[0] = x[0];

				y[N_grid - 1] = x[N_grid - 1];

			}
			else
			{

				dx = dy;

				x[0] = y[0];

				x[N_grid - 1] = y[N_grid - 1];

			}

			for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
			{
				
				for( idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++ )
				{
					y[(idx_grid_1*N_grid) + idx_grid_2] = y[0] + idx_grid_1 * dy;	
					z[(idx_grid_1*N_grid) + idx_grid_2] = z[0] + idx_grid_2 * dz;
					x[(idx_grid_1*N_grid) + idx_grid_2] = x[N_grid - 1] - idx_grid_1 * dx;

				}

			}

			break;


		case 5: // Normal along [-101]

			if ( dx < dz ) // Identifying limiting length in case of non-square supercell in (010) plane
			{

				dz = dx;

				z[0] = x[0];

				z[N_grid - 1] = x[N_grid - 1];

			}
			else
			{

				dx = dz;

				x[0] = z[0];

				x[N_grid - 1] = z[N_grid - 1];

			}

			for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
			{
				
				for( idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++ )
				{
					y[(idx_grid_1*N_grid) + idx_grid_2] = y[0] + idx_grid_2 * dy;	
					z[(idx_grid_1*N_grid) + idx_grid_2] = z[0] + idx_grid_1 * dz;
					x[(idx_grid_1*N_grid) + idx_grid_2] = x[0] + idx_grid_1 * dx;

				}

			}

			break;


		case 6: // Normal along [0-11]

			if ( dy < dz ) // Identifying limiting length in case of non-square supercell in (100) plane
			{

				dz = dy;

				z[0] = y[0];

				z[N_grid - 1] = y[N_grid - 1];

			}
			else
			{

				dy = dz;

				y[0] = z[0];

				y[N_grid - 1] = z[N_grid - 1];

			}

			for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
			{
				
				for( idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++ )
				{
					y[(idx_grid_1*N_grid) + idx_grid_2] = y[0] + idx_grid_1 * dy;	
					z[(idx_grid_1*N_grid) + idx_grid_2] = z[0] + idx_grid_1 * dz;
					x[(idx_grid_1*N_grid) + idx_grid_2] = x[0] + idx_grid_2 * dx;

				}

			}

			break;


		case 7: // Normal along [111]  (Special case that assumes [100], [010], [001] dimensions are the same) An expanded cross-section is output to visualise the effect of the periodic boundary conditions on the (111) plane. 

			if ( dx < dy )
			{

				dy = dx;

				y[0] = x[0];

				y[N_grid - 1] = x[N_grid - 1];

			}
			else
			{

				dx = dy;

				x[0] = y[0];

				x[N_grid - 1] = y[N_grid - 1];

			}

			z[0] = 2*z[0];

			z[N_grid - 1] = 2*z[N_grid - 1];
			

			for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
			{
				
				for( idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++ )
				{
		
					z[(idx_grid_1*N_grid) + idx_grid_2] = z[0] + 2*idx_grid_2 * dz;
					y[(idx_grid_1*N_grid) + idx_grid_2] = -2*sqrt(0.5*supercell[0].L1*0.5*supercell[0].L1+0.5*supercell[0].L2*0.5*supercell[0].L2)/sqrt(2.0) - z[(idx_grid_1*N_grid) + idx_grid_2]/2 + 2*idx_grid_1 * dy;	
					x[(idx_grid_1*N_grid) + idx_grid_2] = -y[(idx_grid_1*N_grid) + idx_grid_2] - z[(idx_grid_1*N_grid) + idx_grid_2];

				}

			}

			break;

	}

}
