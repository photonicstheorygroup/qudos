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

void real_space_grid_linescan( int direction, int N_grid, struct supercell_geometry *supercell, double *x, double *y, double *z )
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


	// Populate real space grid along selected linescan direction

	int idx_grid;

	switch( direction )
	{

		case 0: // Linescan along [100]

			for( idx_grid = 1; idx_grid < N_grid - 1; idx_grid++ )
			{

				x[idx_grid] = x[idx_grid - 1] + dx;

			}

			for( idx_grid = 0; idx_grid < N_grid; idx_grid++ )
			{

				y[idx_grid] = ( z[idx_grid] = 0.0 );

			}

			break;


		case 1: // Linescan along [010]

			for( idx_grid = 1; idx_grid < N_grid - 1; idx_grid++ )
			{

				y[idx_grid] = y[idx_grid - 1] + dy;

			}

			for( idx_grid = 0; idx_grid < N_grid; idx_grid++ )
			{

				x[idx_grid] = ( z[idx_grid] = 0.0 );

			}

			break;


		case 2: // Linescan along [001]

			for( idx_grid = 1; idx_grid < N_grid - 1; idx_grid++ )
			{

				z[idx_grid] = z[idx_grid - 1] + dz;

			}

			for( idx_grid = 0; idx_grid < N_grid; idx_grid++ )
			{

				x[idx_grid] = ( y[idx_grid] = 0.0 );

			}

			break;


		case 3: // Linescan along [110]

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

			for( idx_grid = 1; idx_grid < N_grid - 1; idx_grid++ )
			{

				x[idx_grid] = x[idx_grid - 1] + dx;
				y[idx_grid] = y[idx_grid - 1] + dy;

			}

			for( idx_grid = 0; idx_grid < N_grid; idx_grid++ )
			{

				z[idx_grid] = 0;

			}

			break;


		case 4: // Linescan along [101]

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

			for( idx_grid = 1; idx_grid < N_grid - 1; idx_grid++ )
			{

				x[idx_grid] = x[idx_grid - 1] + dx;
				z[idx_grid] = z[idx_grid - 1] + dz;

			}

			for( idx_grid = 0; idx_grid < N_grid; idx_grid++ )
			{

				y[idx_grid] = 0.0;

			}

			break;


		case 5: // Linescan along [011]

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

			for( idx_grid = 1; idx_grid < N_grid - 1; idx_grid++ )
			{

				y[idx_grid] = y[idx_grid - 1] + dy;
				z[idx_grid] = z[idx_grid - 1] + dz;

			}

			for( idx_grid = 0; idx_grid < N_grid; idx_grid++ )
			{

				x[idx_grid] = 0.0;

			}

			break;


		case 6: // Linescan along [111]

			if ( dx < dy && dx < dz) // Identifying limiting lengths in case of non-cuboidal supercell
			{

				dy = dx;
				dz = dx;

				y[0] = x[0];
				z[0] = x[0];

				y[N_grid - 1] = x[N_grid - 1];
				z[N_grid - 1] = x[N_grid - 1];

			}
			else if ( dy < dz )
			{

				dx = dy;
				dz = dy;

				x[0] = y[0];
				z[0] = y[0];

				x[N_grid - 1] = y[N_grid - 1];
				z[N_grid - 1] = y[N_grid - 1];

			}
			else
			{

				dy = dz;
				dx = dz;

				x[0] = z[0];
				y[0] = z[0];

				x[N_grid - 1] = z[N_grid - 1];
				y[N_grid - 1] = z[N_grid - 1];

			}

			for( idx_grid = 1; idx_grid < N_grid - 1; idx_grid++ )
			{

				x[idx_grid] = x[idx_grid - 1] + dx;
				y[idx_grid] = y[idx_grid - 1] + dy;
				z[idx_grid] = z[idx_grid - 1] + dz;

			}

			break;

	}

}
