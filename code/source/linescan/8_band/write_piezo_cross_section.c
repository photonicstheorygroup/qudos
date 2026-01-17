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

void write_piezo_cross_section( int direction, int N_grid, double *x, double *y, double *z, int piezo_order, complex double *piezo_r )
{

	// Write out cross-section of the piezoelectric potential for the specified normal direction

	int idx_grid_1, idx_grid_2, size = 80;
	char str1[size], str2[size];
	FILE *fp1;


	switch( piezo_order )
	{

		case 1: // Potential due to first order piezoelectric polarisation only

			strcpy( str1, "piezo_cross_section_1st_order_" );

			break;

		case 2: // Potential due to first and second order piezoelectric polarisation

			strcpy( str1, "piezo_cross_section_1st_and_2nd_order_" );

			break;

		case 3: // Potential due to second order piezoelectric polarisation only

			strcpy( str1, "piezo_cross_section_2nd_order_" );

			break;

	}


	switch( direction )
	{

		case 0: // Normal along [100]

			size = sprintf( str2, "100.dat" );

			break;

		case 1: // Normal along [010]

			size = sprintf( str2, "010.dat" );

			break;

		case 2: // Normal along [001]

			size = sprintf( str2, "001.dat" );

			break;

		case 3: // Normal along [-110]

			size = sprintf( str2, "-110.dat" );

			break;

		case 4: // Normal along [110]

			size = sprintf( str2, "110.dat" );

			break;

		case 5: // Normal along [101]

			size = sprintf( str2, "-101.dat" );

			break;

		case 6: // Normal along [011]

			size = sprintf( str2, "0-11.dat" );

			break;

		case 7: // Normal along [111]

			size = sprintf( str2, "111.dat" );

			break;

	}


	strcat( str1, str2 );
	fp1 = fopen( str1, "w" );

	for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
	{

		// Column 1 - x    component of position vector         [nm]
		// Column 2 - y    component of position vector         [nm]
		// Column 3 - z    component of position vector         [nm]
		// Column 4 - Real      part of piezoelectric potential [eV]
		// Column 5 - Imaginary part of piezoelectric potential [eV]
		
		for ( idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++ ){
		
			fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t%17.10e\t%17.10e\n", 0.1*x[idx_grid_1*N_grid+idx_grid_2], 0.1*y[idx_grid_1*N_grid+idx_grid_2], 0.1*z[idx_grid_1*N_grid+idx_grid_2], creal( piezo_r[idx_grid_1*N_grid+idx_grid_2] ), cimag( piezo_r[idx_grid_1*N_grid+idx_grid_2] ) );
		
		}
		
		fprintf (fp1, "\n");

	}
	

	fclose( fp1 );

}
