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

void write_carrier_localisation( int N_shape, int *N_state, complex double *localisation_vb, complex double *localisation_cb )
{

	int idx_state, idx_shape;
	FILE *fp1;

	// Carrier localisation in each shape is written in terms of % and should sum to 100.0 for a given eigenstate. 
	// The rows represent the states, and the columns represent the shapes (first column gives the percentage inside the electronic supercell, but outside any other shape in the supercell), i.e. entries in a given row of the output file should sum to 100.0%, number of columns in a row = number of shapes in supercell (including supercell shape))

	// Hole states

	fp1 = fopen( "carrier_localisation_vb.dat", "w" );

	for( idx_state = 2*N_state[0] - 1; idx_state >= 0; idx_state-- )
	{

		for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
		{

			fprintf( fp1, "%13e\t", 100.0*creal( localisation_vb[idx_shape + N_shape*idx_state] ) );

		}

		fprintf( fp1, "\n" );

	}

	fclose( fp1 );


	// Electron states

	fp1 = fopen( "carrier_localisation_cb.dat", "w" );

	for( idx_state = 0; idx_state < 2*N_state[1]; idx_state++ )
	{

		for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
		{

			fprintf( fp1, "%13e\t", 100.0*creal( localisation_cb[idx_shape + N_shape*idx_state] ) );

		}

		fprintf( fp1, "\n" );

	}

	fclose( fp1 );


}
