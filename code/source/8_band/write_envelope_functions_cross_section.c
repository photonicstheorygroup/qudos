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

void write_envelope_functions_cross_section( int *N_states, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb, int direction )
{

	int idx_state, idx_grid_1, idx_grid_2, size = 80;
	char *str1, *str2;
	str1 = malloc( sizeof(char)*size );
	str2 = malloc( sizeof(char)*size );
	FILE *fp1;


    // 1. Write valence band envelope functions
	for( idx_state = 2*N_states[0] - 1; idx_state > -1; idx_state-- ) // Write out valence bands in order of decreasing energy
	{

        strcpy( str1, "envelope_function_cross_section_vb_" );
        switch( direction )
		{

			case 0: // Normal along [100]

				size = sprintf( str2, "%d_100.dat", 2*N_states[0] - idx_state );

				break;

			case 1: // Normal along [010]

				size = sprintf( str2, "%d_010.dat", 2*N_states[0] - idx_state );

				break;

			case 2: // Normal along [001]

				size = sprintf( str2, "%d_001.dat", 2*N_states[0] - idx_state );

				break;

			case 3: // Normal along [-110]

				size = sprintf( str2, "%d_-110.dat", 2*N_states[0] - idx_state );

				break;

			case 4: // Normal along [110]

				size = sprintf( str2, "%d_110.dat", 2*N_states[0] - idx_state );

				break;

			case 5: // Normal along [-101]

				size = sprintf( str2, "%d_-101.dat", 2*N_states[0] - idx_state );

				break;

			case 6: // Normal along [0-11]

				size = sprintf( str2, "%d_0-11.dat", 2*N_states[0] - idx_state );

				break;

			case 7: // Normal along [111]

				size = sprintf( str2, "%d_111.dat", 2*N_states[0] - idx_state );

				break;

		}
    	
		strcat( str1, str2 );
		fp1 = fopen( str1, "w" );

		for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
		{

			// Column 1  - x position                                          [nm]
            // Column 2  - y position                                          [nm]
            // Column 3  - z position                                          [nm]
            // Column 4  - Energy                                              [eV]
            // Column 5  - Real      part of SO S1 component of envelope function [nm^{-3/2}]
			// Column 6  - Imaginary part of SO S1 component of envelope function [nm^{-3/2}]
			// Column 7  - Real      part of SO S2 component of envelope function [nm^{-3/2}]
			// Column 8  - Imaginary part of SO S2 component of envelope function [nm^{-3/2}]
			// Column 9  - Real      part of LH S1 component of envelope function [nm^{-3/2}]
			// Column 10  - Imaginary part of LH S1 component of envelope function [nm^{-3/2}]
			// Column 11  - Real      part of LH S2 component of envelope function [nm^{-3/2}]
			// Column 12  - Imaginary part of LH S2 component of envelope function [nm^{-3/2}]
			// Column 13  - Real      part of HH S1 component of envelope function [nm^{-3/2}]
			// Column 14  - Imaginary part of HH S1 component of envelope function [nm^{-3/2}]
			// Column 15  - Real      part of HH S2 component of envelope function [nm^{-3/2}]
			// Column 16  - Imaginary part of HH S2 component of envelope function [nm^{-3/2}]
			// Column 17  - Real      part of CB S1 component of envelope function [nm^{-3/2}]
			// Column 18 - Imaginary part of CB S1 component of envelope function [nm^{-3/2}]
			// Column 19  - Real      part of CB S2 component of envelope function [nm^{-3/2}]
			// Column 20 - Imaginary part of CB S2 component of envelope function [nm^{-3/2}]
			
			for (idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++)
			{
		
		       	fprintf( fp1, "%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\n",
		                                                 0.1*x[idx_grid_1*N_grid+idx_grid_2], 0.1*y[idx_grid_1*N_grid+idx_grid_2], 0.1*z[idx_grid_1*N_grid+idx_grid_2], eigenstate_vb[idx_state].E,
		                                            sqrt( 1000.0 )*creal( eigenstate_vb[idx_state].F_4[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_vb[idx_state].F_4[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_vb[idx_state].F_8[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_vb[idx_state].F_8[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_vb[idx_state].F_3[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_vb[idx_state].F_3[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_vb[idx_state].F_7[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_vb[idx_state].F_7[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_vb[idx_state].F_2[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_vb[idx_state].F_2[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_vb[idx_state].F_6[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_vb[idx_state].F_6[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_vb[idx_state].F_1[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_vb[idx_state].F_1[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_vb[idx_state].F_5[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_vb[idx_state].F_5[idx_grid_1*N_grid+idx_grid_2] ) );
		                                                         
			}
			
			fprintf (fp1, "\n");

		}

		fclose( fp1 );

	}


    // 2. Write conduction band envelope functions
	for( idx_state = 0; idx_state < 2*N_states[1]; idx_state++ ) // Write out conduction bands in order of increasing energy
	{

		strcpy( str1, "envelope_function_cross_section_cb_" );
		switch( direction )
		{

			case 0: // Normal along [100]

				size = sprintf( str2, "%d_100.dat", idx_state + 1 );

				break;

			case 1: // Normal along [010]

				size = sprintf( str2, "%d_010.dat", idx_state + 1 );

				break;

			case 2: // Normal along [001]

				size = sprintf( str2, "%d_001.dat", idx_state + 1 );

				break;

			case 3: // Normal along [-110]

				size = sprintf( str2, "%d_-110.dat", idx_state + 1 );

				break;

			case 4: // Normal along [110]

				size = sprintf( str2, "%d_110.dat", idx_state + 1 );

				break;

			case 5: // Normal along [-101]

				size = sprintf( str2, "%d_-101.dat", idx_state + 1 );

				break;

			case 6: // Normal along [0-11]

				size = sprintf( str2, "%d_0-11.dat", idx_state + 1 );

				break;

			case 7: // Normal along [111]

				size = sprintf( str2, "%d_111.dat", idx_state + 1 );

				break;

		}

		strcat( str1, str2 );
		fp1 = fopen( str1, "w" );

		for( idx_grid_1 = 0; idx_grid_1 < N_grid; idx_grid_1++ )
		{

			// Column 1  - x position                                          [nm]
            // Column 2  - y position                                          [nm]
            // Column 3  - z position                                          [nm]
            // Column 4  - Energy                                              [eV]
            // Column 5  - Real      part of SO S1 component of envelope function [nm^{-3/2}]
			// Column 6  - Imaginary part of SO S1 component of envelope function [nm^{-3/2}]
			// Column 7  - Real      part of SO S2 component of envelope function [nm^{-3/2}]
			// Column 8  - Imaginary part of SO S2 component of envelope function [nm^{-3/2}]
			// Column 9  - Real      part of LH S1 component of envelope function [nm^{-3/2}]
			// Column 10  - Imaginary part of LH S1 component of envelope function [nm^{-3/2}]
			// Column 11  - Real      part of LH S2 component of envelope function [nm^{-3/2}]
			// Column 12  - Imaginary part of LH S2 component of envelope function [nm^{-3/2}]
			// Column 13  - Real      part of HH S1 component of envelope function [nm^{-3/2}]
			// Column 14  - Imaginary part of HH S1 component of envelope function [nm^{-3/2}]
			// Column 15  - Real      part of HH S2 component of envelope function [nm^{-3/2}]
			// Column 16  - Imaginary part of HH S2 component of envelope function [nm^{-3/2}]
			// Column 17  - Real      part of CB S1 component of envelope function [nm^{-3/2}]
			// Column 18 - Imaginary part of CB S1 component of envelope function [nm^{-3/2}]
			// Column 19  - Real      part of CB S2 component of envelope function [nm^{-3/2}]
			// Column 20 - Imaginary part of CB S2 component of envelope function [nm^{-3/2}]
               
			for (idx_grid_2 = 0; idx_grid_2 < N_grid; idx_grid_2++)
			{
			
				fprintf( fp1, "%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t%13e\n",
		                                                 0.1*x[idx_grid_1*N_grid+idx_grid_2], 0.1*y[idx_grid_1*N_grid+idx_grid_2], 0.1*z[idx_grid_1*N_grid+idx_grid_2], eigenstate_cb[idx_state].E,
		                                            sqrt( 1000.0 )*creal( eigenstate_cb[idx_state].F_4[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_cb[idx_state].F_4[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_cb[idx_state].F_8[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_cb[idx_state].F_8[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_cb[idx_state].F_3[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_cb[idx_state].F_3[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_cb[idx_state].F_7[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_cb[idx_state].F_7[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_cb[idx_state].F_2[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_cb[idx_state].F_2[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_cb[idx_state].F_6[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_cb[idx_state].F_6[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_cb[idx_state].F_1[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_cb[idx_state].F_1[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*creal( eigenstate_cb[idx_state].F_5[idx_grid_1*N_grid+idx_grid_2] ),
		                                            sqrt( 1000.0 )*cimag( eigenstate_cb[idx_state].F_5[idx_grid_1*N_grid+idx_grid_2] ) );
																	
			}
		
			fprintf (fp1, "\n");

		}

		fclose( fp1 );

	}

	free( str1 );
    free( str2 );

}
