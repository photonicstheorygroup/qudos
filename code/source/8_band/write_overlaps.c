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

void write_overlaps( int *N_states, complex double *overlaps )
{

    int size = 80;
	char *str1;
	str1 = malloc( sizeof(char)*size );
    FILE *fp1;


    // 1. Write overlaps


    strcpy( str1, "overlaps.dat" );
          
	fp1 = fopen( str1, "w" );

    for( int idx_state_cb = 0; idx_state_cb < 2*N_states[1]; idx_state_cb++ )
    {

		for ( int idx_state_vb = 0; idx_state_vb < 2*N_states[0]; idx_state_vb++ )
		{
		  
		       // Column 1  - CB state index 
		       // Column 2  - VB state index      			     	                       
		       // Column 3  - Overlap real
		       // Column 4  - Overlap complex
	      

		    fprintf( fp1, "%d\t%d\t%13e\t%13e\n",
		                                            idx_state_cb,
		                                            idx_state_vb,
		                                            creal( overlaps[idx_state_cb*2*N_states[0] + idx_state_vb] ), cimag( overlaps[idx_state_cb*2*N_states[0] + idx_state_vb] ));
                
        }
                 
        fprintf( fp1, "\n" );

    }

    fclose( fp1 );
          
    free( str1 );

}
