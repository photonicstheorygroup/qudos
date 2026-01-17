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

void write_energies( int *N_states, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb )
{

     int idx_state, size = 80;
	char *str1, *str2;
	str1 = malloc( sizeof(char)*size );
	str2 = malloc( sizeof(char)*size );
     FILE *fp1;


     // 1. Write valaence band energies and bloch characters

     strcpy( str1, "energies_vb" );
     size = sprintf( str2, ".dat" );
     strcat( str1, str2 );
     fp1 = fopen( str1, "w" );

     for( idx_state = 2*N_states[0] - 1; idx_state > -1; idx_state-- )
     {

          // Column 1  - Energy                                              [eV]
          // Column 2  - Bloch character - SO       			     [%]                         
          // Column 3  - Bloch character - LH                                [%]
          // Column 4  - Bloch character - HH                                [%]
          // Column 5  - Bloch character - CB                                [%]
     

          fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t%17.10e\t%17.10e\n",
                                                  eigenstate_vb[idx_state].E,
                                                  100.0*( eigenstate_vb[idx_state].bloch[3] + eigenstate_vb[idx_state].bloch[7] ),
                                                  100.0*( eigenstate_vb[idx_state].bloch[2] + eigenstate_vb[idx_state].bloch[6] ),
                                                  100.0*( eigenstate_vb[idx_state].bloch[1] + eigenstate_vb[idx_state].bloch[5] ),
                                                  100.0*( eigenstate_vb[idx_state].bloch[0] + eigenstate_vb[idx_state].bloch[4] ) );

     }

     fclose( fp1 );


     // 1. Write conduction band energies and bloch characters

     strcpy( str1, "energies_cb" );
     size = sprintf( str2, ".dat" );
     strcat( str1, str2 );
     fp1 = fopen( str1, "w" );

     for( idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
     {

          // Column 1  - Energy                                              [eV]
          // Column 2  - Bloch character - SO       			     [%]                         
          // Column 3  - Bloch character - LH                                [%]
          // Column 4  - Bloch character - HH                                [%]
          // Column 5  - Bloch character - CB                                [%]
     

          fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t%17.10e\t%17.10e\n",
                                                  eigenstate_cb[idx_state].E,
                                                  100.0*( eigenstate_cb[idx_state].bloch[3] + eigenstate_cb[idx_state].bloch[7] ),
                                                  100.0*( eigenstate_cb[idx_state].bloch[2] + eigenstate_cb[idx_state].bloch[6] ),
                                                  100.0*( eigenstate_cb[idx_state].bloch[1] + eigenstate_cb[idx_state].bloch[5] ),
                                                  100.0*( eigenstate_cb[idx_state].bloch[0] + eigenstate_cb[idx_state].bloch[4] ) );

     }

     fclose( fp1 );

     free( str1 );
     free( str2 );

}
