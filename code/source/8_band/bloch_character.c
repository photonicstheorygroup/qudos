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

void bloch_character( int band, int *N_states, struct wave_vectors basis, struct eigenstate *eigenstate )
{

     int N_G = basis.N_Gx * basis.N_Gy * basis.N_Gz;

     switch( band )
     {
	   
          //////////////////////
          //// Valence band ////
          //////////////////////

          case 0:

               for( int idx_band = 0; idx_band < N_BAND; idx_band++ )
               {

                    for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
                    {
                         
                         double result = 0;
			  
			          #pragma omp parallel for reduction(+:result) schedule(guided)
                         for( int idx_G = 0; idx_G < N_G; idx_G++ )
                         {
				
                              if ( ( (N_G * idx_band + idx_G) / block_size ) % np_rows == my_prow && ( (idx_state) / block_size ) % np_cols == my_pcol ) // Identifying if part of envelope function is on the process
                              {	
                              
                              result += cabs( eigenstate[idx_state].fourier[indxg2l( idx_G + N_G*idx_band, block_size, 0, 0, np_rows )] ) * cabs( eigenstate[idx_state].fourier[indxg2l( idx_G + N_G*idx_band, block_size, 0, 0, np_rows )] ); // Calculating contribution from process
                                   

                              }

                         }

                         MPI_Reduce( &result, &eigenstate[idx_state].bloch[idx_band], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ); // Summing contributions from all processes

                    }

               }

               break;


          /////////////////////////
          //// Conduction band ////
          /////////////////////////

          case 1:

               for( int idx_band = 0; idx_band < N_BAND; idx_band++ )
               {

                    for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
                    {

                         double result = 0;
                         
                         
                         #pragma omp parallel for reduction(+:result) schedule(guided)
                         for( int idx_G = 0; idx_G < N_G; idx_G++ )
                         {

                              if ( ( (N_G * idx_band + idx_G) / block_size ) % np_rows == my_prow && ( (idx_state + 2*N_states[0]) / block_size ) % np_cols == my_pcol ) // Identifying if part of envelope function is on the process
                              {	
                              result += cabs( eigenstate[idx_state].fourier[indxg2l( idx_G + N_G*idx_band, block_size, 0, 0, np_rows )] ) * cabs( eigenstate[idx_state].fourier[indxg2l( idx_G + N_G*idx_band, block_size, 0, 0, np_rows )] ); // Calculating contribution from process

                              }

                         }

                         MPI_Reduce( &result, &eigenstate[idx_state].bloch[idx_band], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ); // Summing contributions from all processes
                         
                    }

               }

               break;

     }

}
