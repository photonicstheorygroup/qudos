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

void collect_eigenstates( int *N_states, struct wave_vectors basis, double *eigval, complex double *eigvec, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb )
{

     // 1. Extract valence band energies and local components of envelope function Fourier coefficients from "eigval" and "eigvec"
     #pragma omp parallel for
     for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
     {

          // Energy
          eigenstate_vb[idx_state].E = eigval[idx_state];


          // Envelope function Fourier coefficients
          if ( ( (idx_state) / block_size ) % np_cols == my_pcol ) // Identifying if part of envelope function is on the process
          {

		     for( int idx_coeff = 0; idx_coeff < na_rows; idx_coeff++ )
		     {

		  	     eigenstate_vb[idx_state].fourier[idx_coeff] = eigvec[idx_coeff + na_rows*( ( block_size * ( ( idx_state / block_size ) / np_cols ) ) + ( idx_state % block_size ) )]; // Extracting local part of state from block-cyclic distribution

		     }
          
          }

     }


     // 2. Extract conduction band energies and local components of envelope function Fourier coefficients from "eigval" and "eigvec"
     #pragma omp parallel for
     for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
     {

          // Energy
          eigenstate_cb[idx_state].E = eigval[idx_state + 2*N_states[0]];


          // Envelope function Fourier coefficients
          if ( ( (idx_state + 2*N_states[0]) / block_size ) % np_cols == my_pcol ) // Identifying if part of envelope function is on the process
          {
          
		     for( int idx_coeff = 0; idx_coeff < na_rows; idx_coeff++ )
		     {

		  	     eigenstate_cb[idx_state].fourier[idx_coeff] = eigvec[idx_coeff + na_rows*( ( block_size * ( ( (idx_state + 2*N_states[0]) / block_size ) / np_cols ) ) + ( (idx_state + 2*N_states[0]) % block_size ) )]; // Extracting local part of state from block-cyclic distribution

		     }
          
          }

     }

}
