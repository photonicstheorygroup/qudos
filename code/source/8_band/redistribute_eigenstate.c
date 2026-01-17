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

void redistribute_eigenstate( struct wave_vectors basis, int *N_states, int eigenstate_type, int idx_state, struct eigenstate *eigenstate_stored, struct eigenstate eigenstate_interim )
{

	int count = 0;

	if ( ( ( idx_state + ( eigenstate_type == 1 ? 2*N_states[0] : 0 ) ) / block_size ) % np_cols == my_pcol ) // If eigenstate is on process column
	{

		for( int idx_G = 0; idx_G < N_H; idx_G++ ) // Iterating over wave vector indices
		{

			if ( ( (idx_G) / block_size ) % np_rows == my_prow ) // If Fourier coefficient is on process row
			{	
			
				eigenstate_interim.fourier[idx_G] = eigenstate_stored[idx_state].fourier[count];

				count++;
					
			}

		}
		
	}

	MPI_Allreduce( MPI_IN_PLACE, eigenstate_interim.fourier, N_H, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD ); // Combine interim eigenstate contributions from all processes
	

}
