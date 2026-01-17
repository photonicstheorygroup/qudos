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

void envelope_functions( struct supercell_geometry *supercell, int *N_states, int type, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct eigenstate *eigenstate, int direction )
{
	
	struct eigenstate *eigenstate_temp = malloc( 2*N_states[type]*sizeof( _eigenstate ) );
     	
	for( int idx_state = 0; idx_state < 2*N_states[type]; idx_state++ )
	{

		eigenstate_temp[idx_state].probability_1 = calloc( N_grid, sizeof(         double ) );
		eigenstate_temp[idx_state].probability_2 = calloc( N_grid, sizeof(         double ) );
		eigenstate_temp[idx_state].probability_3 = calloc( N_grid, sizeof(         double ) );
		eigenstate_temp[idx_state].probability_4 = calloc( N_grid, sizeof(         double ) );
		eigenstate_temp[idx_state].probability_5 = calloc( N_grid, sizeof(         double ) );
		eigenstate_temp[idx_state].probability_6 = calloc( N_grid, sizeof(         double ) );
		eigenstate_temp[idx_state].probability_7 = calloc( N_grid, sizeof(         double ) );
		eigenstate_temp[idx_state].probability_8 = calloc( N_grid, sizeof(         double ) );
		eigenstate_temp[idx_state].F_1           = calloc( N_grid, sizeof( complex double ) );
		eigenstate_temp[idx_state].F_2           = calloc( N_grid, sizeof( complex double ) );
		eigenstate_temp[idx_state].F_3           = calloc( N_grid, sizeof( complex double ) );
		eigenstate_temp[idx_state].F_4           = calloc( N_grid, sizeof( complex double ) );
		eigenstate_temp[idx_state].F_5           = calloc( N_grid, sizeof( complex double ) );
		eigenstate_temp[idx_state].F_6           = calloc( N_grid, sizeof( complex double ) );
		eigenstate_temp[idx_state].F_7           = calloc( N_grid, sizeof( complex double ) );
		eigenstate_temp[idx_state].F_8           = calloc( N_grid, sizeof( complex double ) );

	}


	// Compute real space envelope function Bloch components at current real space grid point
	#pragma omp parallel for collapse(2)	
	for( int idx_grid = 0; idx_grid < N_grid; idx_grid++ ) // Looping over real space grid points
	{
		 
		for( int local_col_index = 0; local_col_index < na_colsz; local_col_index++ ) // Local eigenvector columns (states)
		{
		
			int global_col_index = indxl2g(local_col_index, block_size, my_pcol, 0, np_cols); // Identifying state
			
			if ( ( type == 0 && global_col_index < 2*N_states[0] ) || ( type == 1 && global_col_index >= 2*N_states[0] )) // Determining if state is CB or VB
			{
			
				for( int local_row_index = 0; local_row_index < na_rowsz; local_row_index++ ) // Local wave vector index
				{
						
					int global_row_index = indxl2g(local_row_index, block_size, my_prow, 0, np_rows); // Global wave vector index
						
					int idx_state;
					
					if ( type == 1 )
					{
					
						idx_state = global_col_index - 2*N_states[0];
					
					}
					
					else
					{
					
						idx_state = global_col_index;
					
					}
					
					int band = global_row_index / (N_H / 8); // Determining band
					
					// Determinig wave vector component indices
					int idx_Gx = (global_row_index % (N_H / 8)) / (basis.N_Gy*basis.N_Gz); 
					int idx_Gy = ((global_row_index % (N_H / 8)) % (basis.N_Gy*basis.N_Gz)) / basis.N_Gz;
					int idx_Gz = ((global_row_index % (N_H / 8)) % (basis.N_Gy*basis.N_Gz)) % basis.N_Gz;
					
					// Determinig Fourier component index corresponding to wave vector
					int idx_basis_Gx = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_Gx;
					int idx_basis_Gy = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_Gy;
					int idx_basis_Gz = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_Gz;

					complex double exp_part = cexp( I*( ( basis.Gx[idx_basis_Gx] + basis.Kx )*x[idx_grid] + ( basis.Gy[idx_basis_Gy] + basis.Ky )*y[idx_grid] + ( basis.Gz[idx_basis_Gz] + basis.Kz )*z[idx_grid] ) ) / sqrt( supercell[0].L1*supercell[0].L2*supercell[0].L3 ); // Exponential part
					
					
					switch( band )
					{
					
						case 0:
						
							eigenstate_temp[idx_state].F_1[idx_grid] += eigenstate[idx_state].fourier[local_row_index]*exp_part;
							
							break;
							
						case 1:
						
							eigenstate_temp[idx_state].F_2[idx_grid] += eigenstate[idx_state].fourier[local_row_index]*exp_part;
							
							break;
							
						case 2:
						
							eigenstate_temp[idx_state].F_3[idx_grid] += eigenstate[idx_state].fourier[local_row_index]*exp_part;
							
							break;
							
						case 3:
						
							eigenstate_temp[idx_state].F_4[idx_grid] += eigenstate[idx_state].fourier[local_row_index]*exp_part;
							
							break;
							
						case 4:
						
							eigenstate_temp[idx_state].F_5[idx_grid] += eigenstate[idx_state].fourier[local_row_index]*exp_part;
							
							break;
							
						case 5:
						
							eigenstate_temp[idx_state].F_6[idx_grid] += eigenstate[idx_state].fourier[local_row_index]*exp_part;
							
							break;
							
						case 6:
						
							eigenstate_temp[idx_state].F_7[idx_grid] += eigenstate[idx_state].fourier[local_row_index]*exp_part;
							
							break;
							
						case 7:
						
							eigenstate_temp[idx_state].F_8[idx_grid] += eigenstate[idx_state].fourier[local_row_index]*exp_part;
							
							break;
							
					
					}
				
				}
				
			}
		
		}

	}
	
	for( int idx_state = 0; idx_state < 2*N_states[type]; idx_state++ ) // Adding contributions from each process
	{
		
		MPI_Reduce( eigenstate_temp[idx_state].F_1, eigenstate[idx_state].F_1, N_grid, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( eigenstate_temp[idx_state].F_2, eigenstate[idx_state].F_2, N_grid, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( eigenstate_temp[idx_state].F_3, eigenstate[idx_state].F_3, N_grid, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( eigenstate_temp[idx_state].F_4, eigenstate[idx_state].F_4, N_grid, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( eigenstate_temp[idx_state].F_5, eigenstate[idx_state].F_5, N_grid, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( eigenstate_temp[idx_state].F_6, eigenstate[idx_state].F_6, N_grid, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( eigenstate_temp[idx_state].F_7, eigenstate[idx_state].F_7, N_grid, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( eigenstate_temp[idx_state].F_8, eigenstate[idx_state].F_8, N_grid, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );
		
		if ( myid == 0)
		{

			// Determining the probability density for each Bloch band
			for( int idx_grid = 0; idx_grid < N_grid; idx_grid++ )
			{
			
				eigenstate[idx_state].probability_1[idx_grid] = pow( cabs( eigenstate[idx_state].F_1[idx_grid] ), 2.0 );
				eigenstate[idx_state].probability_2[idx_grid] = pow( cabs( eigenstate[idx_state].F_2[idx_grid] ), 2.0 );
				eigenstate[idx_state].probability_3[idx_grid] = pow( cabs( eigenstate[idx_state].F_3[idx_grid] ), 2.0 );
				eigenstate[idx_state].probability_4[idx_grid] = pow( cabs( eigenstate[idx_state].F_4[idx_grid] ), 2.0 );
				eigenstate[idx_state].probability_5[idx_grid] = pow( cabs( eigenstate[idx_state].F_5[idx_grid] ), 2.0 );
				eigenstate[idx_state].probability_6[idx_grid] = pow( cabs( eigenstate[idx_state].F_6[idx_grid] ), 2.0 );
				eigenstate[idx_state].probability_7[idx_grid] = pow( cabs( eigenstate[idx_state].F_7[idx_grid] ), 2.0 );
				eigenstate[idx_state].probability_8[idx_grid] = pow( cabs( eigenstate[idx_state].F_8[idx_grid] ), 2.0 );
			
			}
	
		}
	
	}    
	
	for( int idx_state = 0; idx_state < 2*N_states[type]; idx_state++ )
	{

	       free( eigenstate_temp[idx_state].probability_1 );
	       free( eigenstate_temp[idx_state].probability_2 );
	       free( eigenstate_temp[idx_state].probability_3 );
	       free( eigenstate_temp[idx_state].probability_4 );
	       free( eigenstate_temp[idx_state].probability_5 );
	       free( eigenstate_temp[idx_state].probability_6 );
	       free( eigenstate_temp[idx_state].probability_7 );
	       free( eigenstate_temp[idx_state].probability_8 );
	       free( eigenstate_temp[idx_state].F_1 );
	       free( eigenstate_temp[idx_state].F_2 );
	       free( eigenstate_temp[idx_state].F_3 );
	       free( eigenstate_temp[idx_state].F_4 );
	       free( eigenstate_temp[idx_state].F_5 );
	       free( eigenstate_temp[idx_state].F_6 );
	       free( eigenstate_temp[idx_state].F_7 );
	       free( eigenstate_temp[idx_state].F_8 );

	}
	
	free( eigenstate_temp );

}
