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

void characteristic_functions_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, int strain_region )
{
	
	// Local size of arrays
	int local_size = basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs + ( ( myid < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? 1 : 0 );
	
	// Temporary array to store characteristic functions calculated by the process
	complex double *char_func = ( complex double * ) malloc( local_size*sizeof( complex double ) );
	
	// Will store number of characteristic functions calculated by each process
	int *recvcounts = ( int * ) malloc( nprocs*sizeof( int ) );
	
	// Will store displacement from first value of the "global" characteristic function for the first value of the characteristic functions calculated by each process
	int *displs = ( int * ) malloc( nprocs*sizeof( int ) );
	
	for ( int proc = 0; proc < nprocs; proc++ )
	{
	
		recvcounts[proc] = basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs + ( ( proc < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? 1 : 0 );
		displs[proc] = ( proc == 0 ? 0 : displs[proc - 1] + recvcounts[proc - 1] ); 
	
	}

	int shape_index_start;
	if ( strain_region == 0 || strain_region == 1 ) // Electronic structure supercell or strain supercell with free lengths indicator
	{

		shape_index_start = 1;

	}
	else // Otherwise strain supercell with lengths commensurate with electronic structure supercell lengths
	{

		shape_index_start = 0;

	}

	// Loop over shapes
	for( int idx_shape = shape_index_start; idx_shape < N_shape; idx_shape++ )
	{
	
		// Loop over local set of wave vectors and evaluate the Fourier transforms of the shape characteristic functions
		#pragma omp parallel for
		for( int idx = 0; idx < local_size; idx++ )
		{
	
			int idx_basis = myid * ( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs ) + ( ( myid < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? myid : basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) + idx;
			
			int idx_Gx = idx_basis / ( basis.N_Gy_T*basis.N_Gz_T );
			
			int idx_Gy = ( idx_basis % ( basis.N_Gy_T*basis.N_Gz_T ) ) / basis.N_Gz_T;
			
			int idx_Gz = ( idx_basis % ( basis.N_Gy_T*basis.N_Gz_T ) ) % basis.N_Gz_T;
			
			char_func[idx] = characteristic_function_fourier( basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, idx_shape, strain_region );
			
			
		}
		
		// Gather to all processes
		if ( strain_region == 0 ) // If electronic structure supercell
		{
		
			MPI_Allgatherv( char_func, local_size, MPI_DOUBLE_COMPLEX, characteristic_function_G[idx_shape].chi, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
			
		}
		else // If strain supercell
		{
		
			MPI_Allgatherv( char_func, local_size, MPI_DOUBLE_COMPLEX, characteristic_function_strain_G[idx_shape].chi, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
			
		}
						
	}
	
	free( char_func );
	free( recvcounts );
	free( displs );

}
