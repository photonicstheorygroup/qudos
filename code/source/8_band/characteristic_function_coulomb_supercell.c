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

void characteristic_function_coulomb_supercell( struct wave_vectors coulomb_basis, struct supercell_geometry *supercell, struct characteristic_functions characteristic_function_coulomb_G )
{
	
	// Local size of arrays
	int local_size = coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T / nprocs + ( ( myid < coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T % nprocs ) ? 1 : 0 );
	
	// Temporary array to store characteristic functions calculated by the process
	complex double *char_func = ( complex double * ) malloc( local_size*sizeof( complex double ) );
	
	// Will store number of characteristic functions calculated by each process
	int *recvcounts = ( int * ) malloc( nprocs*sizeof( int ) );
	
	// Will store displacement from first value of the "global" characteristic function for the first value of the characteristic functions calculated by each process
	int *displs = ( int * ) malloc( nprocs*sizeof( int ) );
	
	for ( int proc = 0; proc < nprocs; proc++ )
	{
	
		recvcounts[proc] = coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T / nprocs + ( ( proc < coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T % nprocs ) ? 1 : 0 );
		displs[proc] = ( proc == 0 ? 0 : displs[proc - 1] + recvcounts[proc - 1] ); 
	
	}
	
	// Loop over wave vectors and evaluate the Fourier transforms of the shape characteristic functions
	#pragma omp parallel for
	for( int idx = 0; idx < local_size; idx++ )
	{

		int idx_coulomb_basis = myid * ( coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T / nprocs ) + ( ( myid < coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T % nprocs ) ? myid : coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T % nprocs ) + idx;
		
		int idx_Gx = idx_coulomb_basis / ( coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T );
		
		int idx_Gy = ( idx_coulomb_basis % ( coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T ) ) / coulomb_basis.N_Gz_T;
		
		int idx_Gz = ( idx_coulomb_basis % ( coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T ) ) % coulomb_basis.N_Gz_T;
		
		char_func[idx] = characteristic_function_fourier( coulomb_basis.Gx[idx_Gx], coulomb_basis.Gy[idx_Gy], coulomb_basis.Gz[idx_Gz], supercell, 0, 3 );
		
	}
	
	
	MPI_Allgatherv( char_func, local_size, MPI_DOUBLE_COMPLEX, characteristic_function_coulomb_G.chi, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
		
	
	free( char_func );
	free( recvcounts );
	free( displs );

}
