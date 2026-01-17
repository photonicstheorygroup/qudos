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

complex double spatial_overlap( struct wave_vectors basis, struct eigenstate eigenstate_1, MKL_INT eigenstate_1_idx, struct eigenstate eigenstate_2, MKL_INT eigenstate_2_idx )
{

	complex double result = 0.0 + I*0.0;

	// COmpute overlap between two states
	
	MKL_INT sc_descv1[9]; // Array descriptor of 1st input vector (1 column matrix)
	MKL_INT sc_descv2[9]; // Array descriptor of 2nd input vector (1 column matrix)
	MKL_INT N_COL = 1; // Number of columns in input matrices (1 as they are vectors)
	MKL_INT lrsrc = 0; // Process row over which first row of each matrix containing the input vector is distributed
	MKL_INT lcsrc = ( ( eigenstate_1_idx ) / block_size ) % np_cols; // Process column over which first column of matrix containing input vector 1 is distributed
	MKL_INT ix = 1, jx = 1, iy = 1, jy = 1; // 1st row and column of input vectors (1 column matrices)
	MKL_INT incx = 1, incy = 1;
	
	// Initialise the array descriptor for input vector 1
	DESCINIT( sc_descv1, &N_H, &N_COL, &block_size, &block_size, &lrsrc, &lcsrc, &my_blacs_ctxt, &na_rowsz, &info );
	
	// Process column over which first column of matrix containing input vector 2 is distributed
	lcsrc = ( ( eigenstate_2_idx ) / block_size ) % np_cols;
	
	// Initialise the array descriptor for input vector 2
	DESCINIT( sc_descv2, &N_H, &N_COL, &block_size, &block_size, &lrsrc, &lcsrc, &my_blacs_ctxt, &na_rowsz, &info );
	
	// Vector-vector dot product
	if ( nprocs == 1 ) // Single MPI process
	{
	
		cblas_zdotc_sub( N_H, eigenstate_1.fourier, incx, eigenstate_2.fourier, incy, &result );
	
	}
	else // Multiple MPI processes
	{
	
		pzdotc( &N_H, &result, eigenstate_1.fourier, &ix, &jx, sc_descv1, &incx, eigenstate_2.fourier, &iy, &jy, sc_descv2, &incy );
		
	}
	
	// Broadcast result to all processes (processes not involved in pzdotc may not have result otherwise)
	MPI_Bcast( &result, 1, MPI_DOUBLE_COMPLEX, lcsrc * np_rows, MPI_COMM_WORLD );

	return result;

}
