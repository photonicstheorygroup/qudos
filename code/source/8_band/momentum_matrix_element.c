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

complex double momentum_matrix_element( int polarisation, int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct eigenstate eigenstate_1, MKL_INT eigenstate_1_idx, struct eigenstate eigenstate_2, MKL_INT eigenstate_2_idx, int *N_states )
{

	complex double result = 0.0 + I*0.0;


	// Compute dH/dk matrix for specified polarisation if not already done so
	if ( Polarisation_Indicator != polarisation )
	{
	
		hamiltonian_derivative_8_band( polarisation, N_shape, basis, supercell, parameters, dH_dk );
		Polarisation_Indicator = polarisation;
		
	}


	// Compute momentum matrix element

	char UPLO = 'L'; // Matrix stores lower triangular part
	MKL_Complex16 alpha = 1; // Matrix prefactor
	MKL_Complex16 beta = 0; // Vector (1 column matrix) prefactor
	MKL_INT incx = 1; // Increment of input vector
	MKL_INT incy = 1; // Increment of output vector
	MKL_INT ia = 1; // First row of matrix
	MKL_INT ja = 1; // First column of matrix
	MKL_INT ix = 1; // First row of input vectors
	MKL_INT jx = 1; // First column of input vectors
	MKL_INT iy = 1; // First row of output vector
	MKL_INT jy = 1; // First column of output vector
	
	MKL_INT sc_descv1[9]; // Array descriptor of input vectors
	MKL_INT sc_descv2[9]; // Array descriptor of output vector
	MKL_INT N_COL = 1; // 1 column in matrix (array) for vector
	MKL_INT lrsrc = 0; // Process row over which first row of matrix containing input vectors is distributed
	MKL_INT lcsrc = ( (eigenstate_1_idx) / block_size ) % np_cols; // Process column over which first column of matrix containing input vector 1 is distributed
	
	// Initialise the array descriptor for input vector 1
	DESCINIT(sc_descv1, &N_H, &N_COL, &block_size, &block_size, &lrsrc, &lcsrc, &my_blacs_ctxt, &na_rowsz, &info);
	
	lcsrc = ( (eigenstate_2_idx) / block_size ) % np_cols; // Process column over which first column of matrix containing input vector 2 is distributed
	
	// Initialise the array descriptor for input vector 2
	DESCINIT(sc_descv2, &N_H, &N_COL, &block_size, &block_size, &lrsrc, &lcsrc, &my_blacs_ctxt, &na_rowsz, &info);
	
	// Matrix-vector product
	if ( eigenstate_2_idx == 2*N_states[0] - 1 ) // If first valence band state, require matrix-vector product. Otherwise, matrix-vector product has been stored from calculation involving previous valence band state. 
	{
	
		if ( PACKED && nprocs == 1) // Single MPI process and using PACKED matrix storage
		{
		
			cblas_zhpmv( CblasColMajor, CblasLower, N_H, &alpha, dH_dk, eigenstate_1.fourier, incx, &beta, yV, incy );
		
		}
		else if ( nprocs == 1 ) // Single MPI process and using normal matrix storage
		{
		
			cblas_zhemv( CblasColMajor, CblasLower, N_H, &alpha, dH_dk, N_H, eigenstate_1.fourier, incx, &beta, yV, incy );
		
		}
		else // Multiple MPI processes
		{
		
			pzhemv( &UPLO, &N_H, &alpha, dH_dk, &ia, &ja, sc_desc, eigenstate_1.fourier, &ix, &jx, sc_descv1, &incx, &beta, yV, &iy, &jy, sc_descv1, &incy );
		
		}
	
	}
	
	// Vector-vector dot product
	if ( nprocs == 1 ) // Single MPI process
	{
	
		cblas_zdotc_sub( N_H, eigenstate_2.fourier, incx, yV, incy, &result );
	
	}
	else // Multiple MPI processes
	{
	
		pzdotc( &N_H, &result, eigenstate_2.fourier, &ix, &jx, sc_descv2, &incx, yV, &iy, &jy, sc_descv1, &incy );
	
	}
	
	// Broadcast result to all processes (processes not involved in pzdotc may not have result otherwise)
	MPI_Bcast( &result, 1, MPI_DOUBLE_COMPLEX, lcsrc * np_rows, MPI_COMM_WORLD );

	return result;

} 
