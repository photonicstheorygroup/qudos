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

void eigenstates( int N_shape, int *N_states, struct vector F, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb )
{

	// Indices of eigenstates to be computed

	int idx_vb_min, idx_vb_max, idx_cb_min, idx_cb_max, idx_eig[2];

	idx_cb_min = ( N_BAND - 2 )*basis.N_Gx*basis.N_Gy*basis.N_Gz; // Index of  lowest energy conduction band state
	idx_cb_max = idx_cb_min + 2*N_states[1] - 1;                  // Index of highest energy conduction band state (factor of 2 to account for spin degeneracy)

	idx_vb_max = idx_cb_min - 1;                                  // Index of highest energy valence band state
	idx_vb_min = idx_vb_max - 2*N_states[0] + 1;                  // Index of  lowest energy valence band state (factor of 2 to account for spin degeneracy)

	idx_eig[0] = idx_vb_min + 1; // Index of  lowest energy eigenstate to be calculated (+1 for LAPACK Fortran array indexing)
	idx_eig[1] = idx_cb_max + 1; // Index of highest energy eigenstate to be calculated (+1 for LAPACK Fortran array indexing)

	        double *eigval;// = (         double* ) malloc( N_H    *sizeof(         double ) ); // Array to store Hamiltonian eigenvalues
	complex double *eigvec;// = ( complex double* ) malloc( na_rowsz * na_colsz *sizeof( complex double ) ); // Array to store Hamiltonian eigenvectors

	local_mem += ( N_H ) *sizeof( double ); // Adding mmeory to store eigenvalues

	local_mem += ( na_rowsz * na_colsz ) *sizeof( complex double );	// Adding memory to local components of eigenvector matrix
	


	// Construct and diagonalise supercell Hamiltonian

	complex double *H;
	
	if (PACKED && nprocs == 1)
	{
	
		//H = ( complex double* ) calloc( N_H*(N_H + 1)/2, sizeof( complex double ) );
		local_mem += ( N_H*(N_H + 1)/2 ) *sizeof( complex double ); // Adding memory to local components of Hamiltonian matrix
	
	}
	
	else
	{
	
		//H = ( complex double* ) calloc( na_rows*na_cols, sizeof( complex double ) );

		local_mem += ( na_rows * na_cols ) *sizeof( complex double ); // Adding memory to local components of Hamiltonian matrix
		
	}
    	

	diagonalise_hamiltonian( idx_eig, H, eigval, eigvec );


	local_mem -= ( na_rowsz * na_colsz ) *sizeof( complex double );
	local_mem -= ( N_H ) *sizeof( double );
	
	if (PACKED && nprocs == 1)
	{
	
		local_mem -= ( N_H*(N_H + 1)/2 ) *sizeof( complex double );
	
	}
	
	else
	{
	
		local_mem -= ( na_rows * na_cols ) *sizeof( complex double );
		
	}

}
