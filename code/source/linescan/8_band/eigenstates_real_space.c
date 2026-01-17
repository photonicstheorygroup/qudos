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

#include "8_band_linescan.h"

void eigenstates_real_space( int idx_grid, int idx_shape, double x, double y, double z, double Fx, double Fy, double Fz, struct material_params parameters, struct strain_tensor strain_r, complex double *piezo_r, struct band_edge_energies band_edges )
{

	MKL_INT m, N_H = N_BAND, il = 0, iu = 0; // Number of eigenvalues found by eigensolver; Order of Hamiltonian matrix; lower and upper indices of found eigenvalues
	char jobz = 'V', uplo = 'U', range = 'A'; // Setting eigensolver to find eigenvalues and eigenvectors; Upper triangular part of Hamiltonian is stored; Setting eigensolver to compute all eigenvalues
	double *eigval; // Stores computed eigenvalues
	double vl = 0, vu = 0; // Unreferenced by eigensolver
	int matrix_layout = LAPACK_COL_MAJOR; // Column major Hamiltonian matrix layout

	double abstol = 2*dlamch("S"); // High accuracy setting for error tolerance

	MKL_INT *isuppz; // The support of the eigenvectors (see Intel developer reference for LAPACKE_zheevr)


	// Extract strain tensor at current real space grid point

	struct strain_tensor epsilon;

	epsilon.xx = ( complex double* ) malloc( sizeof( complex double ) );
	epsilon.yy = ( complex double* ) malloc( sizeof( complex double ) );
	epsilon.zz = ( complex double* ) malloc( sizeof( complex double ) );
	epsilon.yz = ( complex double* ) malloc( sizeof( complex double ) );
	epsilon.xz = ( complex double* ) malloc( sizeof( complex double ) );
	epsilon.xy = ( complex double* ) malloc( sizeof( complex double ) );

	epsilon.xx[0] = strain_r.xx[idx_grid];
	epsilon.yy[0] = strain_r.yy[idx_grid];
	epsilon.zz[0] = strain_r.zz[idx_grid];
	epsilon.yz[0] = strain_r.yz[idx_grid];
	epsilon.xz[0] = strain_r.xz[idx_grid];
	epsilon.xy[0] = strain_r.xy[idx_grid];


	// Extract piezoelectric potential at current real space grid point

	complex double phi_piezo;

	phi_piezo = piezo_r[idx_grid];


	// Compute potential due to applied electric field at current real space grid point

	complex double phi_F;

	phi_F = 1.0*( Fx*x + Fy*y + Fz*z );


	// Set up bulk 8-band k.p Hamiltonian at the current real space grid point

	complex double *H = ( complex double* ) malloc( N_BAND*N_BAND*sizeof( complex double ) );
	complex double *eigvec = ( complex double* ) malloc( N_BAND*N_BAND*sizeof( complex double ) );

	for( int idx_row = 0; idx_row < N_BAND; idx_row++ )
	{

		for( int idx_col = 0; idx_col < N_BAND; idx_col++ )
		{

			H[idx_row + N_BAND*idx_col] = hamiltonian_8_band_bulk( idx_row, idx_col, x, y, z, parameters, epsilon, phi_piezo, phi_F );

		}

	}

	free( epsilon.xx );
	free( epsilon.yy );
	free( epsilon.zz );
	free( epsilon.yz );
	free( epsilon.xz );
	free( epsilon.xy );


	isuppz =( MKL_INT* ) malloc( 2*N_H*sizeof( MKL_INT ));


	// Diagonalise bulk 8-band k.p Hamiltonian

	eigval = ( double* ) malloc(     N_H      *sizeof( double ) );

	LAPACKE_zheevr( matrix_layout, jobz, range, uplo, N_H, H, N_H, vl, vu, il, iu, abstol, &m, eigval, eigvec, N_H, isuppz );
	
	// Identifying band of band-edge energies
	for (int i = 0; i < N_BAND / 2; i++){
		int max = 0;

		for (int j = 1; j < N_BAND / 2; j++){ // Finding largest eigenvector for each eigenvalue
			if ( cabs(eigvec[(2*i*N_BAND)+j]) + cabs(eigvec[(2*i*N_BAND)+j+4]) + cabs(eigvec[((2*i+1)*N_BAND)+j]) + cabs(eigvec[((2*i+1)*N_BAND)+j+4]) > cabs(eigvec[(2*i*N_BAND)+max]) + cabs(eigvec[(2*i*N_BAND)+max+4]) + cabs(eigvec[((2*i+1)*N_BAND)+max]) + cabs(eigvec[((2*i+1)*N_BAND)+max+4]))
				max = j;
		}

		if (max == 3)
			band_edges.E_so[idx_grid] = eigval[2*i];
		else if (max == 2)
			band_edges.E_lh[idx_grid] = eigval[2*i];
		else if (max == 1)
			band_edges.E_hh[idx_grid] = eigval[2*i];
		else if (max == 0)
			band_edges.E_cb[idx_grid] = eigval[2*i];

	}
	

	free( H );
	free( eigval );
	free( eigvec ); 
	free( isuppz );

}
