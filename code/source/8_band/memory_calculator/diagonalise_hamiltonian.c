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

void diagonalise_hamiltonian( int *idx_eig, complex double *H, double *eigval, complex double *eigvec )
{

	MKL_INT info;

	int verbose = 0;

	if (PACKED && nprocs == 1)
	{
		MKL_INT* ifail = (MKL_INT *)malloc(N_H * sizeof(MKL_INT));
	
		char JOBVL = 'V'; // Computes eigenvalues and eigenvectors
		char UPLO = 'L'; // Matrix stores lower triangular part
		MKL_INT LDA = N_H; // Leading dimention of output array

		char range = 'I'; // The routine computes eigenvalues with indices il to iu
		double vl = 0; // If range = 'V', vl is the eigenvalue search lower bound
		double vr = 0; // If range = 'V', vr is the eigenvalue search upper bound
		MKL_INT il = idx_eig[0]; // Eigenvalue search indice lower bound
		MKL_INT iu = idx_eig[1]; // Eigenvalue search indice upper bound
		double abstol = dlamch("S"); // High accuracy setting for error tolerance
		MKL_INT matrix_layout = LAPACK_COL_MAJOR; // Column major matrix layout
		MKL_INT M = iu - il + 1; // Number of eigenvalues to search for
	
		//info = LAPACKE_zhpevx(matrix_layout, JOBVL, range, UPLO, N_H, H, vl, vr, il, iu, abstol, &M, eigval, eigvec, N_H, ifail);
		
		free( ifail );
	
	}
	
	else if (nprocs == 1)
	{
		
		MKL_INT* ifail = (MKL_INT *)malloc(N_H * sizeof(MKL_INT));
	
		char JOBVL = 'V'; // Computes eigenvalues and eigenvectors
		char UPLO = 'L'; // Matrix stores lower triangular part
		MKL_INT LDA = N_H; // Leading dimention of output array

		char range = 'I'; // The routine computes eigenvalues with indices il to iu
		double vl = 0; // If range = 'V', vl is the eigenvalue search lower bound
		double vr = 0; // If range = 'V', vr is the eigenvalue search upper bound
		MKL_INT il = idx_eig[0]; // Eigenvalue search indice lower bound
		MKL_INT iu = idx_eig[1]; // Eigenvalue search indice upper bound
		double abstol = dlamch("S"); // High accuracy setting for error tolerance
		MKL_INT matrix_layout = LAPACK_COL_MAJOR; // Column major matrix layout
		MKL_INT M = iu - il + 1; // Number of eigenvalues to search for
	
		//info = LAPACKE_zheevr(matrix_layout, JOBVL, range, UPLO, N_H, H, LDA, vl, vr, il, iu, abstol, &M, eigval, eigvec, N_H, ifail);
		
		free( ifail );
	
	}
	
	else
	{
	
		char jobz = 'V'; // Computes eigenvalues and eigenvectors
		char range = 'I'; // The routine computes eigenvalues with indices il to iu
		char uplo = 'L'; // Matrix stores lower triangular part
		
		MKL_INT ia = 1;
		MKL_INT ja = 1;
		
		double vl = 0; // If range = 'V', vl is the eigenvalue search lower bound
		double vu = 0; // If range = 'V', vr is the eigenvalue search upper bound
		
		MKL_INT il = idx_eig[0]; // Eigenvalue search indice lower bound
		MKL_INT iu = idx_eig[1]; // Eigenvalue search indice upper bound
		
		MKL_INT M = iu - il + 1; // Number of eigenvalues to search for
		
		MKL_INT iz = 1;
		MKL_INT jz = 1;
		
		MKL_INT lwork = -1;
		MKL_Complex16 work_dum;
		MKL_Complex16 *work;
		
		MKL_INT liwork = -1;
		MKL_INT iwork_dum;
		MKL_INT *iwork;
		
		MKL_INT lrwork = -1;
		double rwork_dum;
		double *rwork;
		
		MKL_INT nz;
		
		double abstol = 2*pdlamch(&my_blacs_ctxt, "S"); // High accuracy setting for error tolerance
		
		double orfac = 1e-3; // Default value for ensuring orthogonality of eigenvectors for pzheevx
		
		MKL_INT *ifail = calloc(N_H, sizeof(MKL_INT));

		local_mem += ( N_H ) *sizeof( MKL_INT );	// Adding memory for ifail
		
		MKL_INT *iclustr = calloc(2*np_rows*np_cols, sizeof(MKL_INT));
		
		double *gap = calloc(np_rows*np_cols, sizeof(double));


		// Initial run to determine optimal performance memory workspaces

		// Use of pzheevx in certain cases due to 2 seperate bugs in Intel MKL's pzheevr (27/05/2025):
		// 1. Inability to use memory workspace sizes less than the optimal despite being supported by Netlib's.
		// 2. Orthogonality of energy degenerate state eigenvectors not guarranteed in (what appears to be) cases where the number of degenrate eigenvalues is less than the number of MPI processes
	
		if ( !optimise || ( idx_eig[1] - idx_eig[0] + 1 ) / 2 < nprocs )
		{

			pzheevx( &jobz, &range, &uplo, &N_H, H, &ia, &ja, sc_desc, &vl, &vu, &il, &iu, &abstol, &M, &nz, eigval, &orfac, eigvec, &iz, &jz, sc_descz, &work_dum, &lwork, &rwork_dum, &lrwork, &iwork_dum, &liwork, ifail, iclustr, gap, &info );

		}
		else
		{

			pzheevr( &jobz, &range, &uplo, &N_H, H, &ia, &ja, sc_desc, &vl, &vu, &il, &iu, &M, &nz, eigval, eigvec, &iz, &jz, sc_descz, &work_dum, &lwork, &rwork_dum, &lrwork, &iwork_dum, &liwork, &info );

		}
		
		
		// Determining minimum memory workspace sizes in case of use of non-optimal solver
		
		MKL_INT iproc = 0;
		MKL_INT srcproc = 0;

		MKL_INT np0 = NUMROC(&N_H, &block_size, &iproc, &srcproc, &np_rows);

		MKL_INT nq0 = NUMROC(&N_H, &block_size, &iproc, &srcproc, &np_cols);

		MKL_INT np00 = NUMROC(&N_H, &block_size, &iproc, &srcproc, &np_rows);
		
		MKL_INT mq00;
		if (M > block_size)
		{
		
			mq00 = NUMROC(&M, &block_size, &iproc, &srcproc, &np_cols);
		
		}
		else
		{
		
			mq00 = NUMROC(&block_size, &block_size, &iproc, &srcproc, &np_cols);
		
		}

		lwork = N_H + ( np0 + nq0 + block_size )*block_size;


		if ( optimise )
		{
		
			lwork = (MKL_INT)(creal(work_dum) + 0.5);
			
		}
		//work = calloc(lwork, sizeof(MKL_Complex16));

		local_mem += ( lwork ) *sizeof( complex double );	// Adding memory for lwork

		long long eig_work_mem = ( lwork ) *sizeof( complex double );
		
		
		if (5 * N_H > np00 * mq00 )
		{
		
			lrwork = 5*N_H + 4*N_H + N_H*( ((M + np_rows*np_cols - 1) / (np_rows*np_cols))) + (M - 1)*N_H;
		
		}
		else
		{
		
			lrwork = np00 * mq00 + 4*N_H + N_H*( ((M + np_rows*np_cols - 1) / (np_rows*np_cols)) ) + (M - 1)*N_H;
		
		}
		
		if ( optimise && ( idx_eig[1] - idx_eig[0] + 1 ) / 2 < nprocs )
		{
		
			lrwork = (MKL_INT)(rwork_dum + 0.5) + (M - 1)*N_H;
			
		}
		else if ( optimise )
		{
		
			lrwork = (MKL_INT)(rwork_dum + 0.5);
			
		}
		//rwork = calloc(lrwork, sizeof(double));

		local_mem += ( lrwork ) *sizeof( double );	// Adding memory for rwork
		eig_work_mem += ( lrwork ) *sizeof( double );
		
		liwork = 6*N_H;
		if ( optimise )
		{
		
			liwork = iwork_dum;
			
		}
		//iwork = calloc(liwork, sizeof(MKL_INT));

		local_mem += ( liwork ) *sizeof( MKL_INT );	// Adding memory for iwork
		eig_work_mem += ( liwork ) *sizeof( MKL_INT );

		local_mem_peak = local_mem;
		component_memory_usage.eig_work_mem = eig_work_mem;

		if ( verbose )
		{

			printf("\nFrom process %d:\nlwork memory = %f MB\n rwork memory = %f MB\n iwork memory = %f MB\n\n", myid, 16 * lwork / (1024.0 * 1024), 8 * lrwork / (1024.0 * 1024), 8 * liwork / (1024.0 * 1024));

		}

		
		//pzheevx( &jobz, &range, &uplo, &N_H, H, &ia, &ja, sc_desc, &vl, &vu, &il, &iu, &abstol, &M, &nz, eigval, &orfac, eigvec, &iz, &jz, sc_descz, work, &lwork, rwork, &lrwork, iwork, &liwork, ifail, iclustr, gap, &info );
		
		//pzheevr( &jobz, &range, &uplo, &N_H, H, &ia, &ja, sc_desc, &vl, &vu, &il, &iu, &M, &nz, eigval, eigvec, &iz, &jz, sc_descz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info );


		//free( work );
		//free( iwork );
		//free( rwork );
		free( ifail );
		free( iclustr );
		free( gap );

		local_mem -= ( lwork) *sizeof( complex double );
		local_mem -= ( lrwork) *sizeof( double );
		local_mem -= ( liwork) *sizeof( MKL_INT );
		local_mem -= ( N_H ) *sizeof( MKL_INT );
	
	}
	
	
}
