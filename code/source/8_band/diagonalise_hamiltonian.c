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

	// Solver error indicator
	MKL_INT info;

	// Output verbose eigensolver info. 0(1) disable(enable).
	int verbose = 0;

	if (PACKED && nprocs == 1) // One MPI process and packed matrix storage
	{
		
		MKL_INT* ifail = (MKL_INT *)malloc(N_H * sizeof(MKL_INT)); // Error info
	
		char JOBVL = 'V'; // Computes eigenvalues and eigenvectors
		char UPLO = 'L'; // Matrix stores lower triangular part

		char range = 'I'; // The routine computes eigenvalues with indices il to iu
		double vl = 0; // If range = 'V', vl is the eigenvalue search lower bound
		double vr = 0; // If range = 'V', vr is the eigenvalue search upper bound
		MKL_INT il = idx_eig[0]; // Eigenvalue search indice lower bound
		MKL_INT iu = idx_eig[1]; // Eigenvalue search indice upper bound
		double abstol = dlamch("S"); // High accuracy setting for error tolerance
		MKL_INT matrix_layout = LAPACK_COL_MAJOR; // Column major matrix layout
		MKL_INT M = iu - il + 1; // Number of eigenvalues to search for
	
		info = LAPACKE_zhpevx(matrix_layout, JOBVL, range, UPLO, N_H, H, vl, vr, il, iu, abstol, &M, eigval, eigvec, N_H, ifail);
		
		free( ifail );
	
	}
	else if (nprocs == 1) // One MPI process and normal matrix storage
	{
		
		MKL_INT* ifail = (MKL_INT *)malloc(N_H * sizeof(MKL_INT));  // Error info
	
		char JOBVL = 'V'; // Computes eigenvalues and eigenvectors
		char UPLO = 'L'; // Matrix stores lower triangular part
		MKL_INT LDA = N_H; // Leading dimension of output array

		char range = 'I'; // The routine computes eigenvalues with indices il to iu
		double vl = 0; // If range = 'V', vl is the eigenvalue search lower bound
		double vr = 0; // If range = 'V', vr is the eigenvalue search upper bound
		MKL_INT il = idx_eig[0]; // Eigenvalue search indice lower bound
		MKL_INT iu = idx_eig[1]; // Eigenvalue search indice upper bound
		double abstol = dlamch("S"); // High accuracy setting for error tolerance
		MKL_INT matrix_layout = LAPACK_COL_MAJOR; // Column major matrix layout
		MKL_INT M = iu - il + 1; // Number of eigenvalues to search for
	
		info = LAPACKE_zheevr(matrix_layout, JOBVL, range, UPLO, N_H, H, LDA, vl, vr, il, iu, abstol, &M, eigval, eigvec, N_H, ifail);
		
		free( ifail );
	
	}
	else // Multiple MPI process ScaLAPACK eigensolver
	{
	
		// Employed eigensolver workspace size controlled through 'optimise' input parameter. Set to 1 to allow eigensolver to determine optimal solver workspace. WARNING: Results in large RAM usage by each process. Workspace memory usage in total is equal to roughly 4 times that used to store the Hamiltonian matrix (increases minimally with number of MPI processes). Set to 0 to assign the minimum workspace size required (usually negligible relative to memory used to store matrix elements). Optimal memory usage results in ~3 to 4 times faster calculation.
	
		if ( optimise == 0 && myid == 0)
		{

			printf( "\tUsing minimal eigensolver workspace size...\n" );

		}
		else if ( myid == 0 )
		{

			printf( "\tUsing optimal eigensolver workspace size...\n" );

		}
		
		char jobz = 'V'; // Computes eigenvalues and eigenvectors
		char range = 'I'; // The routine computes eigenvalues with indices il to iu
		char uplo = 'L'; // Matrix stores lower triangular part
		
		MKL_INT ia = 1; // First row in global matrix
		MKL_INT ja = 1; // First column in global matrix
		
		double vl = 0; // If range = 'V', vl is the eigenvalue search lower bound
		double vu = 0; // If range = 'V', vr is the eigenvalue search upper bound
		
		MKL_INT il = idx_eig[0]; // Eigenvalue search indice lower bound
		MKL_INT iu = idx_eig[1]; // Eigenvalue search indice upper bound
		
		MKL_INT M = iu - il + 1; // Number of eigenvalues to search for
		
		MKL_INT iz = 1; // First row in global eigenvector matrix
		MKL_INT jz = 1; // First column in global eigenvector matrix
		
		// Workspace parameters set to initially calculate optimal workspace
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
		
		MKL_INT *ifail = calloc(N_H, sizeof(MKL_INT)); // Error info
		
		MKL_INT *iclustr = calloc(2*np_rows*np_cols, sizeof(MKL_INT)); // Info on eigenvectors belonging to clustered eigenvalues that could not be reorthogonalised 
		
		double *gap = calloc(np_rows*np_cols, sizeof(double)); // Gap between clustered eigenvalues corresponding to eigenvectors that could not be reorthogonalised 
		

		// Initial run to determine optimal performance memory workspaces

		// Use of pzheevx in certain cases due to 2 seperate bugs in Intel MKL's pzheevr (27/05/2025):
		// 1. Inability to use memory workspace sizes less than the optimal despite being supported by Netlib's.
		// 2. Orthogonality of energy degenerate state eigenvectors not guarranteed in (what appears to be) cases where the number of degenrate eigenvalues is less than the number of MPI processes

		if ( !optimise || ( idx_eig[1] - idx_eig[0] + 1 ) / 2 < nprocs ) // In the case of non-optimal workspace size, or where the number of non-degenerate eigenvalues is less than the number of MPI processes
		{

			pzheevx( &jobz, &range, &uplo, &N_H, H, &ia, &ja, sc_desc, &vl, &vu, &il, &iu, &abstol, &M, &nz, eigval, &orfac, eigvec, &iz, &jz, sc_descz, &work_dum, &lwork, &rwork_dum, &lrwork, &iwork_dum, &liwork, ifail, iclustr, gap, &info );

		}
		else // Otherwise
		{

			pzheevr( &jobz, &range, &uplo, &N_H, H, &ia, &ja, sc_desc, &vl, &vu, &il, &iu, &M, &nz, eigval, eigvec, &iz, &jz, sc_descz, &work_dum, &lwork, &rwork_dum, &lrwork, &iwork_dum, &liwork, &info );

		}
		

		// Determine minimum memory workspace sizes in case of use of non-optimal solver. See Intel documentation for details.
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
		work = calloc(lwork, sizeof(MKL_Complex16));
		
		
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
		rwork = calloc(lrwork, sizeof(double));

		
		liwork = 6*N_H;
		if ( optimise )
		{
		
			liwork = iwork_dum;
			
		}
		iwork = calloc(liwork, sizeof(MKL_INT));

		if ( verbose )
		{
		
			printf("\nFrom process %d:\nlwork memory = %f MB\n rwork memory = %f MB\n iwork memory = %f MB\n\n", myid, 16 * lwork / (1024.0 * 1024), 8 * lrwork / (1024.0 * 1024), 8 * liwork / (1024.0 * 1024));

		}

		if ( verbose && myid == 0 )
		{

			mkl_verbose(1);

		}
		
		// Solve step
		if ( !optimise || ( idx_eig[1] - idx_eig[0] + 1 ) / 2 < nprocs ) // In the case of non-optimal workspace size, or where the number of non-degenerate eigenvalues is less than the number of MPI processes
		{

			pzheevx( &jobz, &range, &uplo, &N_H, H, &ia, &ja, sc_desc, &vl, &vu, &il, &iu, &abstol, &M, &nz, eigval, &orfac, eigvec, &iz, &jz, sc_descz, work, &lwork, rwork, &lrwork, iwork, &liwork, ifail, iclustr, gap, &info );

		}
		else // Otherwise
		{

			pzheevr( &jobz, &range, &uplo, &N_H, H, &ia, &ja, sc_desc, &vl, &vu, &il, &iu, &M, &nz, eigval, eigvec, &iz, &jz, sc_descz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info );

		}

		if ( verbose && myid == 0 )
		{

			mkl_verbose(0);

		}

		free( work );
		free( iwork );
		free( rwork );
		free( ifail );
		free( iclustr );
		free( gap );
	
	}
	
	
	if (info != 0)
	{
	
		if ( myid == 0 )
		{
		
			printf( "Eigensolver error. Info = %lld.\n\n", info );

		}

		MPI_Barrier( MPI_COMM_WORLD );

		MPI_Abort( MPI_COMM_WORLD, 1 );
	
	}
	
}
