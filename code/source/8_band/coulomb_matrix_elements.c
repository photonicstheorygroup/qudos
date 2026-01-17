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

void coulomb_matrix_elements( struct wave_vectors basis, struct wave_vectors coulomb_basis, struct supercell_geometry *supercell, struct material_params *parameters, int *N_states, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb, struct coulomb_matrix_element CME, int coulomb_region_type )
{

	// Global number of elements in each B array after upsampling to Coulomb supercell 
	int N_GC_conv = ( 2*coulomb_basis.N_Gx - 1 )*( 2*coulomb_basis.N_Gy - 1 )*( 2*coulomb_basis.N_Gz - 1 );

	
	/////////////////////////////////////////////////////////////////
	//// Calculating B functions to be used in Vijkl calculation ////
	/////////////////////////////////////////////////////////////////
	
	// Process limit for Dfti FFT library
	MKL_INT process_limit = ( 2*basis.N_Gx - 1 <= 2*basis.N_Gy - 1 ? 2*basis.N_Gx - 1 : 2*basis.N_Gy - 1 );
	// Number of processes for Dfti library
	MKL_INT process_num;

	// Sub-communicator splits communicator into 2 for process IDs under limit and above limit
	MPI_Comm MPI_SUB_COMM;
	if ( nprocs > process_limit ) // If number of MPI processes greater than process limit
	{

		MPI_Comm_split( MPI_COMM_WORLD, ( myid < process_limit ? 0 : 1 ), myid, &MPI_SUB_COMM ); // Split communicator
		process_num = process_limit;

	}
	else // If number of MPI processes less than or equal to process limit
	{

		MPI_Comm_dup( MPI_COMM_WORLD, &MPI_SUB_COMM ); // Duplicate communicator
		process_num = nprocs;

	}

	// Will store information on how B arrays should be distributed among processes for Vijkl calculation
	int *sendcounts = ( int * ) malloc( nprocs*sizeof( int ) );
	int *displs_send = ( int * ) malloc( nprocs*sizeof( int ) );

	for ( int proc = 0; proc < nprocs; proc++ )
	{

		sendcounts[proc] = N_GC_conv / nprocs + ( ( proc < N_GC_conv % nprocs ) ? 1 : 0 );

		displs_send[proc] = ( proc == 0 ? 0 : displs_send[proc - 1] + sendcounts[proc - 1] );
	
	}


	// Number of unqiue B functions due to symmetry property
	MKL_INT N_state_conv = ( ( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] + 1 ) ) / 2 );

	// Stores B functions
	struct convolution_state *convolution_states = malloc( N_state_conv*sizeof( _convolution_state ) );
	// Stores wave vector reversed version of B functions
	struct convolution_state *convolution_states_reversed = malloc( N_state_conv*sizeof( _convolution_state ) );

	// Number of process local B wave vector components
	int local_size_conv = N_GC_conv / nprocs + ( ( myid < N_GC_conv % nprocs ) ? 1 : 0 );

	double start = MPI_Wtime();

	if ( myid < process_limit ) // Processes with IDs less than the process limit
	{

		// Total number of elements in each dimension for zero-padded arrays involved in Dfti FFT calculations as part of convolution calculation
		MKL_LONG sizes[] = { 2*basis.N_Gx - 1, 2*basis.N_Gy - 1, 2*basis.N_Gz - 1 };

		// Global size of involved arrays
		MKL_LONG global_size = sizes[0]*sizes[1]*sizes[2];

		// Dfti precision & array type
		MKL_LONG DFTI_PRECISION = DFTI_DOUBLE;
		MKL_LONG DFTI_TYPE = DFTI_COMPLEX;

		// Dfti task descriptor
		DFTI_DESCRIPTOR_DM_HANDLE task;

		// Initialising Dfti task descriptor
		check_dfti( DftiCreateDescriptorDM( MPI_SUB_COMM, &task, DFTI_PRECISION, DFTI_TYPE, 3, sizes ), "DftiCreateDescriptorDM" );

		// Set Dfti backward scale to represent inverse FFT
		check_dfti( DftiSetValueDM( task, DFTI_BACKWARD_SCALE, 1.0 / global_size ), "DftiSetValueDM" );

		// Local memory size for arrays involved in Dfti FFT calculations
		MKL_LONG local_memory_size;
		check_dfti( DftiGetValueDM( task, CDFT_LOCAL_SIZE, &local_memory_size ), "DftiGetValueDM" );
		// Global position of first element of local arrays involved in Dfti FFT calculations
		MKL_LONG local_start;
		check_dfti( DftiGetValueDM( task, CDFT_LOCAL_X_START, &local_start ), "DftiGetValueDM" );
		// Number of global rows in local arrays involved in Dfti FFT calculations
		MKL_LONG local_rows;
		check_dfti( DftiGetValueDM( task, CDFT_LOCAL_NX, &local_rows ), "DftiGetValueDM" );
		// Number of elements in local arrays involved in Dfti FFT calculations
		MKL_LONG local_size = local_rows*sizes[1]*sizes[2];

		// Committing Dfti descriptor
		check_dfti( DftiCommitDescriptorDM( task ), "DftiCommitDescriptorDM" );


		// Will hold local parts of padded envelope functions involved in convolution calculation
		struct eigenstate *eigenstate_vb_padded = malloc( 2*N_states[0]*sizeof( _eigenstate ) );
		struct eigenstate *eigenstate_cb_padded = malloc( 2*N_states[1]*sizeof( _eigenstate ) );
		struct eigenstate *eigenstate_vb_padded_conj = malloc( 2*N_states[0]*sizeof( _eigenstate ) );
		struct eigenstate *eigenstate_cb_padded_conj = malloc( 2*N_states[1]*sizeof( _eigenstate ) );

		for ( int idx_state_vb = 0; idx_state_vb < 2*N_states[0]; idx_state_vb++ )
		{
		
			eigenstate_vb_padded[idx_state_vb].fourier   = calloc( N_BAND*local_memory_size, sizeof( complex double ) );
			eigenstate_vb_padded_conj[idx_state_vb].fourier   = calloc( N_BAND*local_memory_size, sizeof( complex double ) );

		}
		for ( int idx_state_cb = 0; idx_state_cb < 2*N_states[1]; idx_state_cb++ )
		{
		
			eigenstate_cb_padded[idx_state_cb].fourier   = calloc( N_BAND*local_memory_size, sizeof( complex double ) );
			eigenstate_cb_padded_conj[idx_state_cb].fourier   = calloc( N_BAND*local_memory_size, sizeof( complex double ) );

		}

		// Redistributing the block-cyclically distributed envelope functions that were output from the eigensolver so that they are distributed among processes over the first spatial dimension, enabling use of the Dfti library
		state_redistribution( basis, supercell, N_states, eigenstate_vb, eigenstate_cb, eigenstate_vb_padded, eigenstate_cb_padded, eigenstate_vb_padded_conj, eigenstate_cb_padded_conj, process_limit, sizes, local_memory_size, local_start, local_rows );
		
		for ( int idx_state_vb = 0; idx_state_vb < 2*N_states[0]; idx_state_vb++ ) // Iterating over VB states
		{

			for ( int idx_band = 0; idx_band < N_BAND; idx_band++ ) // Iterating over bands
			{

				// Computing FFT using Dfti library of envelope function band components
				check_dfti( DftiComputeForwardDM( task, &eigenstate_vb_padded[idx_state_vb].fourier[idx_band*local_memory_size] ), "DftiComputeForwardDM" );
				check_dfti( DftiComputeForwardDM( task, &eigenstate_vb_padded_conj[idx_state_vb].fourier[idx_band*local_memory_size] ), "DftiComputeForwardDM" );

			}

		}
		for ( int idx_state_cb = 0; idx_state_cb < 2*N_states[1]; idx_state_cb++ ) // Iterating over CB states
		{

			for ( int idx_band = 0; idx_band < N_BAND; idx_band++ ) // Iterating over bands
			{

				// Computing FFT using Dfti library of envelope function band components
				check_dfti( DftiComputeForwardDM( task, &eigenstate_cb_padded[idx_state_cb].fourier[idx_band*local_memory_size] ), "DftiComputeForwardDM" );
				check_dfti( DftiComputeForwardDM( task, &eigenstate_cb_padded_conj[idx_state_cb].fourier[idx_band*local_memory_size] ), "DftiComputeForwardDM" );

			}

		}


		for ( int idx_state = 0; idx_state < N_state_conv; idx_state++ )
		{

			convolution_states[idx_state].conv_state =  calloc( local_size_conv  , sizeof( complex double ) );
			convolution_states_reversed[idx_state].conv_state =  calloc( local_size_conv  , sizeof( complex double ) );

		}

		// Calculating B arrays
		B_convolutions( basis, coulomb_basis, supercell, N_states, eigenstate_vb_padded, eigenstate_cb_padded, eigenstate_vb_padded_conj, eigenstate_cb_padded_conj, process_num, process_limit, task, local_memory_size, local_size, convolution_states, convolution_states_reversed, sendcounts, displs_send, coulomb_region_type );


		for ( int idx_state_vb = 0; idx_state_vb < 2*N_states[0]; idx_state_vb++ )
		{
		
			free( eigenstate_vb_padded[idx_state_vb].fourier );
			free( eigenstate_vb_padded_conj[idx_state_vb].fourier );

		}
		for ( int idx_state_cb = 0; idx_state_cb < 2*N_states[1]; idx_state_cb++ )
		{
		
			free( eigenstate_cb_padded[idx_state_cb].fourier );
			free( eigenstate_cb_padded_conj[idx_state_cb].fourier );

		}
		free( eigenstate_vb_padded );
		free( eigenstate_vb_padded_conj );
		free( eigenstate_cb_padded );
		free( eigenstate_cb_padded_conj );

	}
	else // Processes with IDs greater than or equal to the process limit
	{

		// Function call dummy variables as processes will not take part in Dfti calculations
		MKL_LONG local_memory_size = 0;
		MKL_LONG local_start = 0;
		MKL_LONG local_rows = 0;
		MKL_LONG local_size = 0;

		// Redistributing the block-cyclically distributed envelope functions that were output from the eigensolver so that they are distributed among processes over the first spatial dimension, enabling use of the Dfti library
		state_redistribution( basis, supercell, N_states, eigenstate_vb, eigenstate_cb, NULL, NULL, NULL, NULL, process_limit, NULL, local_memory_size, local_start, local_rows );

		for ( int idx_state = 0; idx_state < N_state_conv; idx_state++ )
		{

			convolution_states[idx_state].conv_state =  calloc( local_size_conv  , sizeof( complex double ) );
			convolution_states_reversed[idx_state].conv_state =  calloc( local_size_conv  , sizeof( complex double ) );

		}

		// Dummy Dfti task
		DFTI_DESCRIPTOR_DM_HANDLE task = NULL;

		// Calculating B arrays
		B_convolutions( basis, coulomb_basis, supercell, N_states, NULL, NULL, NULL, NULL, process_num, process_limit, task, local_memory_size, local_size, convolution_states, convolution_states_reversed, sendcounts, displs_send, coulomb_region_type );
		

	}

	double finish = MPI_Wtime();
    double elapsed = finish - start;
	component_run_times.cme_B = elapsed;



	///////////////////////////
	//// Vijkl calculation ////
	///////////////////////////

	start = MPI_Wtime();

	// Will store components of multipole expansion corrections
	struct coulomb_corrections CME_corrections;

	if ( include_CME >= 2 ) // If monopole correction is included
	{

		// Ewald parameter
		double ewald = sqrt( M_PI / ( supercell[0].CL1*supercell[0].CL2 ) );
		if ( myid == 0 ) // Only one process will contribute the following part of the Madelung term
		{

			CME_corrections.a_mad = -2*ewald / sqrt( M_PI ) - M_PI / ( ewald*ewald*supercell[0].CL1*supercell[0].CL2*supercell[0].CL3 );

		}
		else
		{

			CME_corrections.a_mad = 0;

		}

		// Stores real space sum
		double a_sum = 0;

		// Calculating contribution of each process to real space sum
		#pragma omp parallel for reduction(+:a_sum)
		for ( int idx_wav = 0; idx_wav < local_size_conv; idx_wav++ )
		{

			int idx_Gx = ( displs_send[myid] + idx_wav ) / ( ( 2*coulomb_basis.N_Gy - 1 )*( 2*coulomb_basis.N_Gz - 1 ) );
			int idx_Gy = ( ( displs_send[myid] + idx_wav ) % ( ( 2*coulomb_basis.N_Gy - 1 )*( 2*coulomb_basis.N_Gz - 1 )) ) / ( 2*coulomb_basis.N_Gz - 1 );
			int idx_Gz = ( ( displs_send[myid] + idx_wav ) % ( ( 2*coulomb_basis.N_Gy - 1 )*( 2*coulomb_basis.N_Gz - 1 )) ) % ( 2*coulomb_basis.N_Gz - 1 );

			int idx_basis_Gx = idx_Gx;
			int idx_basis_Gy = idx_Gy;
			int idx_basis_Gz = idx_Gz;

			double Gx = coulomb_basis.Gx[idx_basis_Gx];
			double Gy = coulomb_basis.Gy[idx_basis_Gy];
			double Gz = coulomb_basis.Gz[idx_basis_Gz];


			double G = sqrt( Gx*Gx + Gy*Gy + Gz*Gz );

			if ( idx_basis_Gx != coulomb_basis.N_Gx_T / 2 || idx_basis_Gy != coulomb_basis.N_Gy_T / 2 || idx_basis_Gz != coulomb_basis.N_Gz_T / 2 )
			{

				a_sum += ( 4*M_PI / ( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3 ) ) * exp( -G*G / ( 4*ewald*ewald ) ) / ( G*G );

			}

		}
		

		// Variables to determine process local wave vecttor indices for reciprocal space component of the Madelung term
		int local_size_m = 0;
		int displs_m = 0;

		for ( int proc = 0; proc < myid; proc++ )
		{

			displs_m += local_size_m;
			local_size_m = N_GC_conv / nprocs + ( ( proc < N_GC_conv % nprocs ) ? 1 : 0 );
		
		}

		// Calculating contribution of each process to reciprocal space sum
		#pragma omp parallel for reduction(+:a_sum)
		for ( int idx_m = 0; idx_m < local_size_m; idx_m++ )
		{

			int idx_mx = ( displs_m + idx_m ) / ( ( coulomb_basis.N_Gy_T / 2 )*( coulomb_basis.N_Gz_T / 2 ) );
			int idx_my = ( ( displs_m + idx_m ) % ( ( coulomb_basis.N_Gy_T / 2 )*( coulomb_basis.N_Gz_T / 2 ) ) ) / ( coulomb_basis.N_Gz_T / 2 );
			int idx_mz = ( ( displs_m + idx_m ) % ( ( coulomb_basis.N_Gy_T / 2 )*( coulomb_basis.N_Gz_T / 2 ) ) ) % ( coulomb_basis.N_Gz_T / 2 );

			if ( idx_mx != 0 || idx_my != 0 || idx_mz != 0 )
			{
				double R = sqrt( pow( idx_mx*supercell[0].CL1, 2 ) + pow( idx_my*supercell[0].CL2, 2 ) + pow( idx_mz*supercell[0].CL3, 2 ) );

				a_sum += erfc( R*ewald ) / R;

			}

		}

		CME_corrections.a_mad += a_sum;

	}
	

	if ( include_CME >= 3 ) // If dipole correction is included
	{

		// Local x, y, and z components of I_il wave vector for dipole correction (see Vukmirovic and Tomic 2008 (doi:10.1063/1.2936318))
		CME_corrections.Ix = calloc( local_size_conv  , sizeof( complex double ) );
		CME_corrections.Iy = calloc( local_size_conv  , sizeof( complex double ) );
		CME_corrections.Iz = calloc( local_size_conv  , sizeof( complex double ) );

		if ( myid == 0 )
		{

			// x, y, and z components of d_il vector for dipole correction
			CME_corrections.dx = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );
			CME_corrections.dy = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );
			CME_corrections.dz = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );

		}

	}
	if ( include_CME >= 4 ) // If quadrupole correction is included
	{

		// Stores local sum of x, y, and z components of J vector for quadrupole correction
		CME_corrections.J = calloc( local_size_conv  , sizeof( complex double ) );
		// Stores local Q term for quadrupole correction
		CME_corrections.Q = calloc( N_state_conv  , sizeof( complex double ) );

	}

	// Stores local components of Vijkl factor with the squared wave vector on the denominator
	complex double *wavevector_prefactor =  calloc( local_size_conv  , sizeof( complex double ) );
	
	// Looping over local wave vectors to calculate local components of Vijkl factor with the squared wave vector on the denominator, and calculate array terms for multipole expansion corrections
	#pragma omp parallel for
	for ( int idx_wav = 0; idx_wav < local_size_conv; idx_wav++ )
	{

		// Wave vector component indices
		int idx_Gx = ( displs_send[myid] + idx_wav ) / ( ( 2*coulomb_basis.N_Gy - 1 )*( 2*coulomb_basis.N_Gz - 1 ) );
		int idx_Gy = ( ( displs_send[myid] + idx_wav ) % ( ( 2*coulomb_basis.N_Gy - 1 )*( 2*coulomb_basis.N_Gz - 1 )) ) / ( 2*coulomb_basis.N_Gz - 1 );
		int idx_Gz = ( ( displs_send[myid] + idx_wav ) % ( ( 2*coulomb_basis.N_Gy - 1 )*( 2*coulomb_basis.N_Gz - 1 )) ) % ( 2*coulomb_basis.N_Gz - 1 );

		int idx_basis_Gx = idx_Gx;
		int idx_basis_Gy = idx_Gy;
		int idx_basis_Gz = idx_Gz;

		// Wave vector components
		double Gx = coulomb_basis.Gx[idx_basis_Gx];
		double Gy = coulomb_basis.Gy[idx_basis_Gy];
		double Gz = coulomb_basis.Gz[idx_basis_Gz];

		// Magnitude of wave vector
		double G = sqrt( Gx*Gx + Gy*Gy + Gz*Gz );

		// Calculating Vijkl factor with the squared wave vector on the denominator
		if ( idx_basis_Gx != coulomb_basis.N_Gx_T / 2 || idx_basis_Gy != coulomb_basis.N_Gy_T / 2 || idx_basis_Gz != coulomb_basis.N_Gz_T / 2 )
		{

			wavevector_prefactor[idx_wav] = ESQEPS / ( ( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3 )*G*G*pow( parameters[0].nr, 2.0 ) );

		}
		else
		{

			wavevector_prefactor[idx_wav] = 0;

		}


		// If dipole correction is included, calculate local I_il components
		if ( include_CME >= 3 ) 
		{

			if ( idx_basis_Gx != coulomb_basis.N_Gx_T / 2 && idx_basis_Gy == coulomb_basis.N_Gy_T / 2 && idx_basis_Gz == coulomb_basis.N_Gz_T / 2 )
			{

				CME_corrections.Ix[idx_wav] = supercell[0].CL2*supercell[0].CL3*( ( supercell[0].CL1 / ( I*Gx ) )*( cos( 0.5*Gx*supercell[0].CL1 ) ) + ( 2.0*I / ( Gx*Gx ) )*( sin( 0.5*Gx*supercell[0].CL1 ) ) );

			}

			if ( idx_basis_Gx == coulomb_basis.N_Gx_T / 2 && idx_basis_Gy != coulomb_basis.N_Gy_T / 2 && idx_basis_Gz == coulomb_basis.N_Gz_T / 2 )
			{

				CME_corrections.Iy[idx_wav] = supercell[0].CL1*supercell[0].CL3*( ( supercell[0].CL2 / ( I*Gy ) )*( cos( 0.5*Gy*supercell[0].CL2 ) ) + ( 2.0*I / ( Gy*Gy ) )*( sin( 0.5*Gy*supercell[0].CL2 ) ) );

			}

			if ( idx_basis_Gx == coulomb_basis.N_Gx_T / 2 && idx_basis_Gy == coulomb_basis.N_Gy_T / 2 && idx_basis_Gz != coulomb_basis.N_Gz_T / 2 )
			{

				CME_corrections.Iz[idx_wav] = supercell[0].CL2*supercell[0].CL1*( ( supercell[0].CL3 / ( I*Gz ) )*( cos( 0.5*Gz*supercell[0].CL3 ) ) + ( 2.0*I / ( Gz*Gz ) )*( sin( 0.5*Gz*supercell[0].CL3 ) ) );

			}

		}


		// If quadrupole correction is included, calculate local J components
		if ( include_CME >= 4 )
		{

			if ( idx_basis_Gx != coulomb_basis.N_Gx_T / 2 && idx_basis_Gy == coulomb_basis.N_Gy_T / 2 && idx_basis_Gz == coulomb_basis.N_Gz_T / 2 )
			{

				CME_corrections.J[idx_wav] += supercell[0].CL2*supercell[0].CL3*( ( 0.5*supercell[0].CL1*supercell[0].CL1 / Gx )*sin( 0.5*supercell[0].CL1*Gx ) + ( 2*I / Gx )*( ( supercell[0].CL1 / ( I*Gx ) )*( cos( 0.5*Gx*supercell[0].CL1 ) ) + ( 2.0*I / ( Gx*Gx ) )*( sin( 0.5*Gx*supercell[0].CL1 ) ) ) );

			}
			else if ( idx_basis_Gx == coulomb_basis.N_Gx_T / 2 && idx_basis_Gy == coulomb_basis.N_Gy_T / 2 && idx_basis_Gz == coulomb_basis.N_Gz_T / 2 )
			{

				CME_corrections.J[idx_wav] += ( 1 / 12.0 )*pow( supercell[0].CL1, 3 )*supercell[0].CL2*supercell[0].CL3;

			}
			

			if ( idx_basis_Gx == coulomb_basis.N_Gx_T / 2 && idx_basis_Gy != coulomb_basis.N_Gy_T / 2 && idx_basis_Gz == coulomb_basis.N_Gz_T / 2 )
			{

				CME_corrections.J[idx_wav] += supercell[0].CL1*supercell[0].CL3*( ( 0.5*supercell[0].CL2*supercell[0].CL2 / Gy )*sin( 0.5*supercell[0].CL2*Gy ) + ( 2*I / Gy )*( ( supercell[0].CL2 / ( I*Gy ) )*( cos( 0.5*Gy*supercell[0].CL2 ) ) + ( 2.0*I / ( Gy*Gy ) )*( sin( 0.5*Gy*supercell[0].CL2 ) ) ) );

			}
			else if ( idx_basis_Gx == coulomb_basis.N_Gx_T / 2 && idx_basis_Gy == coulomb_basis.N_Gy_T / 2 && idx_basis_Gz == coulomb_basis.N_Gz_T / 2 )
			{

				CME_corrections.J[idx_wav] += ( 1 / 12.0 )*pow( supercell[0].CL2, 3 )*supercell[0].CL1*supercell[0].CL3;

			}

			if ( idx_basis_Gx == coulomb_basis.N_Gx_T / 2 && idx_basis_Gy == coulomb_basis.N_Gy_T / 2 && idx_basis_Gz != coulomb_basis.N_Gz_T / 2 )
			{

				CME_corrections.J[idx_wav] += supercell[0].CL2*supercell[0].CL1*( ( 0.5*supercell[0].CL3*supercell[0].CL3 / Gz )*sin( 0.5*supercell[0].CL3*Gz ) + ( 2*I / Gz )*( ( supercell[0].CL3 / ( I*Gz ) )*( cos( 0.5*Gz*supercell[0].CL3 ) ) + ( 2.0*I / ( Gz*Gz ) )*( sin( 0.5*Gz*supercell[0].CL3 ) ) ) );

			}
			else if ( idx_basis_Gx == coulomb_basis.N_Gx_T / 2 && idx_basis_Gy == coulomb_basis.N_Gy_T / 2 && idx_basis_Gz == coulomb_basis.N_Gz_T / 2 )
			{

				CME_corrections.J[idx_wav] += ( 1 / 12.0 )*pow( supercell[0].CL3, 3 )*supercell[0].CL2*supercell[0].CL1;

			}

		}

	}

	// Array to store local part of wave vector reversed Bjk function 
	complex double *convolution_state_reversed_wavevector =  calloc( local_size_conv  , sizeof( complex double ) );


	// Local components of d_il terms
	complex double *dx_local;
	complex double *dy_local;
	complex double *dz_local;
	if ( include_CME >= 3)
	{

		dx_local = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );
		dy_local = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );
		dz_local = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );

	}
	
	// Loop over B indices
	for ( int k = 0; k < 2*N_states[0] + 2*N_states[1]; k++ )
	{

		for ( int l = 0; l < 2*N_states[0] + 2*N_states[1]; l++ )
		{

			if ( l >= k )
			{

				// Element-by-element multiplication of wave-vector reversed Bjk terms and factor with the squared wave vector on the denominator
				vzMul( local_size_conv, convolution_states_reversed[(l * (l + 1) / 2) + k].conv_state, wavevector_prefactor, convolution_state_reversed_wavevector );
				// Updating wave vector reversed Bjk to include factor with the squared wave vector on the denominator, so that Vijk can be later simply calculated by a dot product
				cblas_zcopy( local_size_conv, convolution_state_reversed_wavevector, 1, convolution_states_reversed[(l * (l + 1) / 2) + k].conv_state, 1 );

			}

			// If dipole correction is included, calculate local d_il components
			if ( include_CME >= 3 )
			{

				if ( l >= k )
				{

					cblas_zdotu_sub( local_size_conv, convolution_states[(l * (l + 1) / 2) + k].conv_state, 1, CME_corrections.Ix, 1, &dx_local[l * (2*N_states[0] + 2*N_states[1]) + k] );
					cblas_zdotu_sub( local_size_conv, convolution_states[(l * (l + 1) / 2) + k].conv_state, 1, CME_corrections.Iy, 1, &dy_local[l * (2*N_states[0] + 2*N_states[1]) + k] );
					cblas_zdotu_sub( local_size_conv, convolution_states[(l * (l + 1) / 2) + k].conv_state, 1, CME_corrections.Iz, 1, &dz_local[l * (2*N_states[0] + 2*N_states[1]) + k] );

				}
				else
				{

					cblas_zdotc_sub( local_size_conv, convolution_states[(k * (k + 1) / 2) + l].conv_state, 1, CME_corrections.Ix, 1, &dx_local[l * (2*N_states[0] + 2*N_states[1]) + k] );
					cblas_zdotc_sub( local_size_conv, convolution_states[(k * (k + 1) / 2) + l].conv_state, 1, CME_corrections.Iy, 1, &dy_local[l * (2*N_states[0] + 2*N_states[1]) + k] );
					cblas_zdotc_sub( local_size_conv, convolution_states[(k * (k + 1) / 2) + l].conv_state, 1, CME_corrections.Iz, 1, &dz_local[l * (2*N_states[0] + 2*N_states[1]) + k] );

				}

			}

			// If quadrupole correction is included, calculate local Q_il components
			if ( include_CME >= 4 && l >= k )
			{

				cblas_zdotu_sub( local_size_conv, convolution_states[(l * (l + 1) / 2) + k].conv_state, 1, CME_corrections.J, 1, &CME_corrections.Q[(l * (l + 1) / 2) + k] );

			}

		}

	}

	// If dipole correction is included, calculate global d_il components from local
	if ( include_CME >= 3 )
	{

		MPI_Reduce( dx_local, CME_corrections.dx, ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] ), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( dy_local, CME_corrections.dy, ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] ), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );
		MPI_Reduce( dz_local, CME_corrections.dz, ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] ), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );

		free( dx_local );
		free( dy_local );
		free( dz_local );

	}


	// Calculate Vijkl

	// HHHH
	coulomb_integrals( supercell, parameters, N_states, convolution_states, convolution_states_reversed, local_size_conv, CME_corrections, CME, 0 );

	// EEEE
	coulomb_integrals( supercell, parameters, N_states, convolution_states, convolution_states_reversed, local_size_conv, CME_corrections, CME, 1 );

	// HEEH
	coulomb_integrals( supercell, parameters, N_states, convolution_states, convolution_states_reversed, local_size_conv, CME_corrections, CME, 2 );

	// HEHE
	coulomb_integrals( supercell, parameters, N_states, convolution_states, convolution_states_reversed, local_size_conv, CME_corrections, CME, 3 );
	

	finish = MPI_Wtime();
    elapsed = finish - start;
	component_run_times.cme_I = elapsed;


	for ( int idx_state = 0; idx_state < N_state_conv; idx_state++ )
	{

		free( convolution_states[idx_state].conv_state );
		free( convolution_states_reversed[idx_state].conv_state );

	}

	free( convolution_states );
	free( convolution_states_reversed);
	free( wavevector_prefactor );
	free( convolution_state_reversed_wavevector );

	if ( include_CME >= 3 )
	{

		free( CME_corrections.Ix );
		free( CME_corrections.Iy );
		free( CME_corrections.Iz );
		
		if ( myid == 0 )
		{

			free( CME_corrections.dx );
			free( CME_corrections.dy );
			free( CME_corrections.dz );

		}

	}
	if ( include_CME >= 4 )
	{

		free( CME_corrections.J );
		free( CME_corrections.Q );

	}
	
	free( sendcounts );
	free( displs_send );

}
