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

	int N_G = basis.N_Gx*basis.N_Gy*basis.N_Gz;
	int N_GC_conv = ( 2*coulomb_basis.N_Gx - 1 )*( 2*coulomb_basis.N_Gy - 1 )*( 2*coulomb_basis.N_Gz - 1 );

	MKL_INT process_limit = ( 2*basis.N_Gx - 1 <= 2*basis.N_Gy - 1 ? 2*basis.N_Gx - 1 : 2*basis.N_Gy - 1 );
	MKL_INT process_num;

	MPI_Comm MPI_SUB_COMM;
	if ( nprocs > process_limit )
	{

		MPI_Comm_split( MPI_COMM_WORLD, ( myid < process_limit ? 0 : 1 ), myid, &MPI_SUB_COMM );
		process_num = process_limit;

	}
	else
	{

		MPI_Comm_dup( MPI_COMM_WORLD, &MPI_SUB_COMM );
		process_num = nprocs;

	}

	int *sendcounts;// = ( int * ) malloc( nprocs*sizeof( int ) );
	int *displs_send;// = ( int * ) malloc( nprocs*sizeof( int ) );


	MKL_INT N_state_conv = ( ( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] + 1 ) ) / 2 );

	struct convolution_state *convolution_states;// = malloc( N_state_conv*sizeof( _convolution_state ) );
	struct convolution_state *convolution_states_reversed;// = malloc( N_state_conv*sizeof( _convolution_state ) );

	int local_size_conv = N_GC_conv / nprocs + ( ( myid < N_GC_conv % nprocs ) ? 1 : 0 );


	if ( myid < process_limit )
	{

		MKL_LONG status;

		MKL_LONG sizes[] = { 2*basis.N_Gx - 1, 2*basis.N_Gy - 1, 2*basis.N_Gz - 1 };

		MKL_LONG global_size = sizes[0]*sizes[1]*sizes[2];

		MKL_LONG DFTI_PRECISION = DFTI_DOUBLE;
		MKL_LONG DFTI_TYPE = DFTI_COMPLEX;

		DFTI_DESCRIPTOR_DM_HANDLE task;

		status = DftiCreateDescriptorDM( MPI_SUB_COMM, &task, DFTI_PRECISION, DFTI_TYPE, 3, sizes );

		//status = DftiSetValueDM( task, DFTI_BACKWARD_SCALE, 1.0 / global_size );

		MKL_LONG local_memory_size;
		status = DftiGetValueDM( task, CDFT_LOCAL_SIZE, &local_memory_size );
		MKL_LONG local_start;
		status = DftiGetValueDM( task, CDFT_LOCAL_X_START, &local_start );
		MKL_LONG local_rows;
		status = DftiGetValueDM( task, CDFT_LOCAL_NX, &local_rows );
		MKL_LONG local_size = local_rows*sizes[1]*sizes[2];

		//status = DftiCommitDescriptorDM( task );


		struct eigenstate *eigenstate_vb_padded;// = malloc( 2*N_states[0]*sizeof( _eigenstate ) );
		struct eigenstate *eigenstate_cb_padded;// = malloc( 2*N_states[1]*sizeof( _eigenstate ) );
		struct eigenstate *eigenstate_vb_padded_conj;// = malloc( 2*N_states[0]*sizeof( _eigenstate ) );
		struct eigenstate *eigenstate_cb_padded_conj;// = malloc( 2*N_states[1]*sizeof( _eigenstate ) );

		for ( int idx_state_vb = 0; idx_state_vb < 2*N_states[0]; idx_state_vb++ )
		{
		
			//eigenstate_vb_padded[idx_state_vb].fourier   = calloc( N_BAND*local_memory_size, sizeof( complex double ) );
			//eigenstate_vb_padded_conj[idx_state_vb].fourier   = calloc( N_BAND*local_memory_size, sizeof( complex double ) );

		}
		for ( int idx_state_cb = 0; idx_state_cb < 2*N_states[1]; idx_state_cb++ )
		{
		
			//eigenstate_cb_padded[idx_state_cb].fourier   = calloc( N_BAND*local_memory_size, sizeof( complex double ) );
			//eigenstate_cb_padded_conj[idx_state_cb].fourier   = calloc( N_BAND*local_memory_size, sizeof( complex double ) );

		}
		local_mem += ( 4*N_states[0]*N_BAND*local_memory_size )*sizeof( complex double );	// Adding memory to store local parts of Fourier coefficients of VB state envelope functions
		local_mem += ( 4*N_states[1]*N_BAND*local_memory_size )*sizeof( complex double );	// Adding memory to store local parts of Fourier coefficients of CB state envelope functions

		//state_redistribution( basis, supercell, N_states, eigenstate_vb, eigenstate_cb, eigenstate_vb_padded, eigenstate_cb_padded, eigenstate_vb_padded_conj, eigenstate_cb_padded_conj, process_limit, sizes, local_memory_size, local_start, local_rows );

		local_mem += ( N_H )*sizeof( complex double );	// Adding memory to store global interim eigenstate for state redistribution

		if ( local_mem > local_mem_peak )
		{
			
			local_mem_peak = local_mem; // Store if peak memory

		}

		local_mem -= ( N_H )*sizeof( complex double );


		for ( int idx_state = 0; idx_state < N_state_conv; idx_state++ )
		{

			//convolution_states[idx_state].conv_state =  calloc( local_size_conv  , sizeof( complex double ) );
			//convolution_states_reversed[idx_state].conv_state =  calloc( local_size_conv  , sizeof( complex double ) );

		}

		local_mem += ( 2*N_state_conv*local_size_conv )*sizeof( complex double ); // Adding memory for local parts of B arrays


		B_convolutions( basis, coulomb_basis, supercell, N_states, eigenstate_vb_padded, eigenstate_cb_padded, eigenstate_vb_padded_conj, eigenstate_cb_padded_conj, process_num, process_limit, task, local_memory_size, local_size, convolution_states, convolution_states_reversed, sendcounts, displs_send, coulomb_region_type );


		for ( int idx_state_vb = 0; idx_state_vb < 2*N_states[0]; idx_state_vb++ )
		{
		
			//free( eigenstate_vb_padded[idx_state_vb].fourier );
			//free( eigenstate_vb_padded_conj[idx_state_vb].fourier );
			

		}
		for ( int idx_state_cb = 0; idx_state_cb < 2*N_states[1]; idx_state_cb++ )
		{
		
			//free( eigenstate_cb_padded[idx_state_cb].fourier );
			//free( eigenstate_cb_padded_conj[idx_state_cb].fourier );

		}
		//free( eigenstate_vb_padded );
		//free( eigenstate_vb_padded_conj );
		//free( eigenstate_cb_padded );
		//free( eigenstate_cb_padded_conj );

		local_mem -= ( 4*N_states[0]*N_BAND*local_memory_size )*sizeof( complex double );
		local_mem -= ( 4*N_states[1]*N_BAND*local_memory_size )*sizeof( complex double );

	}
	else
	{

		MKL_LONG local_memory_size = 0;
		MKL_LONG local_start = 0;
		MKL_LONG local_rows = 0;
		MKL_LONG local_size = 0;

		//state_redistribution( basis, supercell, N_states, eigenstate_vb, eigenstate_cb, NULL, NULL, NULL, NULL, process_limit, NULL, local_memory_size, local_start, local_rows );

		local_mem += ( N_H )*sizeof( complex double ); // Adds memory to store interim eigenstate for eigenstate redistribution among MPI processes

		if ( local_mem > local_mem_peak )
		{
			
			local_mem_peak = local_mem; // Store if peak memory

		}

		local_mem -= ( N_H )*sizeof( complex double );

		for ( int idx_state = 0; idx_state < N_state_conv; idx_state++ )
		{

			//convolution_states[idx_state].conv_state =  calloc( local_size_conv  , sizeof( complex double ) );
			//convolution_states_reversed[idx_state].conv_state =  calloc( local_size_conv  , sizeof( complex double ) );

		}
		local_mem += ( 2*N_state_conv*local_size_conv )*sizeof( complex double );

		DFTI_DESCRIPTOR_DM_HANDLE task;

		B_convolutions( basis, coulomb_basis, supercell, N_states, NULL, NULL, NULL, NULL, process_num, process_limit, task, local_memory_size, local_size, convolution_states, convolution_states_reversed, sendcounts, displs_send, coulomb_region_type );
		

	}



	struct coulomb_corrections CME_corrections;
	

	if ( include_CME >= 3 )
	{

		//CME_corrections.Ix = calloc( local_size_conv  , sizeof( complex double ) );
		//CME_corrections.Iy = calloc( local_size_conv  , sizeof( complex double ) );
		//CME_corrections.Iz = calloc( local_size_conv  , sizeof( complex double ) );

		local_mem += ( 3*local_size_conv )*sizeof( complex double ); // Add memory to store local x, y, and z components of I_il wave vector for dipole correction

		if ( myid == 0 )
		{

			//CME_corrections.dx = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );
			//CME_corrections.dy = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );
			//CME_corrections.dz = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );

			local_mem += 3*(  2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )*sizeof( complex double ); // Add memory to store x, y, and z components of d_il vector for dipole correction

		}

	}
	if ( include_CME >= 4 )
	{

		//CME_corrections.J = calloc( local_size_conv  , sizeof( complex double ) );
		//CME_corrections.Q = calloc( N_state_conv  , sizeof( complex double ) );

		local_mem += ( local_size_conv )*sizeof( complex double ); // Add memory to store local sum of x, y, and z components of J vector for quadrupole correction
		local_mem += ( N_state_conv )*sizeof( complex double ); // Add memory to store local Q term for quadrupole correction

	}

	complex double *wavevector_prefactor;// =  calloc( local_size_conv  , sizeof( complex double ) );
	local_mem += ( local_size_conv )*sizeof( complex double ); // Add memory to store local components of Vijkl factor with the squared wave vector on the denominator


	complex double *convolution_state_reversed_wavevector;// =  calloc( local_size_conv  , sizeof( complex double ) );
	local_mem += ( local_size_conv )*sizeof( complex double ); // Add memory to store to store local part of wave vector reversed Bjk function 


	complex double *dx_local;
	complex double *dy_local;
	complex double *dz_local;
	if ( include_CME >= 3)
	{

		//dx_local = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );
		//dy_local = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );
		//dz_local = calloc( ( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )  , sizeof( complex double ) );

		local_mem += 3*( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )*sizeof( complex double ); // Add memory to store local components of d_il terms

	}

	if ( local_mem > local_mem_peak )
	{
		
		local_mem_peak = local_mem; // Store if peak memory

	}
	

	if ( include_CME >= 3 )
	{

		//free( dx_local );
		//free( dy_local );
		//free( dz_local );

		local_mem -= 3*( 2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )*sizeof( complex double );

	}


	// HHHH
	//coulomb_integrals( supercell, parameters, N_states, convolution_states, convolution_states_reversed, local_size_conv, CME_corrections, CME, 0 );

	// EEEE
	//coulomb_integrals( supercell, parameters, N_states, convolution_states, convolution_states_reversed, local_size_conv, CME_corrections, CME, 1 );

	// HHEE
	//coulomb_integrals( supercell, parameters, N_states, convolution_states, convolution_states_reversed, local_size_conv, CME_corrections, CME, 2 );

	// HEEH
	//coulomb_integrals( supercell, parameters, N_states, convolution_states, convolution_states_reversed, local_size_conv, CME_corrections, CME, 3 );


	int N_CME_ee = 2*N_states[1]*2*N_states[1]*2*N_states[1];
	int N_CME_hh = 2*N_states[0]*2*N_states[0]*2*N_states[0];
	int N_CME_eh = 2*N_states[1]*2*N_states[0]*2*N_states[1];



	local_mem += ( N_CME_ee )*sizeof( complex double ); // Add memory to store Vijkl involving electron states only
	if ( local_mem > local_mem_peak )
	{
		
		local_mem_peak = local_mem; // Store if peak memory

	}
	local_mem -= ( N_CME_ee )*sizeof( complex double );

	local_mem += ( N_CME_hh )*sizeof( complex double ); // Add memory to store Vijkl involving hole states only
	if ( local_mem > local_mem_peak )
	{
		
		local_mem_peak = local_mem; // Store if peak memory

	}
	local_mem -= ( N_CME_hh )*sizeof( complex double );

	local_mem += ( N_CME_eh )*sizeof( complex double ); // Add memory to store Vijkl involving electron and hole states
	if ( local_mem > local_mem_peak )
	{
		
		local_mem_peak = local_mem; // Store if peak memory

	}
	local_mem -= ( N_CME_eh )*sizeof( complex double );




	for ( int idx_state = 0; idx_state < N_state_conv; idx_state++ )
	{

		//free( convolution_states[idx_state].conv_state );
		//free( convolution_states_reversed[idx_state].conv_state );

	}

	local_mem -= ( 2*N_state_conv*local_size_conv )*sizeof( complex double );

	//free( convolution_states );
	//free( convolution_states_reversed);
	//free( wavevector_prefactor );

	local_mem -= ( local_size_conv )*sizeof( complex double );

	//free( convolution_state_reversed_wavevector );

	local_mem -= ( local_size_conv )*sizeof( complex double );

	if ( include_CME >= 3 )
	{

		//free( CME_corrections.Ix );
		//free( CME_corrections.Iy );
		//free( CME_corrections.Iz );

		local_mem -= ( 3*local_size_conv )*sizeof( complex double );
		
		if ( myid == 0 )
		{

			//free( CME_corrections.dx );
			//free( CME_corrections.dy );
			//free( CME_corrections.dz );

			local_mem -= 3*(  2*N_states[0] + 2*N_states[1] )*( 2*N_states[0] + 2*N_states[1] )*sizeof( complex double );

		}

	}
	if ( include_CME >= 4 )
	{

		//free( CME_corrections.J );
		//free( CME_corrections.Q );

		local_mem -= ( local_size_conv )*sizeof( complex double );
		local_mem -= ( N_state_conv )*sizeof( complex double );

	}
	
	//free( sendcounts );
	//free( displs_send );

}
