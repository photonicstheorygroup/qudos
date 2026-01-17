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

void piezo_region_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	MKL_INT process_limit = ( 2*basis.N_Gx_T - 1 <= 2*basis.N_Gy_T - 1 ? 2*basis.N_Gx_T - 1 : 2*basis.N_Gy_T - 1 );
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

	MKL_INT local_start_reduced;
	MKL_INT local_rows_reduced;
	MKL_INT local_size_reduced;


	if ( myid < process_limit )
	{

		MKL_LONG status;

		MKL_LONG sizes[] = { 2*basis.N_Gx_T - 1, 2*basis.N_Gy_T - 1, 2*basis.N_Gz_T - 1 };

		MKL_LONG global_size = sizes[0]*sizes[1]*sizes[2];

		MKL_LONG DFTI_PRECISION = DFTI_DOUBLE;
		MKL_LONG DFTI_TYPE = DFTI_COMPLEX;

		DFTI_DESCRIPTOR_DM_HANDLE task;

		status = DftiCreateDescriptorDM( MPI_SUB_COMM, &task, DFTI_PRECISION, DFTI_TYPE, 3, sizes );

		status = DftiSetValueDM( task, DFTI_BACKWARD_SCALE, 1.0 / global_size );

		MKL_LONG local_memory_size;
		status = DftiGetValueDM( task, CDFT_LOCAL_SIZE, &local_memory_size );
		MKL_LONG local_start;
		status = DftiGetValueDM( task, CDFT_LOCAL_X_START, &local_start );
		MKL_LONG local_rows;
		status = DftiGetValueDM( task, CDFT_LOCAL_NX, &local_rows );
		MKL_LONG local_size = local_rows*sizes[1]*sizes[2];

		local_mem += ( 7 )*( local_memory_size ) *sizeof( complex double ); // Adding memory to store local parts of strain tensor components and shape characteristic function


		if ( local_start >= basis.N_Gx_T / 2 && local_start < basis.N_Gx_T + basis.N_Gx_T / 2 )
		{

			local_start_reduced = 0;
			if ( local_start + local_rows < basis.N_Gx_T + basis.N_Gx_T / 2 )
			{

				local_rows_reduced = local_rows;

			}
			else
			{

				local_rows_reduced = basis.N_Gx_T + basis.N_Gx_T / 2 - local_start;

			}

		}
		else if ( local_start < basis.N_Gx_T / 2)
		{

			if ( local_start + local_rows >= basis.N_Gx_T / 2 && local_start + local_rows < basis.N_Gx_T + basis.N_Gx_T / 2 )
			{

				local_start_reduced = basis.N_Gx_T / 2 - local_start;
				local_rows_reduced = local_rows - local_start_reduced;

			}
			else if ( local_start + local_rows >= basis.N_Gx_T + basis.N_Gx_T / 2 )
			{

				local_start_reduced = basis.N_Gx_T / 2 - local_start;
				local_rows_reduced = basis.N_Gx_T;

			}
			else
			{

				local_start_reduced = 0;
				local_rows_reduced = 0;

			}

		}
		else
		{

			local_start_reduced = 0;
			local_rows_reduced = 0;

		}
		local_size_reduced = local_rows_reduced*basis.N_Gy_T*basis.N_Gz_T;


		local_mem += ( 6 )*( local_memory_size ) *sizeof( complex double ); 	// Adding memory to store local parts of convolutions between strain tensor components and shape characteristic function
		local_mem += ( 7 )*( local_size_reduced ) *sizeof( complex double );	// Adding memory to store local parts of central portion of convolutions between strain tensor components and shape characteristic function, and the local part of the piezoelectric potential


		local_mem_peak = local_mem;


		if ( include_piezo == 2 || include_piezo == 3 )
		{

			local_mem += ( 12 )*( local_memory_size ) *sizeof( complex double );	// Adding memory to store local parts of convolutions between strain tensor components
			local_mem += ( 12 )*( local_size_reduced ) *sizeof( complex double );	// Adding memory to store local parts of central portion of convolutions between strain tensor components

			local_mem_peak = local_mem;

			local_mem -= ( 12 )*( local_memory_size ) *sizeof( complex double );

			local_mem += ( 12 )*( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T ) *sizeof( complex double ); // Adding memory to store global convolutions between strain tensor components
			if ( local_mem > local_mem_peak )
			{

				local_mem_peak = local_mem; // Storing peak memory

			}

			local_mem -= ( 12 )*( local_size_reduced ) *sizeof( complex double );


			local_mem += ( 12 )*( local_memory_size ) *sizeof( complex double );	// Adding memory to store local parts of convolutions between strain tensor components for future calculation of convolutions with shape characteristic functions
			if ( local_mem > local_mem_peak )
			{

				local_mem_peak = local_mem; // Storing peak memory

			}

			local_mem -= ( 12 )*( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T ) *sizeof( complex double );


			local_mem += ( 12 )*( local_memory_size ) *sizeof( complex double );	// Adding memory to store local parts of convolutions of convolutions between strain tensor components with shape characteristic functions


			local_mem += ( 12 )*( local_size_reduced ) *sizeof( complex double );	// Adding memory to store local parts of central portion of convolutions of convolutions between strain tensor components with shape characteristic functions

			if ( local_mem > local_mem_peak )
			{

				local_mem_peak = local_mem; // Storing peak memory

			}


		}


		local_mem -= ( 7 )*( local_memory_size ) *sizeof( complex double );

		local_mem -= ( 6 )*( local_memory_size ) *sizeof( complex double );

		local_mem -= ( 7 )*( local_size_reduced ) *sizeof( complex double );

		status = DftiFreeDescriptorDM( &task );

		
		if ( include_piezo == 2 || include_piezo == 3 )
		{

			local_mem -= ( 24 )*( local_memory_size ) *sizeof( complex double );

			local_mem -= ( 12 )*( local_size_reduced ) *sizeof( complex double );
			
		}

	}

	MPI_Comm_free( &MPI_SUB_COMM );

}
