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

void characteristic_functions_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, int strain_region )
{
	
	// Local size of arrays
	int local_size = basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs + ( ( myid < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? 1 : 0 );
	
	// Temporary array to store characteristic functions calculated by the process
	complex double *char_func;// = ( complex double * ) malloc( local_size*sizeof( complex double ) );

	local_mem += local_size*sizeof( complex double );  // Adding memory to store local portion of shape characteristic function
	local_mem_peak = local_mem;							// Peak process memory usage for shape characteristic functions calculation
	
	local_mem -= local_size*sizeof( complex double );

}
