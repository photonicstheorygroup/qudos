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

void piezo_region_reciprocal_space( int piezo_order, int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct characteristic_functions *characteristic_function_G, struct strain_tensor strain_G, complex double *piezo_G )
{

	struct convolution_char_strain convolutions_char_strain_region;

	// Arrays for convolutions between characteristic functions and strain tensor components
	convolutions_char_strain_region.char_xx = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
	convolutions_char_strain_region.char_yy = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
	convolutions_char_strain_region.char_zz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
	convolutions_char_strain_region.char_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
	convolutions_char_strain_region.char_xz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
	convolutions_char_strain_region.char_xy = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
	
	if ( piezo_order == 2 || piezo_order == 3 ) // If 2nd order piezoelectric potential included in calculation
	{
	
		// Arrays for convolutions between different strain tensor components
		convolutions_char_strain_region.xx_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.yy_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.zz_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.xy_xz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.yy_xz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.xx_xz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.zz_xz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.xy_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.zz_xy = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.xx_xy = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.yy_xy = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.xz_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );

		// Arrays for convolutions between characteristic functions and convolutions between different strain tensor components
		convolutions_char_strain_region.char_xx_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_yy_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_zz_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_xy_xz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_yy_xz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_xx_xz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_zz_xz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_xy_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_zz_xy = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_xx_xy = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_yy_xy = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		convolutions_char_strain_region.char_xz_yz = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );
		
	}

	
	VSLConvTaskPtr task; // vsl task pointer

	MKL_INT VSL_CONV_MODE = VSL_CONV_MODE_FFT; // Using FFT mode of vsl to compute linear convolution

	// Sizes of dimensions of two input and one output array (output array sized to only extract central portion of output)
	const MKL_INT xshape[] = {basis.N_Gx, basis.N_Gy, basis.N_Gz};
	const MKL_INT yshape[] = {basis.N_Gx, basis.N_Gy, basis.N_Gz};
	const MKL_INT zshape[] = {basis.N_Gx, basis.N_Gy, basis.N_Gz};

	const MKL_INT start[] = { basis.N_Gx/2, basis.N_Gy/2, basis.N_Gz/2 }; // Where to start taking data from output array (set to only extract central portion)

	// Strides of stored data (standard data layout)
	const MKL_INT xstride[] = {( basis.N_Gy )*( basis.N_Gz ), basis.N_Gz, 1};
	const MKL_INT ystride[] = {( basis.N_Gy )*( basis.N_Gz ), basis.N_Gz, 1};
	const MKL_INT zstride[] = {( basis.N_Gy )*( basis.N_Gz ), basis.N_Gz, 1};


	if ( piezo_order == 2 || piezo_order == 3 ) // If 2nd order piezoelectric potential included in calculation
	{

		// Create new vsl task (takes array as argument that will be involved in multiple convolution calculations)
		check_vsl( vslzConvNewTaskX( &task, VSL_CONV_MODE, 3, xshape, yshape, zshape, strain_G.xx, xstride ), "vslzConvNewTaskX" );
		// Set first value of stored output to be at "start"
		check_vsl( vslConvSetStart( task, start ), "vslConvSetStart" );

		// Compute convolutions between xx and yz,xz,xy strain tensor components
		check_vsl( vslzConvExecX( task, strain_G.yz, ystride, convolutions_char_strain_region.xx_yz, zstride ), "vslzConvExecX" );
		check_vsl( vslzConvExecX( task, strain_G.xz, ystride, convolutions_char_strain_region.xx_xz, zstride ), "vslzConvExecX" );
		check_vsl( vslzConvExecX( task, strain_G.xy, ystride, convolutions_char_strain_region.xx_xy, zstride ), "vslzConvExecX" );

		// Delete vsl task
		check_vsl( vslConvDeleteTask( &task ), "vslConvDeleteTask" );


		// Create new vsl task (takes array as argument that will be involved in multiple convolution calculations)
		check_vsl( vslzConvNewTaskX( &task, VSL_CONV_MODE, 3, xshape, yshape, zshape, strain_G.yy, xstride ), "vslzConvNewTaskX" );
		// Set first value of stored output to be at "start"
		check_vsl( vslConvSetStart( task, start ), "vslConvSetStart" );

		// Compute convolutions between yy and yz,xz,xy strain tensor components
		check_vsl( vslzConvExecX( task, strain_G.yz, ystride, convolutions_char_strain_region.yy_yz, zstride ), "vslzConvExecX" );
		check_vsl( vslzConvExecX( task, strain_G.xz, ystride, convolutions_char_strain_region.yy_xz, zstride ), "vslzConvExecX" );
		check_vsl( vslzConvExecX( task, strain_G.xy, ystride, convolutions_char_strain_region.yy_xy, zstride ), "vslzConvExecX" );

		// Delete vsl task
		check_vsl( vslConvDeleteTask( &task ), "vslConvDeleteTask" );


		// Create new vsl task (takes array as argument that will be involved in multiple convolution calculations)
		check_vsl( vslzConvNewTaskX( &task, VSL_CONV_MODE, 3, xshape, yshape, zshape, strain_G.zz, xstride ), "vslzConvNewTaskX" );
		// Set first value of stored output to be at "start"
		check_vsl( vslConvSetStart( task, start ), "vslConvSetStart" );

		// Compute convolutions between zz and yz,xz,xy strain tensor components
		check_vsl( vslzConvExecX( task, strain_G.yz, ystride, convolutions_char_strain_region.zz_yz, zstride ), "vslzConvExecX" );
		check_vsl( vslzConvExecX( task, strain_G.xz, ystride, convolutions_char_strain_region.zz_xz, zstride ), "vslzConvExecX" );
		check_vsl( vslzConvExecX( task, strain_G.xy, ystride, convolutions_char_strain_region.zz_xy, zstride ), "vslzConvExecX" );

		// Delete vsl task
		check_vsl( vslConvDeleteTask( &task ), "vslConvDeleteTask" );


		// Create new vsl task (takes array as argument that will be involved in multiple convolution calculations)
		check_vsl( vslzConvNewTaskX( &task, VSL_CONV_MODE, 3, xshape, yshape, zshape, strain_G.xy, xstride ), "vslzConvNewTaskX" );
		// Set first value of stored output to be at "start"
		check_vsl( vslConvSetStart( task, start ), "vslConvSetStart" );

		// Compute convolutions between xy and yz,xz strain tensor components
		check_vsl( vslzConvExecX( task, strain_G.yz, ystride, convolutions_char_strain_region.xy_yz, zstride ), "vslzConvExecX" );
		check_vsl( vslzConvExecX( task, strain_G.xz, ystride, convolutions_char_strain_region.xy_xz, zstride ), "vslzConvExecX" );

		// Delete vsl task
		check_vsl( vslConvDeleteTask( &task ), "vslConvDeleteTask" );


		// Create new vsl task (takes array as argument that will be involved in multiple convolution calculations)
		check_vsl( vslzConvNewTaskX( &task, VSL_CONV_MODE, 3, xshape, yshape, zshape, strain_G.xz, xstride ), "vslzConvNewTaskX" );
		// Set first value of stored output to be at "start"
		check_vsl( vslConvSetStart( task, start ), "vslConvSetStart" );

		// Compute convolutions between yz and xz strain tensor components
		check_vsl( vslzConvExecX( task, strain_G.yz, ystride, convolutions_char_strain_region.xz_yz, zstride ), "vslzConvExecX" );

		// Delete vsl task
		check_vsl( vslConvDeleteTask( &task ), "vslConvDeleteTask" );

	}


	// Loop over shapes
	for ( int idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		if ( idx_shape > 0 ) // If not supercell shape
		{

			// Create new vsl task (takes array as argument that will be involved in multiple convolution calculations)
			check_vsl( vslzConvNewTaskX( &task, VSL_CONV_MODE, 3, xshape, yshape, zshape, characteristic_function_G[idx_shape].chi, xstride ), "vslzConvNewTaskX" );
			// Set first value of stored output to be at "start"
			check_vsl( vslConvSetStart( task, start ), "vslConvSetStart" );

			// Compute convolutions between shape characteristic function and strain tensor components
			check_vsl( vslzConvExecX( task, strain_G.xx, ystride, convolutions_char_strain_region.char_xx, zstride ), "vslzConvExecX" );
			check_vsl( vslzConvExecX( task, strain_G.xy, ystride, convolutions_char_strain_region.char_xy, zstride ), "vslzConvExecX" );
			check_vsl( vslzConvExecX( task, strain_G.xz, ystride, convolutions_char_strain_region.char_xz, zstride ), "vslzConvExecX" );
			check_vsl( vslzConvExecX( task, strain_G.yy, ystride, convolutions_char_strain_region.char_yy, zstride ), "vslzConvExecX" );
			check_vsl( vslzConvExecX( task, strain_G.yz, ystride, convolutions_char_strain_region.char_yz, zstride ), "vslzConvExecX" );
			check_vsl( vslzConvExecX( task, strain_G.zz, ystride, convolutions_char_strain_region.char_zz, zstride ), "vslzConvExecX" );

			

			if ( piezo_order == 2 || piezo_order == 3 ) // If 2nd order piezoelectric potential included in calculation
			{

				// Compute convolutions between shape characteristic function and convolutions between strain tensor components
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.xx_yz, ystride, convolutions_char_strain_region.char_xx_yz, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.yy_yz, ystride, convolutions_char_strain_region.char_yy_yz, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.zz_yz, ystride, convolutions_char_strain_region.char_zz_yz, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.xy_xz, ystride, convolutions_char_strain_region.char_xy_xz, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.yy_xz, ystride, convolutions_char_strain_region.char_yy_xz, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.xx_xz, ystride, convolutions_char_strain_region.char_xx_xz, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.zz_xz, ystride, convolutions_char_strain_region.char_zz_xz, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.xy_yz, ystride, convolutions_char_strain_region.char_xy_yz, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.zz_xy, ystride, convolutions_char_strain_region.char_zz_xy, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.xx_xy, ystride, convolutions_char_strain_region.char_xx_xy, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.yy_xy, ystride, convolutions_char_strain_region.char_yy_xy, zstride ), "vslzConvExecX" );
				check_vsl( vslzConvExecX( task, convolutions_char_strain_region.xz_yz, ystride, convolutions_char_strain_region.char_xz_yz, zstride ), "vslzConvExecX" );

			}

			// Delete vsl task
			check_vsl( vslConvDeleteTask( &task ), "vslConvDeleteTask" );

		}

		#pragma omp parallel
		{

			// Loop over basis set
			#pragma omp for collapse( 3 )
			for( int idx_Gx = 0; idx_Gx < basis.N_Gx; idx_Gx++ )
			{

				for( int idx_Gy = 0; idx_Gy < basis.N_Gy; idx_Gy++ )
				{

					for( int idx_Gz = 0; idx_Gz < basis.N_Gz; idx_Gz++ )
					{

						// Set Fourier coefficient (plane wave) index

						int idx_basis = idx_Gz + basis.N_Gz*idx_Gy + basis.N_Gy*basis.N_Gz*idx_Gx;

						// Compute Fourier coefficient of the piezoelectric potential for the specified wave vector

						if( piezo_order == 1 ) // Potential due to first order piezoelectric polarisation only
						{

							piezo_G[idx_basis] += piezo_fourier_coefficient( 1, idx_shape, idx_Gx, idx_Gy, idx_Gz, basis, supercell, parameters, characteristic_function_G, strain_G, convolutions_char_strain_region );

						}
						else if( piezo_order == 2 ) // Potential due to first and second order piezoelectric polarisation
						{

							piezo_G[idx_basis] += piezo_fourier_coefficient( 1, idx_shape, idx_Gx, idx_Gy, idx_Gz, basis, supercell, parameters, characteristic_function_G, strain_G, convolutions_char_strain_region )
											+ piezo_fourier_coefficient( 2, idx_shape, idx_Gx, idx_Gy, idx_Gz, basis, supercell, parameters, characteristic_function_G, strain_G, convolutions_char_strain_region );

						}
						else if( piezo_order == 3 ) // Potential due to second order piezoelectric polarisation only
						{

							piezo_G[idx_basis] += piezo_fourier_coefficient( 2, idx_shape, idx_Gx, idx_Gy, idx_Gz, basis, supercell, parameters, characteristic_function_G, strain_G, convolutions_char_strain_region );

						}

					}

				}

			}

		}

	}

	
	free( convolutions_char_strain_region.char_xx );
	free( convolutions_char_strain_region.char_yy );
	free( convolutions_char_strain_region.char_zz );
	free( convolutions_char_strain_region.char_yz );
	free( convolutions_char_strain_region.char_xz );
	free( convolutions_char_strain_region.char_xy );
	
	if ( piezo_order == 2 || piezo_order == 3 )
	{

		free( convolutions_char_strain_region.xx_yz );
		free( convolutions_char_strain_region.yy_yz );
		free( convolutions_char_strain_region.zz_yz );
		free( convolutions_char_strain_region.xy_xz );
		free( convolutions_char_strain_region.yy_xz );
		free( convolutions_char_strain_region.xx_xz );
		free( convolutions_char_strain_region.zz_xz );
		free( convolutions_char_strain_region.xy_yz );
		free( convolutions_char_strain_region.zz_xy );
		free( convolutions_char_strain_region.xx_xy );
		free( convolutions_char_strain_region.yy_xy );
		free( convolutions_char_strain_region.xz_yz );

		free( convolutions_char_strain_region.char_xx_yz );
		free( convolutions_char_strain_region.char_yy_yz );
		free( convolutions_char_strain_region.char_zz_yz );
		free( convolutions_char_strain_region.char_xy_xz );
		free( convolutions_char_strain_region.char_yy_xz );
		free( convolutions_char_strain_region.char_xx_xz );
		free( convolutions_char_strain_region.char_zz_xz );
		free( convolutions_char_strain_region.char_xy_yz );
		free( convolutions_char_strain_region.char_zz_xy );
		free( convolutions_char_strain_region.char_xx_xy );
		free( convolutions_char_strain_region.char_yy_xy );
		free( convolutions_char_strain_region.char_xz_yz );
		
	}

}
