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

void convolutions_strain_correction( int idx_shape, struct wave_vectors basis, struct characteristic_functions *characteristic_function_strain_G, struct strain_tensor *strain_correction_region_G, struct convolution_char_strain *convolutions_char_strain_correction )
{
	
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


	// Create new vsl task (takes array as argument that will be involved in multiple convolution calculations)
	check_vsl( vslzConvNewTaskX( &task, VSL_CONV_MODE, 3, xshape, yshape, zshape, characteristic_function_strain_G[idx_shape].chi, xstride ), "vslzConvNewTaskX" );
	// Set first value of stored output to be at "start"
	check_vsl( vslConvSetStart( task, start ), "vslConvSetStart" );

	// Compute convolutions
	check_vsl( vslzConvExecX( task, strain_correction_region_G->xx, ystride, convolutions_char_strain_correction->char_xx, zstride ), "vslzConvExecX" );
	check_vsl( vslzConvExecX( task, strain_correction_region_G->xy, ystride, convolutions_char_strain_correction->char_xy, zstride ), "vslzConvExecX" );
	check_vsl( vslzConvExecX( task, strain_correction_region_G->xz, ystride, convolutions_char_strain_correction->char_xz, zstride ), "vslzConvExecX" );
	check_vsl( vslzConvExecX( task, strain_correction_region_G->yy, ystride, convolutions_char_strain_correction->char_yy, zstride ), "vslzConvExecX" );
	check_vsl( vslzConvExecX( task, strain_correction_region_G->yz, ystride, convolutions_char_strain_correction->char_yz, zstride ), "vslzConvExecX" );
	check_vsl( vslzConvExecX( task, strain_correction_region_G->zz, ystride, convolutions_char_strain_correction->char_zz, zstride ), "vslzConvExecX" );

	// Delete vsl task
	check_vsl( vslConvDeleteTask( &task ), "vslConvDeleteTask" );

}
