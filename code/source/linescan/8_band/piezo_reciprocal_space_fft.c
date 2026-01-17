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

void piezo_reciprocal_space_fft( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell, complex double *phi_piezo_G, complex double *phi_piezo_region_G, struct characteristic_functions *characteristic_function_strain_G )
{

	VSLConvTaskPtr task; // vsl task pointer

	MKL_INT VSL_CONV_MODE = VSL_CONV_MODE_FFT; // Using FFT mode of vsl to compute linear convolution

	// Sizes of dimensions of two input and one output array (output array sized to only extract central portion of output)
	const MKL_INT xshape[] = { strain_basis.N_Gx, strain_basis.N_Gy, strain_basis.N_Gz };
	const MKL_INT yshape[] = { strain_basis.N_Gx, strain_basis.N_Gy, strain_basis.N_Gz };
	const MKL_INT zshape[] = { basis.N_Gx, basis.N_Gy, basis.N_Gz };

	// Where to start taking data from output array (set to only extract central portion)
	const MKL_INT start[] = { strain_basis.N_Gx - 1 - supercell[0].RL1*( basis.N_Gx / 2 ), strain_basis.N_Gy - 1 - supercell[0].RL2*( basis.N_Gy / 2 ), strain_basis.N_Gz - 1 - supercell[0].RL3*( basis.N_Gz / 2 ) };

	// Strides of stored data (standard data layout)
	const MKL_INT xstride[] = { ( strain_basis.N_Gy )*( strain_basis.N_Gz ), strain_basis.N_Gz, 1 };
	const MKL_INT ystride[] = { ( strain_basis.N_Gy )*( strain_basis.N_Gz ), strain_basis.N_Gz, 1 };
	const MKL_INT zstride[] = { ( basis.N_Gy )*( basis.N_Gz ), basis.N_Gz, 1 };

	// Every (supercell[0].RL(1/2/3))th value of the output array is extracted 
	const MKL_INT decimation[] = { supercell[0].RL1, supercell[0].RL2, supercell[0].RL3 };

	// Create new vsl task (takes array as argument that will be involved in multiple convolution calculations)
	check_vsl( vslzConvNewTaskX( &task, VSL_CONV_MODE, 3, xshape, yshape, zshape, characteristic_function_strain_G[0].chi, xstride ), "vslzConvNewTaskX" );
	// Set first value of stored output to be at "start"
	check_vsl( vslConvSetStart( task, start ), "vslConvSetStart" );
	// Set decimation of output array
	check_vsl( vslConvSetDecimation( task, decimation ), "vslConvSetDecimation" );

	// Compute convolution
	check_vsl( vslzConvExecX( task, phi_piezo_region_G, ystride, phi_piezo_G, zstride ), "vslzConvExecX" );

	// Delete vsl task
	check_vsl( vslConvDeleteTask( &task ), "vslConvDeleteTask" );

}
