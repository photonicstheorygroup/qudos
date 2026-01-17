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

//////////////
//// E_CB ////
//////////////

complex double dE_CB_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result
	
	// Row wave vector index
	int idx_G_1_basis[3];
	// Column wave vector index
	int idx_G_2_basis[3];
	
	idx_G_1_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_1[0];
	idx_G_1_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_1[1];
	idx_G_1_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_1[2];
	
	idx_G_2_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_2[0];
	idx_G_2_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_2[1];
	idx_G_2_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_2[2];


	// Parameter arrays
	double *sc   = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		  sc[idx_shape] = parameters[idx_shape].sc;

	}
	

	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result = 0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, sc )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx );

			break;

		case 1: // y polarisation

			result = 0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, sc )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky );

			break;

		case 2: // z polarisation

			result = 0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, sc )*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz );

			break;

	}


	// Free dynamically allocated memory
	free( sc );


	return result;

}


//////////////
//// E_HH ////
//////////////

complex double dE_HH_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result
	
	// Row wave vector index
	int idx_G_1_basis[3];
	// Column wave vector index
	int idx_G_2_basis[3];
	
	idx_G_1_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_1[0];
	idx_G_1_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_1[1];
	idx_G_1_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_1[2];
	
	idx_G_2_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_2[0];
	idx_G_2_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_2[1];
	idx_G_2_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_2[2];


	// Parameter arrays
	double *g1_plus_g2    = ( double* ) malloc( N_shape*sizeof( double ) );
	double *g1_minus_2_g2 = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		g1_plus_g2   [idx_shape] = parameters[idx_shape].g1 +     parameters[idx_shape].g2;
		g1_minus_2_g2[idx_shape] = parameters[idx_shape].g1 - 2.0*parameters[idx_shape].g2;

	}


	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result = -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1_plus_g2    )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx );

			break;

		case 1: // y polarisation

			result = -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1_plus_g2    )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky );

			break;

		case 2: // z polarisation

			result = -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1_minus_2_g2 )*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz );

			break;

	}


	// Free dynamically allocated memory
	free( g1_plus_g2    );
	free( g1_minus_2_g2 );


	return result;

}


//////////////
//// E_LH ////
//////////////

complex double dE_LH_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result
	
	// Row wave vector index
	int idx_G_1_basis[3];
	// Column wave vector index
	int idx_G_2_basis[3];
	
	idx_G_1_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_1[0];
	idx_G_1_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_1[1];
	idx_G_1_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_1[2];
	
	idx_G_2_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_2[0];
	idx_G_2_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_2[1];
	idx_G_2_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_2[2];


	// Parameter arrays
	double *g1_minus_g2  = ( double* ) malloc( N_shape*sizeof( double ) );
	double *g1_plus_2_g2 = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		g1_minus_g2 [idx_shape] = parameters[idx_shape].g1 -     parameters[idx_shape].g2;
		g1_plus_2_g2[idx_shape] = parameters[idx_shape].g1 + 2.0*parameters[idx_shape].g2;

	}


	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result = -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1_minus_g2  )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx );

			break;

		case 1: // y polarisation

			result = -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1_minus_g2  )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky );

			break;

		case 2: // z polarisation

			result = -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1_plus_2_g2 )*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz );

			break;

	}


	// Free dynamically allocated memory
	free( g1_minus_g2  );
	free( g1_plus_2_g2 );


	return result;

}


//////////////
//// E_SO ////
//////////////

complex double dE_SO_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result
	
	// Row wave vector index
	int idx_G_1_basis[3];
	// Column wave vector index
	int idx_G_2_basis[3];
	
	idx_G_1_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_1[0];
	idx_G_1_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_1[1];
	idx_G_1_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_1[2];
	
	idx_G_2_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_2[0];
	idx_G_2_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_2[1];
	idx_G_2_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_2[2];


	// Parameter arrays
	double *g1 = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		g1[idx_shape] = parameters[idx_shape].g1;

	}


	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result = -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1 )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx );

			break;

		case 1: // y polarisation

			result = -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1 )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky );

			break;

		case 2: // z polarisation

			result = -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1 )*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz );

			break;

	}


	// Free dynamically allocated memory
	free( g1 );


	return result;

}


///////////
//// Q ////
///////////

complex double dQ_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result
	
	// Row wave vector index
	int idx_G_1_basis[3];
	// Column wave vector index
	int idx_G_2_basis[3];
	
	idx_G_1_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_1[0];
	idx_G_1_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_1[1];
	idx_G_1_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_1[2];
	
	idx_G_2_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_2[0];
	idx_G_2_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_2[1];
	idx_G_2_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_2[2];


	// Parameter arrays
	double *g2 = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		g2[idx_shape] = parameters[idx_shape].g2;

	}


	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result = -1.0*sqrt( 0.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx );

			break;

		case 1: // y polarisation

			result = -1.0*sqrt( 0.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky );

			break;

		case 2: // z polarisation

			result =      sqrt( 2.0 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 )*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz );

			break;

	}


	// Free dynamically allocated memory
	free( g2 );


	return result;

}


///////////
//// R ////
///////////

complex double dR_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result
	
	// Row wave vector index
	int idx_G_1_basis[3];
	// Column wave vector index
	int idx_G_2_basis[3];
	
	idx_G_1_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_1[0];
	idx_G_1_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_1[1];
	idx_G_1_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_1[2];
	
	idx_G_2_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_2[0];
	idx_G_2_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_2[1];
	idx_G_2_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_2[2];


	// Parameter arrays
	double *g2 = ( double* ) malloc( N_shape*sizeof( double ) );
	double *g3 = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		g2[idx_shape] = parameters[idx_shape].g2;
		g3[idx_shape] = parameters[idx_shape].g3;

	}


	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result =  0.5*sqrt( 3.0 )*HBSQM*( parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx )
		 							   - 1.0*I*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) );

			break;

		case 1: // y polarisation

			result = -0.5*sqrt( 3.0 )*HBSQM*( parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky )
									   +     I*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx ) );

			break;

		case 2: // z polarisation

			result = 0.0 + I*0.0;

			break;

	}


	// Free dynamically allocated memory
	free( g2 );
	free( g3 );


	return result;

}


////////////////
//// R_conj ////
////////////////

complex double dR_conj_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result
	
	// Row wave vector index
	int idx_G_1_basis[3];
	// Column wave vector index
	int idx_G_2_basis[3];
	
	idx_G_1_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_1[0];
	idx_G_1_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_1[1];
	idx_G_1_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_1[2];
	
	idx_G_2_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_2[0];
	idx_G_2_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_2[1];
	idx_G_2_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_2[2];


	// Parameter arrays
	double *g2 = ( double* ) malloc( N_shape*sizeof( double ) );
	double *g3 = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		g2[idx_shape] = parameters[idx_shape].g2;
		g3[idx_shape] = parameters[idx_shape].g3;

	}


	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result =  0.5*sqrt( 3.0 )*HBSQM*( parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx )
		 							   + 1.0*I*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) );

			break;

		case 1: // y polarisation

			result = -0.5*sqrt( 3.0 )*HBSQM*( parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky )
									   -     I*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx ) );

			break;

		case 2: // z polarisation

			result = 0.0 + I*0.0;

			break;

	}


	// Free dynamically allocated memory
	free( g2 );
	free( g3 );


	return result;

}


///////////
//// S ////
///////////

complex double dS_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result
	
	// Row wave vector index
	int idx_G_1_basis[3];
	// Column wave vector index
	int idx_G_2_basis[3];
	
	idx_G_1_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_1[0];
	idx_G_1_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_1[1];
	idx_G_1_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_1[2];
	
	idx_G_2_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_2[0];
	idx_G_2_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_2[1];
	idx_G_2_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_2[2];


	// Parameter arrays
	double *g3 = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		g3[idx_shape] = parameters[idx_shape].g3;

	}


	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result =  0.5  *sqrt( 1.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz );

			break;

		case 1: // y polarisation

			result = -0.5*I*sqrt( 1.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz );

			break;

		case 2: // z polarisation

			result = 0.5*  sqrt( 1.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*(   ( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx )
																								    - I*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) );

			break;

	}


	// Free dynamically allocated memory
	free( g3 );


	return result;

}


////////////////
//// S_conj ////
////////////////

complex double dS_conj_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result
	
	// Row wave vector index
	int idx_G_1_basis[3];
	// Column wave vector index
	int idx_G_2_basis[3];
	
	idx_G_1_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_1[0];
	idx_G_1_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_1[1];
	idx_G_1_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_1[2];
	
	idx_G_2_basis[0] = basis.N_Gx_T/2 - basis.N_Gx/2 + idx_G_2[0];
	idx_G_2_basis[1] = basis.N_Gy_T/2 - basis.N_Gy/2 + idx_G_2[1];
	idx_G_2_basis[2] = basis.N_Gz_T/2 - basis.N_Gz/2 + idx_G_2[2];


	// Parameter arrays
	double *g3 = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		g3[idx_shape] = parameters[idx_shape].g3;

	}


	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result =  0.5  *sqrt( 1.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz );

			break;

		case 1: // y polarisation

			result =  0.5*I*sqrt( 1.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz );

			break;

		case 2: // z polarisation

			result = 0.5*  sqrt( 1.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*(   ( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx )
																								    + I*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) );

			break;

	}


	// Free dynamically allocated memory
	free( g3 );


	return result;

}


////////////////
//// T_plus ////
////////////////

complex double dT_plus_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	
	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result


	// Parameter arrays
	double *P = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		P[idx_shape] = parameters[idx_shape].P;

	}


	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result =     parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, P )  /sqrt( 6.0 );
			
			// Strain
			if( include_strain != 0 )
			{
	
				// Wave vector index
				int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
				 // Matrix term
				result += ( -1/sqrt( 6 ) )*parameters[0].P*( ( strain_G.xx[idx_basis] + I*strain_G.xy[idx_basis] ) );
				
				// Shape terms
				for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
				{

					result += ( -1/sqrt( 6 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*( convolutions_char_strain[idx_shape].char_xx[idx_basis] + I*convolutions_char_strain[idx_shape].char_xy[idx_basis] );

				}

			}

			break;

		case 1: // y polarisation

			result = ( I*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, P ) )/sqrt( 6.0 );

			// Strain
			if( include_strain != 0 )
			{
	
				// Wave vector index
				int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
				// Matrix term
				result += ( -1/sqrt( 6 ) )*parameters[0].P*( ( strain_G.xy[idx_basis] + I*strain_G.yy[idx_basis] ) ); 
				
				// Shape terms
				for( idx_shape = 1; idx_shape < N_shape; idx_shape++ ) 
				{

					result += ( -1/sqrt( 6 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*( convolutions_char_strain[idx_shape].char_xy[idx_basis] + I*convolutions_char_strain[idx_shape].char_yy[idx_basis] );

				}

			}

			break;

		case 2: // z polarisation

			result = 0.0 + I*0.0;

			// Strain
			if( include_strain != 0 )
			{
	
				// Wave vector index
				int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
				// Matrix term
				result += ( -1/sqrt( 6 ) )*parameters[0].P*( ( strain_G.xz[idx_basis] + I*strain_G.yz[idx_basis] ) );
				
				// Shape terms
				for( idx_shape = 1; idx_shape < N_shape; idx_shape++ ) 
				{

					result += ( -1/sqrt( 6 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*( convolutions_char_strain[idx_shape].char_xz[idx_basis] + I*convolutions_char_strain[idx_shape].char_yz[idx_basis] );

				}

			}

			break;

	}


	// Free dynamically allocated memory
	free( P );


	return result;

}


/////////////////
//// T_minus ////
/////////////////

complex double dT_minus_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result = 0.0 + I*0.0; // Will store overall result


	// Parameter arrays
	double *P = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		P[idx_shape] = parameters[idx_shape].P;

	}


	// Evaluate matrix element
	switch( polarisation )
	{

		case 0: // x polarisation

			result =     parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, P )  /sqrt( 6.0 );
			
			// Strain
			if( include_strain != 0 )
			{
	
				// Wave vector index
				int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
				// Matrix term
				result += ( -1/sqrt( 6 ) )*parameters[0].P*( ( strain_G.xx[idx_basis] - I*strain_G.xy[idx_basis] ) ); 
				
				// Shape terms
				for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
				{

					result += ( -1/sqrt( 6 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*( convolutions_char_strain[idx_shape].char_xx[idx_basis] - I*convolutions_char_strain[idx_shape].char_xy[idx_basis] );

				}

			}

			break;

		case 1: // y polarisation

			result = ( -I*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, P ) )/sqrt( 6.0 );

			// Strain
			if( include_strain != 0 )
			{
	
				// Wave vector index
				int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
				// Matrix term
				result += ( -1/sqrt( 6 ) )*parameters[0].P*( ( strain_G.xy[idx_basis] - I*strain_G.yy[idx_basis] ) ); 
				
				// Shape terms
				for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
				{

					result += ( -1/sqrt( 6 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*( convolutions_char_strain[idx_shape].char_xy[idx_basis] - I*convolutions_char_strain[idx_shape].char_yy[idx_basis] );

				}

			}

			break;

		case 2: // z polarisation

			result = 0.0 + I*0.0;

			// Strain
			if( include_strain != 0 )
			{
	
				// Wave vector index
				int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
				// Matrix term
				result += ( -1/sqrt( 6 ) )*parameters[0].P*( ( strain_G.xz[idx_basis] - I*strain_G.yz[idx_basis] ) );
				
				// Shape terms
				for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
				{

					result += ( -1/sqrt( 6 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*( convolutions_char_strain[idx_shape].char_xz[idx_basis] - I*convolutions_char_strain[idx_shape].char_yz[idx_basis] );

				}

			}

			break;

	}


	// Free dynamically allocated memory
	free( P );

	return result;

}


///////////
//// U ////
///////////

complex double dU_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{
	
	int idx_shape;
	complex double result = 0.0 + I*0.0;


	// Parameter arrays
	double *P = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		P[idx_shape] = parameters[idx_shape].P;

	}


	// Evaluate matrix element

	switch( polarisation )
	{

		case 0: // x polarisation

			result = 0.0 + I*0.0;
			
			// Strain
			if( include_strain != 0 )
			{
	
				// Wave vector index
				int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
				// Matrix term
				result += ( -1/sqrt( 3 ) )*parameters[0].P*( ( strain_G.xz[idx_basis] ) );
				
				// Shape terms
				for( idx_shape = 1; idx_shape < N_shape; idx_shape++ ) 
				{

					result += ( -1/sqrt( 3 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*convolutions_char_strain[idx_shape].char_xz[idx_basis];

				}

			}

			break;

		case 1: // y polarisation

			result = 0.0 + I*0.0;

			// Strain


			if( include_strain != 0 )
			{
	
				// Wave vector index
				int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
				// Matrix term
				result += ( -1/sqrt( 3 ) )*parameters[0].P*( ( strain_G.yz[idx_basis] ) );
				
				// Shape terms
				for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
				{

					result += ( -1/sqrt( 3 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*convolutions_char_strain[idx_shape].char_yz[idx_basis];

				}

			}

			break;

		case 2: // z polarisation

			result = parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, P )/sqrt( 3.0 );

			// Strain


			if( include_strain != 0 )
			{
	
				// Wave vector index
				int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
				// Matrix term
				result += ( -1/sqrt( 3 ) )*parameters[0].P*( ( strain_G.zz[idx_basis] ) );
				
				// Shape terms
				for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
				{

					result += ( -1/sqrt( 3 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*convolutions_char_strain[idx_shape].char_zz[idx_basis];

				}

			}

			break;

	}


	// Free dynamically allocated memory
	free( P );


	return result;

}
