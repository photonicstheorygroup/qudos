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

complex double E_CB( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


	// Parameter arrays
	double *E_cb = ( double* ) malloc( N_shape*sizeof( double ) );
	double *sc   = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		E_cb[idx_shape] = parameters[idx_shape].E_cb;
		  sc[idx_shape] = parameters[idx_shape].sc  ;

	}
	
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

	// Band edge energy
	result = parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, E_cb );

	// Kinetic
	result += 0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, sc )*( basis.Gx[idx_G_1_basis[0]]*basis.Gx[idx_G_2_basis[0]] + basis.Kx*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] )
																				 + basis.Gy[idx_G_1_basis[1]]*basis.Gy[idx_G_2_basis[1]] + basis.Ky*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] )
																				 + basis.Gz[idx_G_1_basis[2]]*basis.Gz[idx_G_2_basis[2]] + basis.Kz*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] )
																				 + basis.Kx*basis.Kx + basis.Ky*basis.Ky + basis.Kz*basis.Kz );


	// Strain
	if( include_strain != 0 )
	{
		
		// Wave vector index
		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
		// Matrix term
		result += parameters[0].ac*( strain_G.xx[idx_basis]
							  + strain_G.yy[idx_basis]
							  + strain_G.zz[idx_basis] );
							  
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ ) // Shape terms
		{

			result += ( parameters[idx_shape].ac - parameters[supercell[idx_shape].embedded].ac )*convolutions_char_strain[idx_shape].char_hydrostatic[idx_basis];


		}

		// Piezoelectric potential
		if( include_piezo != 0 )
		{

			result += phi_piezo_G.pz[idx_basis];

		}

	}


	// Applied electric field
	result += field_fourier_coefficient( F, idx_G_1_basis, idx_G_2_basis, basis );


	// Free dynamically allocated memory
	free( E_cb );
	free( sc   );


	return result;

}


//////////////
//// E_HH ////
//////////////

complex double E_HH( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


	// Parameter arrays
	double *E_hh          = ( double* ) malloc( N_shape*sizeof( double ) );
	double *g1_plus_g2    = ( double* ) malloc( N_shape*sizeof( double ) );
	double *g1_minus_2_g2 = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		E_hh         [idx_shape] = parameters[idx_shape].E_hh;
		g1_plus_g2   [idx_shape] = parameters[idx_shape].g1 +     parameters[idx_shape].g2;
		g1_minus_2_g2[idx_shape] = parameters[idx_shape].g1 - 2.0*parameters[idx_shape].g2;

	}
	
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


	// Band edge energy
	result = parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, E_hh );


	// Kinetic
	result += -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1_plus_g2    )*( basis.Gx[idx_G_1_basis[0]]*basis.Gx[idx_G_2_basis[0]] + basis.Kx*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] ) + basis.Gy[idx_G_1_basis[1]]*basis.Gy[idx_G_2_basis[1]] + basis.Ky*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] ) + basis.Kx*basis.Kx + basis.Ky*basis.Ky )
		   -  0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1_minus_2_g2 )*( basis.Gz[idx_G_1_basis[2]]*basis.Gz[idx_G_2_basis[2]] + basis.Kz*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] ) + basis.Kz*basis.Kz );


	// Strain
	if( include_strain != 0 )
	{
		
		// Wave vector index
		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
		// Matrix term
		result +=     parameters[0].av*( strain_G.xx[idx_basis]
							  + strain_G.yy[idx_basis]
							  + strain_G.zz[idx_basis] )
			   + 0.5*parameters[0].b *( strain_G.xx[idx_basis] + strain_G.yy[idx_basis] - 2.0*strain_G.zz[idx_basis] );
			   
			   
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ ) // Shape terms
		{

			result += ( parameters[idx_shape].av - parameters[supercell[idx_shape].embedded].av )*convolutions_char_strain[idx_shape].char_hydrostatic[idx_basis] + 0.5*( parameters[idx_shape].b - parameters[supercell[idx_shape].embedded].b )*convolutions_char_strain[idx_shape].char_biaxial[idx_basis];


		}

		// Piezoelectric potential
		if( include_piezo != 0 )
		{

			result += phi_piezo_G.pz[idx_basis];

		}

	}


	// Applied electric field
	result += field_fourier_coefficient( F, idx_G_1_basis, idx_G_2_basis, basis );


	// Free dynamically allocated memory
	free( E_hh          );
	free( g1_plus_g2    );
	free( g1_minus_2_g2 );


	return result;

}


//////////////
//// E_LH ////
//////////////

complex double E_LH( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


	// Parameter arrays
	double *E_lh         = ( double* ) malloc( N_shape*sizeof( double ) );
	double *g1_minus_g2  = ( double* ) malloc( N_shape*sizeof( double ) );
	double *g1_plus_2_g2 = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		E_lh        [idx_shape] = parameters[idx_shape].E_lh;
		g1_minus_g2 [idx_shape] = parameters[idx_shape].g1 -     parameters[idx_shape].g2;
		g1_plus_2_g2[idx_shape] = parameters[idx_shape].g1 + 2.0*parameters[idx_shape].g2;

	}
	
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


	// Band edge energy
	result = parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, E_lh );


	// Kinetic
	result += -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1_minus_g2  )*( basis.Gx[idx_G_1_basis[0]]*basis.Gx[idx_G_2_basis[0]] + basis.Kx*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] ) + basis.Gy[idx_G_1_basis[1]]*basis.Gy[idx_G_2_basis[1]] + basis.Ky*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] ) + basis.Kx*basis.Kx + basis.Ky*basis.Ky )
		   -  0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1_plus_2_g2 )*( basis.Gz[idx_G_1_basis[2]]*basis.Gz[idx_G_2_basis[2]] + basis.Kz*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] ) + basis.Kz*basis.Kz );


	// Strain
	if( include_strain != 0 )
	{

		// Wave vector index
		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
		// Matrix term
		result +=     parameters[0].av*( strain_G.xx[idx_basis]
							  + strain_G.yy[idx_basis]
							  + strain_G.zz[idx_basis] )
			   - 0.5*parameters[0].b *( strain_G.xx[idx_basis] + strain_G.yy[idx_basis] - 2.0*strain_G.zz[idx_basis] ); 
			   
		// Shape terms
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			result += ( parameters[idx_shape].av - parameters[supercell[idx_shape].embedded].av )*convolutions_char_strain[idx_shape].char_hydrostatic[idx_basis] - 0.5*( parameters[idx_shape].b - parameters[supercell[idx_shape].embedded].b )*convolutions_char_strain[idx_shape].char_biaxial[idx_basis];


		}

		// Piezoelectric potential
		if( include_piezo != 0 )
		{

			result += phi_piezo_G.pz[idx_basis];

		}

	}


	// Applied electric field
	result += field_fourier_coefficient( F, idx_G_1_basis, idx_G_2_basis, basis );


	// Free dynamically allocated memory
	free( E_lh         );
	free( g1_minus_g2  );
	free( g1_plus_2_g2 );


	return result;

}


//////////////
//// E_SO ////
//////////////

complex double E_SO( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


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
	double *E_so = ( double* ) malloc( N_shape*sizeof( double ) );
	double *g1   = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		E_so[idx_shape] = parameters[idx_shape].E_so;
		g1  [idx_shape] = parameters[idx_shape].g1;

	}


	// Band edge energy
	result = parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, E_so );


	// Kinetic
	result += -0.5*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g1 )*( basis.Gx[idx_G_1_basis[0]]*basis.Gx[idx_G_2_basis[0]] + basis.Kx*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] )
																									+ basis.Gy[idx_G_1_basis[1]]*basis.Gy[idx_G_2_basis[1]] + basis.Ky*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] )
																									+ basis.Gz[idx_G_1_basis[2]]*basis.Gz[idx_G_2_basis[2]] + basis.Kz*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] )
																									+ basis.Kx*basis.Kx + basis.Ky*basis.Ky + basis.Kz*basis.Kz );


	// Strain
  	if( include_strain != 0 )
  	{

		// Wave vector index
  		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
		// Matrix term
		result +=     parameters[0].av*( strain_G.xx[idx_basis]
							  + strain_G.yy[idx_basis]
							  + strain_G.zz[idx_basis] );
							  
		// Shape terms
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			result += ( parameters[idx_shape].av - parameters[supercell[idx_shape].embedded].av )*convolutions_char_strain[idx_shape].char_hydrostatic[idx_basis];

		}

		// Piezoelectric potential
		if( include_piezo != 0 )
		{

			result += phi_piezo_G.pz[idx_basis];

		}

  	}


	// Applied electric field
	result += field_fourier_coefficient( F, idx_G_1_basis, idx_G_2_basis, basis );


	// Free dynamically allocated memory
	free( E_so );
	free( g1   );


	return result;

}


///////////
//// Q ////
///////////

complex double Q( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


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


	// Kinetic
	result = -1.0*sqrt( 0.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 )*( basis.Gx[idx_G_1_basis[0]]*basis.Gx[idx_G_2_basis[0]] + basis.Kx*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] ) + basis.Gy[idx_G_1_basis[1]]*basis.Gy[idx_G_2_basis[1]] + basis.Ky*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] ) + basis.Kx*basis.Kx + basis.Ky*basis.Ky )
		  +      sqrt( 2.0 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 )*( basis.Gz[idx_G_1_basis[2]]*basis.Gz[idx_G_2_basis[2]] + basis.Kz*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] ) + basis.Kz*basis.Kz );


	// Strain
	if( include_strain != 0 )
	{

		// Wave vector index
		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
		// Matrix term
		result +=  sqrt(0.5)*parameters[0].b *( strain_G.xx[idx_basis] + strain_G.yy[idx_basis] - 2.0*strain_G.zz[idx_basis] );
		
		// Shape terms
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			result += sqrt(0.5)*( parameters[idx_shape].b - parameters[supercell[idx_shape].embedded].b )*convolutions_char_strain[idx_shape].char_biaxial[idx_basis];

		}

	}


	// Free dynamically allocated memory
	free( g2 );


	return result;

}


///////////
//// R ////
///////////

complex double R( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


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


	// Kinetic
	result =   0.5*sqrt( 3.0 )*HBSQM*(parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 ))*( basis.Gx[idx_G_1_basis[0]]*basis.Gx[idx_G_2_basis[0]] + basis.Kx*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] ) - basis.Gy[idx_G_1_basis[1]]*basis.Gy[idx_G_2_basis[1]] - basis.Ky*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] ) + basis.Kx*basis.Kx - basis.Ky*basis.Ky )
		  - I*0.5*sqrt( 3.0 )*HBSQM*(parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 ))*( ( basis.Gx[idx_G_1_basis[0]] + basis.Kx )*( basis.Gy[idx_G_2_basis[1]] + basis.Ky ) + ( basis.Gy[idx_G_1_basis[1]] + basis.Ky )*( basis.Gx[idx_G_2_basis[0]] + basis.Kx ) );


	// Strain
	if( include_strain != 0 )
	{
		
		// Wave vector index
		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
		// Matrix term
		result +=  -0.5*sqrt(3.0)*parameters[0].b*( strain_G.xx[idx_basis] - strain_G.yy[idx_basis]) +                I*parameters[0].d*strain_G.xy[idx_basis];

		// Shape terms
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			result += -0.5*sqrt(3.0)*( parameters[idx_shape].b - parameters[supercell[idx_shape].embedded].b )*( convolutions_char_strain[idx_shape].char_xx[idx_basis] - convolutions_char_strain[idx_shape].char_yy[idx_basis] ) + I*( parameters[idx_shape].d - parameters[supercell[idx_shape].embedded].d )*convolutions_char_strain[idx_shape].char_xy[idx_basis];

		}

	}


	// Free dynamically allocated memory
	free( g2 );
	free( g3 );


	return result;

}


////////////////
//// R_conj ////
////////////////

complex double R_conj( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


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


	// Kinetic
	result =   0.5*sqrt( 3.0 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g2 )*( basis.Gx[idx_G_1_basis[0]]*basis.Gx[idx_G_2_basis[0]] + basis.Kx*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] ) - basis.Gy[idx_G_1_basis[1]]*basis.Gy[idx_G_2_basis[1]] - basis.Ky*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] ) + basis.Kx*basis.Kx - basis.Ky*basis.Ky )
		  + I*0.5*sqrt( 3.0 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( ( basis.Gx[idx_G_1_basis[0]] + basis.Kx )*( basis.Gy[idx_G_2_basis[1]] + basis.Ky ) + ( basis.Gy[idx_G_1_basis[1]] + basis.Ky )*( basis.Gx[idx_G_2_basis[0]] + basis.Kx ) );


	// Strain
	if( include_strain != 0 )
	{

		// Wave vector index
		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);
		
		// Matrix term
		result +=  -0.5*sqrt(3.0)*parameters[0].b*( strain_G.xx[idx_basis] - strain_G.yy[idx_basis]) -                I*parameters[0].d*strain_G.xy[idx_basis];
		
		// Shape terms
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			result += -0.5*sqrt( 3.0)*( parameters[idx_shape].b - parameters[supercell[idx_shape].embedded].b )*( convolutions_char_strain[idx_shape].char_xx[idx_basis] - convolutions_char_strain[idx_shape].char_yy[idx_basis] ) - I*( parameters[idx_shape].d - parameters[supercell[idx_shape].embedded].d )*convolutions_char_strain[idx_shape].char_xy[idx_basis];

		}

	}


	// Free dynamically allocated memory
	free( g2 );
	free( g3 );


	return result;

}


///////////
//// S ////
///////////

complex double S( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


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


	// Kinetic
	result =   0.5*sqrt( 1.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( ( basis.Gx[idx_G_1_basis[0]] + basis.Kx )*( basis.Gz[idx_G_2_basis[2]] + basis.Kz ) + ( basis.Gz[idx_G_1_basis[2]] + basis.Kz )*( basis.Gx[idx_G_2_basis[0]] + basis.Kx ) )
		  - I*0.5*sqrt( 1.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( ( basis.Gy[idx_G_1_basis[1]] + basis.Ky )*( basis.Gz[idx_G_2_basis[2]] + basis.Kz ) + ( basis.Gz[idx_G_1_basis[2]] + basis.Kz )*( basis.Gy[idx_G_2_basis[1]] + basis.Ky ) );


	// Strain
	if( include_strain != 0 )
	{
	
		// Wave vector index
		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);

		// Matrix term
		result += -1.0*sqrt( 0.5 )*parameters[0].d*( strain_G.xz[idx_basis] - I*strain_G.yz[idx_basis] );
		
		// Shape terms
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			result += -1.0*sqrt( 0.5 )*( parameters[idx_shape].d - parameters[supercell[idx_shape].embedded].d )*( convolutions_char_strain[idx_shape].char_xz[idx_basis] - I*convolutions_char_strain[idx_shape].char_yz[idx_basis] );

		}

	}


	// Free dynamically allocated memory
	free( g3 );


	return result;

}


////////////////
//// S_conj ////
////////////////

complex double S_conj( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


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


	// Kinetic
	result =   0.5*sqrt( 1.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( ( basis.Gx[idx_G_1_basis[0]] + basis.Kx )*( basis.Gz[idx_G_2_basis[2]] + basis.Kz ) + ( basis.Gz[idx_G_1_basis[2]] + basis.Kz )*( basis.Gx[idx_G_2_basis[0]] + basis.Kx ) )
		  + I*0.5*sqrt( 1.5 )*HBSQM*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, g3 )*( ( basis.Gy[idx_G_1_basis[1]] + basis.Ky )*( basis.Gz[idx_G_2_basis[2]] + basis.Kz ) + ( basis.Gz[idx_G_1_basis[2]] + basis.Kz )*( basis.Gy[idx_G_2_basis[1]] + basis.Ky ) );


	// Strain
  	if( include_strain != 0 )
  	{

		// Wave vector index
  		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);

		// Matrix term
		result += -1.0*sqrt( 0.5 )*parameters[0].d*( strain_G.xz[idx_basis] + I*strain_G.yz[idx_basis] );
		
		// Shape terms
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ ) 
		{

			result += -1.0*sqrt( 0.5 )*( parameters[idx_shape].d - parameters[supercell[idx_shape].embedded].d )*( convolutions_char_strain[idx_shape].char_xz[idx_basis] + I*convolutions_char_strain[idx_shape].char_yz[idx_basis] );

		}

  	}


	// Free dynamically allocated memory
	free( g3 );


	return result;

}


////////////////
//// T_plus ////
////////////////

complex double T_plus( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


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
	double *P = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		P[idx_shape] = parameters[idx_shape].P;

	}


	// Kinetic
	result = (   parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, P )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx ) )/( 2.0*sqrt( 6.0 ) )
	 	  + ( I*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, P )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) )/( 2.0*sqrt( 6.0 ) );


	// Strain
	if( include_strain != 0 )
	{
	
		// Wave vector index
		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);

		// Matrix term
		result += ( -0.5/sqrt( 6 ) )*parameters[0].P*( ( strain_G.xx[idx_basis] + I*strain_G.xy[idx_basis] ) * ( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx ) + ( strain_G.xy[idx_basis] + I*strain_G.yy[idx_basis] ) * ( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) + ( strain_G.xz[idx_basis] + I*strain_G.yz[idx_basis] ) * ( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz ) ); 

		// Shape terms
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			result += ( -0.5/sqrt( 6 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*( ( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx ) * ( convolutions_char_strain[idx_shape].char_xx[idx_basis] + I*convolutions_char_strain[idx_shape].char_xy[idx_basis] ) + ( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) * ( convolutions_char_strain[idx_shape].char_xy[idx_basis] + I*convolutions_char_strain[idx_shape].char_yy[idx_basis] ) + ( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz ) * ( convolutions_char_strain[idx_shape].char_xz[idx_basis] + I*convolutions_char_strain[idx_shape].char_yz[idx_basis] ) );

		}

	}

	// Free dynamically allocated memory
	free( P );


	return result;

}


/////////////////
//// T_minus ////
/////////////////

complex double T_minus( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


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
	double *P = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		P[idx_shape] = parameters[idx_shape].P;

	}


	// Kinetic
	result = (   parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, P )*( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx ) )/( 2.0*sqrt( 6.0 ) )
	 	  - ( I*parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, P )*( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) )/( 2.0*sqrt( 6.0 ) );
	 	  
	
	// Strain
	if( include_strain != 0 )
	{
	
		// Wave vector index
		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);

		// Matrix term
		result += ( -0.5/sqrt( 6 ) )*parameters[0].P*( ( strain_G.xx[idx_basis] - I*strain_G.xy[idx_basis] ) * ( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx ) + ( strain_G.xy[idx_basis] - I*strain_G.yy[idx_basis] ) * ( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) + ( strain_G.xz[idx_basis] - I*strain_G.yz[idx_basis] ) * ( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz ) ); 
		
		// Shape terms
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			result += ( -0.5/sqrt( 6 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*( ( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx ) * ( convolutions_char_strain[idx_shape].char_xx[idx_basis] - I*convolutions_char_strain[idx_shape].char_xy[idx_basis] ) + ( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) * ( convolutions_char_strain[idx_shape].char_xy[idx_basis] - I*convolutions_char_strain[idx_shape].char_yy[idx_basis] ) + ( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz ) * ( convolutions_char_strain[idx_shape].char_xz[idx_basis] - I*convolutions_char_strain[idx_shape].char_yz[idx_basis] ) );

		}

	}


	// Free dynamically allocated memory
	free( P );


	return result;

}


///////////
//// U ////
///////////

complex double U( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
{

	int idx_shape; // Shape index
	complex double result; // Will store overall result


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
	double *P = ( double* ) malloc( N_shape*sizeof( double ) );

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		P[idx_shape] = parameters[idx_shape].P;

	}


	// Kinetic
	result = ( parameter_matrix_element( N_shape, idx_G_1, idx_G_2, basis, supercell, P )*( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz ) )/( 2.0*sqrt( 3.0 ) );
	
	
	// Strain
	if( include_strain != 0 )
	{
	
		// Wave vector index
		int idx_basis = (basis.N_Gz_T/2 + idx_G_1[2] - idx_G_2[2]) + basis.N_Gz_T*(basis.N_Gy_T/2 + idx_G_1[1] - idx_G_2[1]) + basis.N_Gy_T*basis.N_Gz_T*(basis.N_Gx_T/2 + idx_G_1[0] - idx_G_2[0]);

		// Matrix term
		result += ( -0.5/sqrt( 3 ) )*parameters[0].P*( ( strain_G.xz[idx_basis] ) * ( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx ) + ( strain_G.yz[idx_basis] ) * ( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) + ( strain_G.zz[idx_basis] ) * ( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz ) ); 
		
		// Shape terms
		for ( idx_shape = 1; idx_shape < N_shape; idx_shape++ ) 
		{

			result += ( -0.5/sqrt( 3 ) )*( parameters[idx_shape].P - parameters[supercell[idx_shape].embedded].P )*( ( basis.Gx[idx_G_1_basis[0]] + basis.Gx[idx_G_2_basis[0]] + 2*basis.Kx )  * ( convolutions_char_strain[idx_shape].char_xz[idx_basis] ) + ( basis.Gy[idx_G_1_basis[1]] + basis.Gy[idx_G_2_basis[1]] + 2*basis.Ky ) * ( convolutions_char_strain[idx_shape].char_yz[idx_basis] ) + ( basis.Gz[idx_G_1_basis[2]] + basis.Gz[idx_G_2_basis[2]] + 2*basis.Kz ) * ( convolutions_char_strain[idx_shape].char_zz[idx_basis] ) );

		}

	}


	// Free dynamically allocated memory
	free( P );


	return result;

}
