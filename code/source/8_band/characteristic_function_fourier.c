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
//
//
// Shapes:
//
//		supercell[idx_shape].shape == 0   ==>  Cuboid
//		supercell[idx_shape].shape == 1   ==>  Sphere
//		supercell[idx_shape].shape == 2   ==>  Ellipsoid
//		supercell[idx_shape].shape == 3   ==>  Hemisphere
//		supercell[idx_shape].shape == 4   ==>  Hemiellipsoid (lens)
//		supercell[idx_shape].shape == 5   ==>  Rectangular-base pyramid
//		supercell[idx_shape].shape == 6   ==>  Rectangular-base truncated pyramid
//		supercell[idx_shape].shape == 7   ==>  N-sided pyramid
//		supercell[idx_shape].shape == 8   ==>  N-sided truncated pyramid
//		supercell[idx_shape].shape == 9   ==>  Cylinder
//		supercell[idx_shape].shape == 10  ==>  Cylindrical ring
//		supercell[idx_shape].shape == 11  ==>  Cone
//		supercell[idx_shape].shape == 12  ==>  Truncated cone
//		supercell[idx_shape].shape == 13  ==>  Truncated conic ring
//		supercell[idx_shape].shape == 14  ==>  Elliptical cone
//		supercell[idx_shape].shape == 15  ==>  Truncated elliptical cone
//		supercell[idx_shape].shape == 16  ==>  N-sided prism
//		supercell[idx_shape].shape == 17  ==>  Quantum well
//		supercell[idx_shape].shape == 18  ==>  Quantum wire
//
// The supercell is a cuboid centred at (0, 0, 0) with:
//
//	   	supercell[0].L1 = Length of supercell along x, i.e. -L1/2 <= x <= +L1/2
//      supercell[0].L2 = Length of supercell along y, i.e. -L2/2 <= x <= +L2/2
//      supercell[0].L3 = Length of supercell along z, i.e. -L3/2 <= x <= +L3/2
//
// Each "shape" placed in the supercell is centered at (Ox, Oy, Oz) with:
//
//      supercell[idx_shape].Ox = Centre position of "shape" along x
//      supercell[idx_shape].Oy = Centre position of "shape" along y
//      supercell[idx_shape].Oz = Centre position of "shape" along z

#include "8_band.h"

complex double characteristic_function_fourier( double Gx0, double Gy0, double Gz0, struct supercell_geometry *supercell, int idx_shape, int supercell_type )
{

	int supercell_size = N_shape*sizeof( _supercell_geometry ); // Memory size to later create new supercell struct, required for some shapes


	///////////////////////////////
	//// Set location of shape ////
	///////////////////////////////

	//	This prefactor specifies the location of the centre of the shape in the supercell, at position ( Ox, Oy, Oz )
	//
	//	The factor of 1/(supercell volume) is included in our definition of the Fourier transform, and hence in the Fourier transform of the QD characteristic function computed here

	complex double result;
	
	if ( supercell_type == 3 ) // Coulomb supercell
	{
	
		result = ( cexp( -1.0*I*( Gx0*supercell[idx_shape].Ox + Gy0*supercell[idx_shape].Oy + Gz0*supercell[idx_shape].Oz ) ) )/( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3 );
		
	}
	else if ( supercell_type == 0 || idx_shape == 0 ) // Electronic structure supercell
	{
	
		result = ( cexp( -1.0*I*( Gx0*supercell[idx_shape].Ox + Gy0*supercell[idx_shape].Oy + Gz0*supercell[idx_shape].Oz ) ) )/( supercell[0].L1*supercell[0].L2*supercell[0].L3 );
		
	}
	else // Strain supercell
	{
	
		result = ( cexp( -1.0*I*( Gx0*supercell[idx_shape].Ox + Gy0*supercell[idx_shape].Oy + Gz0*supercell[idx_shape].Oz ) ) )/( supercell[0].SL1*supercell[0].SL2*supercell[0].SL3 );
		
	}
	

	///////////////////////////////
	//// Set rotation of shape ////
	///////////////////////////////
	
	//	This prefactor specifies the rotation about the centre of the shape in the supercell, at position ( Ox, Oy, Oz )
	//
	//	The shapes rotate intrinsically about z, y', z'' with angles alpha, beta, gamma respectively (Proper Euler angles).
	
	double sin_alpha = sin( M_PI*supercell[idx_shape].alpha/180 );
	double cos_alpha = cos( M_PI*supercell[idx_shape].alpha/180 );
	double sin_beta = sin( M_PI*supercell[idx_shape].beta/180 );
	double cos_beta = cos( M_PI*supercell[idx_shape].beta/180 );
	double sin_gamma = sin( M_PI*supercell[idx_shape].gamma/180 );
	double cos_gamma = cos( M_PI*supercell[idx_shape].gamma/180 );
	
	double Gx = ( ( cos_beta * cos_alpha * cos_gamma ) - ( sin_alpha * sin_gamma ))*Gx0 + ( ( cos_gamma * cos_beta * sin_alpha ) + ( sin_gamma * cos_alpha ) )*Gy0 + ( -sin_beta * cos_gamma )*Gz0;
	
	double Gy = ( ( cos_beta * -cos_alpha * sin_gamma ) - ( sin_alpha * cos_gamma ) )*Gx0 + ( ( -sin_gamma * cos_beta * sin_alpha ) + ( cos_gamma * cos_alpha ) )*Gy0 + ( sin_beta * sin_gamma )*Gz0;
	
	double Gz = ( sin_beta * cos_alpha )*Gx0 + ( sin_alpha * sin_beta )*Gy0 + ( cos_beta )*Gz0;
	

	double G_mag = sqrt( Gx*Gx + Gy*Gy + Gz*Gz ); // Magnitude of wave vector
	double G_par = sqrt( Gx*Gx + Gy*Gy );         // Magnitude of wave vector in growth plane


	/////////////////////////////
	//// Interface softening ////
	/////////////////////////////

	// Soften interfaces by convolving the real space characteristic function with a real space Gaussian distribution with "softening length" defined by the standard deviation
	//
	// Fourier transform of a Gaussian is a Gaussian, so by convolution theorem simply multiply Fourier transform of characteristic function by a Gaussian

	if( softening != 0 )
	{

		result *= exp( -0.5*pow( sigma_softening, 2.0 )*( pow( Gx, 2.0 ) + pow( Gy, 2.0 ) + pow( Gz, 2.0 ) ) );

	}



	////////////////
	//// Cuboid ////
	////////////////

	if( supercell[idx_shape].shape == 0 )
	{

		// Note:
		//
		//	The body centre of the cuboid is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The cuboid then extends between: supercell[idx_shape].Ox - supercell[idx_shape].L1/2 <= x <= supercell[idx_shape].Ox + supercell[idx_shape].L1/2
		//							   supercell[idx_shape].Oy - supercell[idx_shape].L2/2 <= y <= supercell[idx_shape].Oy + supercell[idx_shape].L2/2
		//							   supercell[idx_shape].Oz - supercell[idx_shape].L3/2 <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L3/2
		//
		//	supercell[idx_shape].L1 = Length of cuboid along x
		//	supercell[idx_shape].L2 = Length of cuboid along y
		//	supercell[idx_shape].L3 = Length of cuboid along z
		//	supercell[idx_shape].L4 = Not used
		//	supercell[idx_shape].L5 = Not used
		//
		//	GSL defines sinc(x) = sin(pi*x)/(pi*x), hence the factor of pi in the denominator, which does not appear in the analytical expressions

		result *=  supercell[idx_shape].L1*gsl_sf_sinc( ( 0.5*Gx*supercell[idx_shape].L1 )/M_PI )
				*supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gy*supercell[idx_shape].L2 )/M_PI )
				*supercell[idx_shape].L3*gsl_sf_sinc( ( 0.5*Gz*supercell[idx_shape].L3 )/M_PI );
			

	}


	////////////////
	//// Sphere ////
	////////////////

	if( supercell[idx_shape].shape == 1 )
	{

		// Note:
		//
		//	The centre of the sphere is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		// 	The sphere then extends out to a radius r = supercell[idx_shape].L1 about this centre point
		//
		//	supercell[idx_shape].L1 = Radius of sphere
		//	supercell[idx_shape].L2 = Not used
		//	supercell[idx_shape].L3 = Not used
		//	supercell[idx_shape].L4 = Not used
		//	supercell[idx_shape].L5 = Not used

		if( fabs( G_mag ) < EPS ) // |G| = 0
		{

			result *= ( 4.0*M_PI*pow( supercell[idx_shape].L1, 3.0 ) )/3.0 + I*0.0; // Volume

		}
		else
		{

			result *= ( 4.0*M_PI*( sin( G_mag*supercell[idx_shape].L1 ) - G_mag*supercell[idx_shape].L1*cos( G_mag*supercell[idx_shape].L1 ) ) )/( pow( G_mag, 3.0 ) );

		}

	}


	///////////////////
	//// Ellipsoid ////
	///////////////////

	if( supercell[idx_shape].shape == 2 )
	{

		// Note:
		//
		//	The centre of the ellipsoid is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		// 	The ellipsoid then extends between = -supercell[idx_shape].L1 <= x <= supercell[idx_shape].L1, -supercell[idx_shape].L2 <= y <= supercell[idx_shape].L2, -supercell[idx_shape].L3 <= z <= supercell[idx_shape].L3
		//
		//	supercell[idx_shape].L1 = x semi-axis
		//	supercell[idx_shape].L2 = y semi-axis
		//	supercell[idx_shape].L3 = z semi-axis
		//	supercell[idx_shape].L4 = Not used
		//	supercell[idx_shape].L5 = Not used

		G_mag = sqrt( pow( Gx*supercell[idx_shape].L1, 2.0 ) + pow( Gy*supercell[idx_shape].L2, 2.0 ) + pow( Gz*supercell[idx_shape].L3, 2.0 ) ); // Magnitude of re-scaled wave vector

		if( fabs( G_mag ) < EPS ) // |G| = 0
		{

			result *= ( 4.0*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3 )/3.0 + I*0.0; // Volume

		}
		else
		{

			result *= ( 4.0*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3*( sin( G_mag ) - G_mag*cos( G_mag ) ) )/( pow( G_mag, 3.0 ) );

		}

	}

	////////////////////
	//// Hemisphere ////
	////////////////////

	if( supercell[idx_shape].shape == 3 )
	{

		// Note:
		//
		//	The centre of the base of the hemisphere is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The hemisphere then extends out to a radius r = supercell[idx_shape].L1 upward in z about this centre point
		//
		//	supercell[idx_shape].L1 = Radius of hemisphere
		//	supercell[idx_shape].L2 = Not used
		//	supercell[idx_shape].L3 = Not used
		//	supercell[idx_shape].L4 = Not used
		//	supercell[idx_shape].L5 = Not used

		if( fabs( G_mag ) < EPS ) // |G| = 0
		{

			result *= ( 2.0*M_PI*pow( supercell[idx_shape].L1, 3.0 ) )/3.0; // Volume

		}
		else if( ( fabs( G_par ) < EPS ) && ( fabs( Gz ) > EPS ) ) // G_par = 0, Gz != 0
		{

			result *= (   2.0*M_PI*( sin( G_mag*supercell[idx_shape].L1 ) - G_mag*supercell[idx_shape].L1*cos( G_mag*supercell[idx_shape].L1 ) ) )/( pow( G_mag, 3.0 ) )
				   - ( ( 2.0*M_PI*I*pow(       supercell[idx_shape].L1, 2.0 ) )/Gz )*( 0.5 - ( supercell[idx_shape].L1*Gz*sin( supercell[idx_shape].L1*Gz ) + cos( supercell[idx_shape].L1*Gz ) - 1.0 )/( pow( Gz*supercell[idx_shape].L1, 2.0 ) ) );

		}
		else if( ( fabs( Gz ) < EPS ) && ( fabs( G_par ) > EPS ) ) // Gz = 0, G_par != 0
		{

			result *= ( 2.0*M_PI*( sin( G_mag*supercell[idx_shape].L1 ) - G_mag*supercell[idx_shape].L1*cos( G_mag*supercell[idx_shape].L1 ) ) )/( pow( G_mag, 3.0 ) );

		}
		else
		{

			result *= ( 2.0*M_PI  *( sin( G_mag*supercell[idx_shape].L1 ) - G_mag*supercell[idx_shape].L1*       cos( G_mag*supercell[idx_shape].L1 ) ) )/( pow( G_mag, 3.0 ) )
				   - ( 2.0*M_PI*I*             supercell[idx_shape].L1                            *gsl_sf_bessel_J1( G_par*supercell[idx_shape].L1 )   )/(      G_par*Gz     )
				   + ( 2.0*M_PI*I*        pow( supercell[idx_shape].L1, 2.0 ) *I_J0_GSL( Gz*supercell[idx_shape].L1, G_par*supercell[idx_shape].L1 )   )/(            Gz     );
		}

	}


	//////////////////////////////
	//// Hemiellipsoid (lens) ////
	//////////////////////////////

	if( supercell[idx_shape].shape == 4 )
	{

		// Note:
		//
		//	The centre of the base of the lens is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		// 	The lens then extends between = -supercell[idx_shape].L1 <= x <= supercell[idx_shape].L1, -supercell[idx_shape].L2 <= y <= supercell[idx_shape].L2, 0 <= z <= supercell[idx_shape].L3
		//
		//	supercell[idx_shape].L1 = x semi-axis
		//	supercell[idx_shape].L2 = y semi-axis
		//	supercell[idx_shape].L3 = z semi-axis
		//	supercell[idx_shape].L4 = Not used
		//	supercell[idx_shape].L5 = Not used

		G_mag = sqrt( pow( Gx*supercell[idx_shape].L1, 2.0 ) + pow( Gy*supercell[idx_shape].L2, 2.0 ) + pow( Gz*supercell[idx_shape].L3, 2.0 ) ); // Magnitude of re-scaled          wave vector
		G_par = sqrt( pow( Gx*supercell[idx_shape].L1, 2.0 ) + pow( Gy*supercell[idx_shape].L2, 2.0 )                                          ); // Magnitude of re-scaled in-plane wave vector

		if( fabs( G_mag ) < EPS ) // |G| = 0
		{

			result *= ( 2.0*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3 )/3.0 + I*0.0; // Volume

		}
		else if( fabs( G_par ) < EPS )  // G_par = 0, Gz != 0
		{

			result *= ( 2.0*M_PI  *supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3*( sin( G_mag ) - G_mag*cos( G_mag ) ) )/( pow( G_mag, 3.0 ) )
				   - ( 2.0*M_PI*I*supercell[idx_shape].L1*supercell[idx_shape].L2/Gz )*( 0.5 - ( cos( Gz*supercell[idx_shape].L3 ) + Gz*supercell[idx_shape].L3*sin( supercell[idx_shape].L3*Gz ) - 1.0 )/( pow( Gz*supercell[idx_shape].L3, 2.0 ) ) );

		}
		else if( fabs( Gz ) < EPS ) //  G_par != 0, Gz = 0
		{

			result *= ( 2.0*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3*( sin( G_mag ) - G_mag*cos( G_mag ) ) )/( pow( G_mag, 3.0 ) );

		}
		else
		{

			result *= ( 2.0*M_PI  *supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3*( sin( G_mag ) - G_mag*cos( G_mag ) ) )/( pow( G_mag, 3.0 ) )
				   - ( 2.0*M_PI*I*supercell[idx_shape].L1*supercell[idx_shape].L2/Gz )*( gsl_sf_bessel_J1( G_par )/G_par - I_J0_GSL( Gz*supercell[idx_shape].L3, G_par ) );

		}

	}


	/////////////////
	//// Pyramid ////
	/////////////////

	if( supercell[idx_shape].shape == 5 )
	{

		// Note:
		//
		//	The centre of the base of the pyramid is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The pyramid then extends between: supercell[idx_shape].Ox - supercell[idx_shape].L1 <= x <= supercell[idx_shape].Ox + supercell[idx_shape].L1
		//							    supercell[idx_shape].Oy - supercell[idx_shape].L2 <= y <= supercell[idx_shape].Oy + supercell[idx_shape].L2
		//							    supercell[idx_shape].Oz                            <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L3
		//
		//	supercell[idx_shape].L1 = Half the length of pyramid base along x
		//  supercell[idx_shape].L2 = Half the length of pyramid base along y
		//  supercell[idx_shape].L3 = Height          of pyramid      along z
		//  supercell[idx_shape].L4 = Not used
		//  supercell[idx_shape].L5 = Not used

		int N_side = 4;

		complex double seg_result = 0 + 0*I;

		for( int idx_side = 0; idx_side < N_side; idx_side++ ) // Considers pyramid as split up into 4 segments
		{

			seg_result += characteristic_function_fourier_pyramid_segment( Gx, Gy, Gz, idx_side, 0, supercell, idx_shape );

		}

		result *= seg_result;
	}


	///////////////////////////
	//// Truncated pyramid ////
	///////////////////////////

	if( supercell[idx_shape].shape == 6 )
	{

		// Note:
		//
		//	The centre of the base of the pyramid is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The truncated pyramid base extends between: supercell[idx_shape].Ox - supercell[idx_shape].L1 <= x <= supercell[idx_shape].Ox + supercell[idx_shape].L1 & supercell[idx_shape].Oy - supercell[idx_shape].L2 <= y <= supercell[idx_shape].Oy + supercell[idx_shape].L2
		//	The truncated pyramid top extends between: supercell[idx_shape].Ox - supercell[idx_shape].L3 <= x <= supercell[idx_shape].Ox + supercell[idx_shape].L3
		//	The truncated pyramid height extends between: supercell[idx_shape].Oz <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L4
		//
		//	supercell[idx_shape].L1 = Half the length of truncated pyramid base along x
		//  supercell[idx_shape].L2 = Half the length of truncated pyramid base along y
		//	supercell[idx_shape].L3 = Half the length of truncated pyramid top  along x
		//  supercell[idx_shape].L4 = Height          of truncated pyramid      along z
		//  supercell[idx_shape].L5 = Not used

		//  Considers truncated pyramid as a larger pyramid with its top cut off

		int N_side = 4;

		complex double seg_result = 0 + 0*I;

		struct supercell_geometry *supercell_copy = ( _supercell_geometry* ) malloc( supercell_size );

		memcpy( supercell_copy, supercell, supercell_size );

		supercell_copy[idx_shape].L3 = ( supercell[idx_shape].L4*supercell[idx_shape].L1 )/( supercell[idx_shape].L1 - supercell[idx_shape].L3 );

		for( int idx_side = 0; idx_side < N_side; idx_side++ )
		{

			seg_result += characteristic_function_fourier_pyramid_segment( Gx, Gy, Gz, idx_side, 0, supercell_copy, idx_shape );

		}

		supercell_copy[idx_shape].L1 =   supercell[idx_shape].L3;
		supercell_copy[idx_shape].L2 = ( supercell[idx_shape].L2*supercell[idx_shape].L3 )/  supercell[idx_shape].L1;
		supercell_copy[idx_shape].L3 = ( supercell[idx_shape].L4*supercell[idx_shape].L3 )/( supercell[idx_shape].L1 - supercell[idx_shape].L3 );

		for( int idx_side = 0; idx_side < N_side; idx_side++)
		{

			seg_result -= cexp( -1.0*I*Gz*supercell[idx_shape].L4 )*characteristic_function_fourier_pyramid_segment( Gx, Gy, Gz, idx_side, 0, supercell_copy, idx_shape );

		}

		result *= seg_result;

		free( supercell_copy );

	}


	/////////////////////////
	//// N-sided pyramid ////
	/////////////////////////

	if( supercell[idx_shape].shape == 7 )
	{

		// Note:
		//
		//	The centre of the base of the pyramid is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The pyramid then extends between: supercell[idx_shape].Ox - supercell[idx_shape].L1/2 <= x <= supercell[idx_shape].Ox + supercell[idx_shape].L1/2
		//							    supercell[idx_shape].Oy - supercell[idx_shape].L2/2 <= y <= supercell[idx_shape].Oy + supercell[idx_shape].L2/2
		//							    supercell[idx_shape].Oz                             <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L3
		//
		//	supercell[idx_shape].L1 = Distance from centre of pyramid base to side
		//  supercell[idx_shape].L2 = Height of pyramid along z
		//  supercell[idx_shape].L3 = Number of sides
		//  supercell[idx_shape].L4 = Not used
		//  supercell[idx_shape].L5 = Not used

		int N_side = lround( 0.1*supercell[idx_shape].L3 );

		complex double seg_result = 0 + 0*I;

		for( int idx_side = 0; idx_side < N_side; idx_side++ ) // Considers pyramid as split up into N segments
		{

			seg_result += characteristic_function_fourier_pyramid_segment( Gx, Gy, Gz, idx_side, N_side, supercell, idx_shape );

		}

		result *= seg_result;

	}

	///////////////////////////////////
	//// N-sided truncated pyramid ////
	///////////////////////////////////

	if( supercell[idx_shape].shape == 8 )
	{

		// Note:
		//
		//	The centre of the base of the pyramid is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The truncated pyramid base   extends between: supercell[idx_shape].Ox - supercell[idx_shape].L1/2 <= x <= supercell[idx_shape].Ox + supercell[idx_shape].L1/2 & supercell[idx_shape].Oy - supercell[idx_shape].L2/2 <= y <= supercell[idx_shape].Oy + supercell[idx_shape].L2/2
		//	The truncated pyramid top    extends between: supercell[idx_shape].Ox - supercell[idx_shape].L3/2 <= x <= supercell[idx_shape].Ox + supercell[idx_shape].L3/2
		//	The truncated pyramid height extends between: supercell[idx_shape].Oz <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L4
		//
		//	supercell[idx_shape].L1 = Distance from centre of truncated pyramid base to side
		//  supercell[idx_shape].L2 = Distance from centre of truncated pyramid top to side
		//	supercell[idx_shape].L3 = Height of truncated pyramid along z
		//  supercell[idx_shape].L4 = Number of sides
		//  supercell[idx_shape].L5 = Not used

		//  Considers truncated pyramid as a larger pyramid with its top cut off

		int N_side = lround( 0.1*supercell[idx_shape].L4 );

		complex double seg_result = 0 + 0*I;

		struct supercell_geometry *supercell_copy = ( _supercell_geometry* ) malloc( supercell_size );

		memcpy( supercell_copy, supercell, supercell_size );

		supercell_copy[idx_shape].L2 = ( supercell[idx_shape].L1*supercell[idx_shape].L3 )/( supercell[idx_shape].L1 - supercell[idx_shape].L2 );

		for( int idx_side = 0; idx_side < N_side; idx_side++ )
		{

			seg_result += characteristic_function_fourier_pyramid_segment( Gx, Gy, Gz, idx_side, N_side, supercell_copy, idx_shape );

		}


		supercell_copy[idx_shape].L1 =   supercell[idx_shape].L2;
		supercell_copy[idx_shape].L2 = ( supercell[idx_shape].L3*supercell[idx_shape].L2 )/( supercell[idx_shape].L1 - supercell[idx_shape].L2 );

		for( int idx_side = 0; idx_side < N_side; idx_side++ )
		{

			seg_result -= cexp( -1.0*I*Gz*supercell[idx_shape].L3 )*characteristic_function_fourier_pyramid_segment( Gx, Gy, Gz, idx_side, N_side, supercell_copy, idx_shape );

		}

		result *= seg_result;

		free( supercell_copy );

	}


	//////////////////
	//// Cylinder ////
	//////////////////

	if( supercell[idx_shape].shape == 9 )
	{

		// Note:
		//
		//	The centre of the cylinder is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The cylinder then extends our to a radius r = supercell[idx_shape].L1 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy ) in the x-y plane,
		//	and between supercell[idx_shape].Oz - supercell[idx_shape].L2/2 <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L2/2 along z
		//
		//	supercell[idx_shape].L1 = Radius of cylinder in the x-y plane
		//	supercell[idx_shape].L2 = Height of cylinder along z
		//	supercell[idx_shape].L3 = Not used
		//	supercell[idx_shape].L4 = Not used
		//	supercell[idx_shape].L5 = Not used

		if( fabs( G_mag ) < EPS ) // |G| = 0
		{

			result *= M_PI*pow( supercell[idx_shape].L1, 2.0 )*supercell[idx_shape].L2; // Volume

		}
		else if( ( fabs( G_par ) < EPS ) && ( fabs( Gz ) > EPS ) ) // G_par = 0, Gz != 0
		{

			result *= 2.0*M_PI*pow( supercell[idx_shape].L1, 2.0 )*sin( 0.5*Gz*supercell[idx_shape].L2 ) / Gz;

		}
		else if( ( fabs( Gz ) < EPS ) && ( fabs( G_par ) > EPS ) ) // Gz = 0, G_par != 0
		{

			result *= ( 2.0*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*gsl_sf_bessel_J1( G_par*supercell[idx_shape].L1 ) )/G_par;

		}
		else
		{

			result *= ( 4.0*M_PI*supercell[idx_shape].L1*sin( 0.5*Gz*supercell[idx_shape].L2 )*gsl_sf_bessel_J1( G_par*supercell[idx_shape].L1 ) )/( G_par*Gz );

		}

	}


	//////////////////////////
	//// Cylindrical ring ////
	//////////////////////////

	if( supercell[idx_shape].shape == 10 )
	{

		// Note:
		//
		//	The centre of the cylindtrical ring is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The outside cylinder then extends our to a radius r = supercell[idx_shape].L1 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy ) in the x-y plane.
		//	The inside cylinder then extends our to a radius r = supercell[idx_shape].L2 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy ) in the x-y plane,
		//	The height extends between supercell[idx_shape].Oz - supercell[idx_shape].L3/2 <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L3/2 along z
		//
		//	supercell[idx_shape].L1 = Outside radius of cylinder in the x-y plane
		//	supercell[idx_shape].L2 = Inside radius of cylinder in the x-y plane
		//	supercell[idx_shape].L3 = Height of ring along z
		//	supercell[idx_shape].L4 = Not used
		//	supercell[idx_shape].L5 = Not used

		// Considers ring as cylinder with inner cylinder removed

		struct supercell_geometry *supercell_copy = ( _supercell_geometry* ) malloc( supercell_size );

		memcpy( supercell_copy, supercell, supercell_size );

		supercell_copy[idx_shape].shape = 9;
		supercell_copy[idx_shape].L2    = supercell[idx_shape].L3;

		result = characteristic_function_fourier( Gx0, Gy0, Gz0, supercell_copy, idx_shape, supercell_type );


		supercell_copy[idx_shape].L1 = supercell[idx_shape].L2;

		result -= characteristic_function_fourier( Gx0, Gy0, Gz0, supercell_copy, idx_shape, supercell_type );

		free( supercell_copy );

	}


	//////////////
	//// Cone ////
	//////////////

	if( supercell[idx_shape].shape == 11 )
	{

		// Note:
		//
		//	The centre of the base of the cone is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The cone extends our to a radius r = supercell[idx_shape].L1 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy )
		//  in the x-y plane at its base and between supercell[idx_shape].Oz <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L2
		//
		//	supercell[idx_shape].L1 = Radius of cone
		//  supercell[idx_shape].L2 = Height of cone
		//  supercell[idx_shape].L3 = Not used
		//  supercell[idx_shape].L4 = Not used
		//  supercell[idx_shape].L5 = Not used

		if( fabs( G_mag ) < EPS ) // |G| = 0
		{

			result *= ( M_PI*pow( supercell[idx_shape].L1, 2.0 )*supercell[idx_shape].L2 )/3.0 + I*0.0; // Volume

		}
		else if( fabs( G_par ) < EPS ) // G_par = 0, Gz != 0
		{

			result *= ( 2.0*I*M_PI*pow( supercell[idx_shape].L1, 2.0 )*( ( 1.0 - I*Gz*supercell[idx_shape].L2 - cexp( -1.0*I*Gz*supercell[idx_shape].L2 ) )/pow( Gz*supercell[idx_shape].L2, 2.0 ) - 0.5 ) )/Gz;

		}
		else if( fabs( Gz ) < EPS ) // G_par != 0, Gz = 0
		{

			result *= ( 2.0*M_PI  *pow( supercell[idx_shape].L1, 2.0 )*supercell[idx_shape].L2 )*( gsl_sf_bessel_J1( G_par*supercell[idx_shape].L1 )/( G_par*supercell[idx_shape].L1 ) - I_2_Gz_GSL( G_par*supercell[idx_shape].L1 ) );

		}
		else
		{

			result *= ( 2.0*M_PI  *supercell[idx_shape].L1*supercell[idx_shape].L1*I*cexp( -1.0*I*Gz*supercell[idx_shape].L2 )*I_2_GSL( G_par*supercell[idx_shape].L1, Gz*supercell[idx_shape].L2 ) )/(       Gz )
				   - ( 2.0*M_PI*I*supercell[idx_shape].L1                                                       *gsl_sf_bessel_J1( G_par*supercell[idx_shape].L1 )                             )/( G_par*Gz );

		}

	}


	////////////////////////
	//// Truncated cone ////
	////////////////////////

	if( supercell[idx_shape].shape == 12 )
	{

		// Note:
		//
		//	The centre of the base of the cone is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The cone extends our to a radius r = supercell[idx_shape].L1 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy) in the x-y plane at its base.
		//	The cone extends our to a radius r = supercell[idx_shape].L2 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].0z + supercell[idx_shape].L3) in the x-y plane at its truncated height. 			//	The cone extends between supercell[idx_shape].Oz <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L3
		//
		//	supercell[idx_shape].L1 = Outer radius at 0 height
		//  supercell[idx_shape].L2 = Inner radius at maximum height
		//  supercell[idx_shape].L3 = Truncated height
		//  supercell[idx_shape].L4 = Not used
		//  supercell[idx_shape].L5 = Not used

		if( fabs( G_mag ) < EPS ) // |G| = 0
		{

			result *= ( M_PI*supercell[idx_shape].L3*( supercell[idx_shape].L1*supercell[idx_shape].L1 + supercell[idx_shape].L1*supercell[idx_shape].L2 + supercell[idx_shape].L2*supercell[idx_shape].L2 ) )/3.0; // Volume

		}
		else if( fabs( G_par ) < EPS ) // G_par = 0, Gz != 0
		{

			result *= ( ( 2.0*I*M_PI )/Gz )*( ( 0.5*  pow( supercell[idx_shape].L2, 2.0 )*( cexp( -1.0*I*Gz*supercell[idx_shape].L3 ) - 1.0 ) )
							   		  + (       pow( supercell[idx_shape].L1 - supercell[idx_shape].L2, 2.0 )/ pow( Gz*supercell[idx_shape].L3, 2.0 )*( cexp( -1.0*I*Gz*supercell[idx_shape].L3 )*( ( I*Gz*supercell[idx_shape].L2*supercell[idx_shape].L3 )/( supercell[idx_shape].L1 - supercell[idx_shape].L2 ) - 1.0 ) + ( 1.0 - ( I*supercell[idx_shape].L1*Gz*supercell[idx_shape].L3 )/( supercell[idx_shape].L1 - supercell[idx_shape].L2 ) ) ) )
							   		  - ( 0.5*( pow( supercell[idx_shape].L1, 2.0 ) - pow( supercell[idx_shape].L2, 2.0 ) ) ) );

		}
		else if( fabs( Gz ) < EPS ) // G_par != 0, Gz = 0
		{

			result *= ( 2.0*M_PI*supercell[idx_shape].L3 )*( (      supercell[idx_shape].L2                                 *gsl_sf_bessel_J1( G_par *supercell[idx_shape].L2 )/G_par )
												  - ( pow( supercell[idx_shape].L1 - supercell[idx_shape].L2, 2.0 )*       I_2_TC_Gz_GSL( G_par, supercell[idx_shape].L2, supercell[idx_shape].L1 ) )
												  + (    ( supercell[idx_shape].L1/( supercell[idx_shape].L1 - supercell[idx_shape].L2 ) )*( supercell[idx_shape].L1*gsl_sf_bessel_J1( G_par*supercell[idx_shape].L1 ) - supercell[idx_shape].L2*gsl_sf_bessel_J1( G_par*supercell[idx_shape].L2 ) )/G_par ) );

		}
		else
		{

			result *= ( 2.0*I*M_PI/Gz )*( (   supercell[idx_shape].L2*( cexp( -1.0*I*Gz*supercell[idx_shape].L3 ) - 1.0 )*gsl_sf_bessel_J1( G_par*supercell[idx_shape].L2 )/G_par )
								   + ( pow( supercell[idx_shape].L1 - supercell[idx_shape].L2, 2.0 )*cexp( ( -1.0*I*Gz*supercell[idx_shape].L3*supercell[idx_shape].L1 )/( supercell[idx_shape].L1 - supercell[idx_shape].L2 ) )*I_2_TC_GSL( G_par, supercell[idx_shape].L2, supercell[idx_shape].L1, Gz*supercell[idx_shape].L3 ) )
								   - ( ( supercell[idx_shape].L1*gsl_sf_bessel_J1( G_par*supercell[idx_shape].L1 ) - supercell[idx_shape].L2*gsl_sf_bessel_J1( G_par*supercell[idx_shape].L2 ) )/G_par ) );

		}

	}


	//////////////////////////////
	//// Truncated conic ring ////
	//////////////////////////////

	if( supercell[idx_shape].shape == 13 )
	{

		// Note:
		//
		//	The centre of the base of the conic ring is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The outside truncated cone extends our to a radius r = supercell[idx_shape].L1 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy) in the x-y plane at its base.
		//	The outside truncated cone extends our to a radius r = supercell[idx_shape].L2 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy in the x-y plane at its top.
		//	The inside cylinder extends to a radius r = supercell[idx_shape].L3 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy in the x-y plane.
		//	The cone extends between supercell[idx_shape].Oz <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L4
		//
		//	supercell[idx_shape].L1 = Outer radius at 0 height
		//  supercell[idx_shape].L2 = Inner radius at maximum height
		//  supercell[idx_shape].L3 = Inside cylinder radius
		//  supercell[idx_shape].L4 = Truncated height
		//  supercell[idx_shape].L5 = Not used

		//  Considers ring as truncated cone with inner cylinder removed. Inner cylinder radius must be less than truncated cone top radius

		struct supercell_geometry *supercell_copy = ( _supercell_geometry* ) malloc( supercell_size );

		memcpy( supercell_copy, supercell, supercell_size );

		supercell_copy[idx_shape].shape = 12;
		supercell_copy[idx_shape].L3    = supercell[idx_shape].L4;

		result = characteristic_function_fourier( Gx0, Gy0, Gz0, supercell_copy, idx_shape, supercell_type );

		supercell_copy[idx_shape].shape = 9;
		supercell_copy[idx_shape].L1    = supercell[idx_shape].L3;
		supercell_copy[idx_shape].L2    = supercell[idx_shape].L4;

		result -= cexp( -0.5*I*Gz*supercell[idx_shape].L4 )*characteristic_function_fourier( Gx0, Gy0, Gz0, supercell_copy, idx_shape, supercell_type );

		free( supercell_copy );

	}


	/////////////////////////
	//// Elliptical cone ////
	/////////////////////////

	if( supercell[idx_shape].shape == 14 )
	{

		// Note:
		//
		//	The centre of the base of the cone is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The cone extends our to a distance r = supercell[idx_shape].L1 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy ) along the semi-major axis.
		//	The cone extends our to a distance r = supercell[idx_shape].L2 about ( supercell[idx_shape].Ox, supercell[idx_shape].Oy ) along the semi-minor axis.
		//  The cone extends between supercell[idx_shape].Oz <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L3
		//
		//	supercell[idx_shape].L1 = Semi-major axis of elliptical cone
		//  supercell[idx_shape].L2 = Semi-minor axis of elliptical cone
		//  supercell[idx_shape].L3 = Heigth of cone
		//  supercell[idx_shape].L4 = Not used
		//  supercell[idx_shape].L5 = Not used

		G_par = sqrt( pow( Gx*supercell[idx_shape].L1, 2.0 ) + pow( Gy*supercell[idx_shape].L2, 2.0 ) ); // Magnitude of re-scaled wave vector in [001] plane

		if( fabs( G_mag ) < EPS ) // |G| = 0
		{

			result *= ( M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3 )/3.0 + I*0.0; // Volume

		}

		else if( fabs( G_par ) < EPS ) // G_par = 0, Gz != 0
		{

			result *= ( 1.0/Gz )*( 2.0*I*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*( ( 1.0 - I*Gz*supercell[idx_shape].L3 - cexp( -1.0*I*Gz*supercell[idx_shape].L3 ) )/pow( Gz*supercell[idx_shape].L3, 2.0 ) - 0.5 ) );

		}

		else if( fabs( Gz ) < EPS ) // G_par != 0, Gz = 0
		{

			result *= ( 2.0*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3 )*( gsl_sf_bessel_J1( G_par )/G_par - I_2_Gz_GSL( G_par ) );

		}


		else
		{

			result *= ( 2.0*M_PI  *supercell[idx_shape].L1*supercell[idx_shape].L2*I*cexp( -1.0*I*Gz*supercell[idx_shape].L3 )*I_2_GSL( G_par, Gz*supercell[idx_shape].L3 ) )/(       Gz )
				   - ( 2.0*M_PI*I*supercell[idx_shape].L1*supercell[idx_shape].L2                               *gsl_sf_bessel_J1( G_par )                             )/( G_par*Gz );

		}

	}


	///////////////////////////////////
	//// Truncated elliptical cone ////
	///////////////////////////////////

	if( supercell[idx_shape].shape == 15 )
	{

		// Note:
		//
		//	The centre of the base of the cone is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The cone extends out to semi-major axes a = supercell[idx_shape].L1 along ( supercell[idx_shape].Ox) at its base.
		//	The cone extends out to semi-minor axes b = supercell[idx_shape].L2 along ( supercell[idx_shape].Oy) at its base.
		//	The cone extends out to semi-major axis f*a = ( supercell[idx_shape].L3 / supercell[idx_shape].L1 )*supercell[idx_shape].L1 along ( supercell[idx_shape].Ox) at (supercell[idx_shape].0z + supercell[idx_shape].L4) at its truncated height.
		//	The cone extends out to semi-minor axis f*b = ( supercell[idx_shape].L3 / supercell[idx_shape].L1 )*supercell[idx_shape].L2 along ( supercell[idx_shape].Oy) at (supercell[idx_shape].0z + supercell[idx_shape].L4) at its truncated height.
		//	The cone extends between supercell[idx_shape].Oz <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L4
		//
		//	supercell[idx_shape].L1 = Semi-major axis at 0 height
		//  supercell[idx_shape].L2 = Semi-minor axis at 0 height
		//  supercell[idx_shape].L3 = Semi-major axis at truncated height
		//  supercell[idx_shape].L4 = Truncated height
		//  supercell[idx_shape].L5 = Not used

		G_par = sqrt( pow( Gx*supercell[idx_shape].L1, 2.0 ) + pow( Gy*supercell[idx_shape].L2, 2.0 ) ); // Magnitude of re-scaled wave vector in [001] plane

		double f = supercell[idx_shape].L3/supercell[idx_shape].L1;

		if( fabs( G_mag ) < EPS ) // |G| = 0
		{

			result *= ( M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L4*( 1.0 + f + pow( f, 2.0 ) ) )/3.0; // Volume

		}

		else if( fabs( G_par ) < EPS ) // G_par = 0, Gz != 0
		{

			result *= ( 2.0*M_PI*I*supercell[idx_shape].L1*supercell[idx_shape].L2*( ( 0.5*pow( f, 2.0 )*( cexp( -1.0*I*Gz*supercell[idx_shape].L4 ) - 1.0 ) )
																	 + ( ( f - 1.0 )/( Gz*Gz*supercell[idx_shape].L4*supercell[idx_shape].L4 )*( cexp( -1.0*I*Gz*supercell[idx_shape].L4 )*( -1.0*I*Gz*supercell[idx_shape].L4*f - f + 1.0 ) + ( I*Gz*supercell[idx_shape].L4 + f - 1.0 ) ) )
																	 - ( 0.5*( 1.0 - pow( f, 2.0 ) ) ) ) / Gz );

		}

		else if( fabs( Gz ) < EPS ) // G_par != 0, Gz = 0
		{

			result *= ( 2.0*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L4 )*( ( f*gsl_sf_bessel_J1( f*G_par )/G_par ) - ( 1.0*I_2_ETC_Gz_GSL( G_par, f ) )/( 1.0 - f ) + ( ( gsl_sf_bessel_J1( G_par ) - f*gsl_sf_bessel_J1( f*G_par ) )/( ( 1.0 - f )*G_par ) ) );

		}


		else
		{

			result *= ( ( 2.0*I*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2 )/Gz )*( ( f*( cexp( -1.0*I*Gz*supercell[idx_shape].L4 ) - 1.0 )*gsl_sf_bessel_J1( f*G_par )/G_par ) + ( cexp( -1.0*I*Gz*supercell[idx_shape].L4/( 1.0 - f ) )*I_2_ETC_GSL( G_par, Gz*supercell[idx_shape].L4, f ) ) - ( ( gsl_sf_bessel_J1( G_par ) - f*gsl_sf_bessel_J1( f*G_par ) ) / G_par ) );

		}

	}


	///////////////////////
	//// N-sided prism ////
	///////////////////////

	if( supercell[idx_shape].shape == 16 )
	{

		// Note:
		//
		//	The centre of the prism is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//	The prism then extends vertically between:     supercell[idx_shape].Oz - supercell[idx_shape].L2 / 2  <= z <= supercell[idx_shape].Oz + supercell[idx_shape].L2 / 2
		//
		//  Without rotation, the first side of the prism base is defined to be parallel to the y-axis with equal extent either side of the x-axis.
		//
		//	supercell[idx_shape].L1 = Distance from centre of prism base to side (Side length equal to 2*L1*tan( pi/N ), where N is the number of sides )
		//  supercell[idx_shape].L2 = Height of prism along z
		//  supercell[idx_shape].L3 = Number of sides
		//  supercell[idx_shape].L4 = Not used
		//  supercell[idx_shape].L5 = Not used

		int N_side = lround( 0.1*supercell[idx_shape].L3 );

		complex double seg_result = 0 + 0*I;

		for( int idx_side = 0; idx_side < N_side; idx_side++ ) // Considers shape as consisting of N segments
		{

			seg_result += characteristic_function_fourier_prism_segment( Gx, Gy, Gz, idx_side, N_side, supercell, idx_shape );

		}

		result *= seg_result;

	}


	//////////////////////
	//// Quantum well ////
	//////////////////////

	// Insert a quantum well with quantum confinement along either the [100], [010], or [001] directions. Differs from a cuboid that extends in-plane across the supercell in that the quantum well shape extends to the larger strain supercell boundaries for the relevant dimensions when calculating the strain.

	if( supercell[idx_shape].shape == 17 )
	{

		// Note:
		//
		//	The body centre of the quantum well is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//  The quantum confinement direction is given by supercell[idx_shape].L1
		//
		//	The quantum well then extends along the quantum confinement direction between: supercell[idx_shape].Ox/y/z - supercell[idx_shape].L2/2 <= x/y/z <= supercell[idx_shape].Ox/y/z + supercell[idx_shape].L2/2
		//
		//	supercell[idx_shape].L1 = Quantum confinement direction (0: [100],  1: [010],  2: [001])
		//	supercell[idx_shape].L2 = Length of quantum well
		//	supercell[idx_shape].L3 = Not used
		//	supercell[idx_shape].L4 = Not used
		//	supercell[idx_shape].L5 = Not used
		//
		//	GSL defines sinc(x) = sin(pi*x)/(pi*x), hence the factor of pi in the denominator, which does not appear in the analytical expressions
		//
		//	NOTE: Shape body centre must always be positioned at 0 for two dimensions orthogonal to the quantum confinement dimension.
		//
		//	NOTE: Shape rotation is not supported in this case. Ensure shape alpha, beta, and gamma angles are set to 0.

		if ( supercell_type == 0 )
		{

			if ( lround( 0.1*supercell[idx_shape].L1 ) == 0 )
			{

				result *=  supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gx*supercell[idx_shape].L2 )/M_PI )
				*supercell[0].L2*gsl_sf_sinc( ( 0.5*Gy*supercell[0].L2 )/M_PI )
				*supercell[0].L3*gsl_sf_sinc( ( 0.5*Gz*supercell[0].L3 )/M_PI );

			}
			else if ( lround( 0.1*supercell[idx_shape].L1 ) == 1 )
			{

				result *=  supercell[0].L1*gsl_sf_sinc( ( 0.5*Gx*supercell[0].L1 )/M_PI )
				*supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gy*supercell[idx_shape].L2 )/M_PI )
				*supercell[0].L3*gsl_sf_sinc( ( 0.5*Gz*supercell[0].L3 )/M_PI );

			}
			else if ( lround( 0.1*supercell[idx_shape].L1 ) == 2 )
			{

				result *=  supercell[0].L1*gsl_sf_sinc( ( 0.5*Gx*supercell[0].L1 )/M_PI )
				*supercell[0].L2*gsl_sf_sinc( ( 0.5*Gy*supercell[0].L2 )/M_PI )
				*supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gz*supercell[idx_shape].L2 )/M_PI );

			}

		}
		else
		{

			if ( lround( 0.1*supercell[idx_shape].L1 ) == 0 )
			{

				result *=  supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gx*supercell[idx_shape].L2 )/M_PI )
				*supercell[0].SL2*gsl_sf_sinc( ( 0.5*Gy*supercell[0].SL2 )/M_PI )
				*supercell[0].SL3*gsl_sf_sinc( ( 0.5*Gz*supercell[0].SL3 )/M_PI );

			}
			else if ( lround( 0.1*supercell[idx_shape].L1 ) == 1 )
			{

				result *=  supercell[0].SL1*gsl_sf_sinc( ( 0.5*Gx*supercell[0].SL1 )/M_PI )
				*supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gy*supercell[idx_shape].L2 )/M_PI )
				*supercell[0].SL3*gsl_sf_sinc( ( 0.5*Gz*supercell[0].SL3 )/M_PI );

			}
			else if ( lround( 0.1*supercell[idx_shape].L1 ) == 2 )
			{

				result *=  supercell[0].SL1*gsl_sf_sinc( ( 0.5*Gx*supercell[0].SL1 )/M_PI )
				*supercell[0].SL2*gsl_sf_sinc( ( 0.5*Gy*supercell[0].SL2 )/M_PI )
				*supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gz*supercell[idx_shape].L2 )/M_PI );

			}

		}	

	}


	//////////////////////
	//// Quantum wire ////
	//////////////////////

	// Insert a quantum wire with quantum confinement in either the (100), (010), or (001) planes. Differs from a cuboid that extends across the supercell along a direction in that the quantum wire shape extends to the larger strain supercell boundaries for the relevant dimension when calculating the strain.

	if( supercell[idx_shape].shape == 18 )
	{

		// Note:
		//
		//	The body centre of the quantum wire is at position ( x, y, z ) = ( supercell[idx_shape].Ox, supercell[idx_shape].Oy, supercell[idx_shape].Oz )
		//
		//  The quantum confinement plane is given by supercell[idx_shape].L1
		//
		//	The quantum wire then extends along the quantum confinement plane between: 
		//  supercell[idx_shape].Ox/y/z - supercell[idx_shape].L2/2 <= x/y/z <= supercell[idx_shape].Ox/y/z + supercell[idx_shape].L2/2
		//  supercell[idx_shape].Ox/y/z - supercell[idx_shape].L3/2 <= x/y/z <= supercell[idx_shape].Ox/y/z + supercell[idx_shape].L3/2
		//
		//	supercell[idx_shape].L1 = Quantum confinement plane (0: (100),  1: (010),  2: (001) )
		//	supercell[idx_shape].L2 = Length of quantum wire along ([010] if L1 = 0, [100] if L1 = 1, [100] if L1 = 2)
		//	supercell[idx_shape].L3 = Length of quantum wire along ([001] if L1 = 0, [001] if L1 = 1, [010] if L1 = 2)
		//	supercell[idx_shape].L4 = Not used
		//	supercell[idx_shape].L5 = Not used
		//
		//	GSL defines sinc(x) = sin(pi*x)/(pi*x), hence the factor of pi in the denominator, which does not appear in the analytical expressions
		//
		//	NOTE: Shape body centre must always be positioned at 0 for the dimension normal to the quantum confinement plane.
		//
		//	NOTE: Shape rotation is not supported in this case. Ensure shape alpha, beta, and gamma angles are set to 0.

		if ( supercell_type == 0 )
		{

			if ( lround( 0.1*supercell[idx_shape].L1 ) == 0 )
			{

				result *=  supercell[0].L1*gsl_sf_sinc( ( 0.5*Gx*supercell[0].L1 )/M_PI )
				*supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gy*supercell[idx_shape].L2 )/M_PI )
				*supercell[idx_shape].L3*gsl_sf_sinc( ( 0.5*Gz*supercell[idx_shape].L3 )/M_PI );

			}
			else if ( lround( 0.1*supercell[idx_shape].L1 ) == 1 )
			{

				result *=  supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gx*supercell[idx_shape].L2 )/M_PI )
				*supercell[0].L2*gsl_sf_sinc( ( 0.5*Gy*supercell[0].L2 )/M_PI )
				*supercell[idx_shape].L3*gsl_sf_sinc( ( 0.5*Gz*supercell[idx_shape].L3 )/M_PI );

			}
			else if ( lround( 0.1*supercell[idx_shape].L1 ) == 2 )
			{

				result *=  supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gx*supercell[idx_shape].L2 )/M_PI )
				*supercell[idx_shape].L3*gsl_sf_sinc( ( 0.5*Gy*supercell[idx_shape].L3 )/M_PI )
				*supercell[0].L3*gsl_sf_sinc( ( 0.5*Gz*supercell[0].L3 )/M_PI );

			}

		}
		else
		{

			if ( lround( 0.1*supercell[idx_shape].L1 ) == 0 )
			{

				result *=  supercell[0].SL1*gsl_sf_sinc( ( 0.5*Gx*supercell[0].SL1 )/M_PI )
				*supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gy*supercell[idx_shape].L2 )/M_PI )
				*supercell[idx_shape].L3*gsl_sf_sinc( ( 0.5*Gz*supercell[idx_shape].L3 )/M_PI );

			}
			else if ( lround( 0.1*supercell[idx_shape].L1 ) == 1 )
			{

				result *=  supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gx*supercell[idx_shape].L2 )/M_PI )
				*supercell[0].SL2*gsl_sf_sinc( ( 0.5*Gy*supercell[0].SL2 )/M_PI )
				*supercell[idx_shape].L3*gsl_sf_sinc( ( 0.5*Gz*supercell[idx_shape].L3 )/M_PI );

			}
			else if ( lround( 0.1*supercell[idx_shape].L1 ) == 2 )
			{

				result *=  supercell[idx_shape].L2*gsl_sf_sinc( ( 0.5*Gx*supercell[idx_shape].L2 )/M_PI )
				*supercell[idx_shape].L3*gsl_sf_sinc( ( 0.5*Gy*supercell[idx_shape].L3 )/M_PI )
				*supercell[0].SL3*gsl_sf_sinc( ( 0.5*Gz*supercell[0].SL3 )/M_PI );

			}

		}	

	}


	return result;

}


/////////////////////
//// Subroutines ////
/////////////////////

complex double characteristic_function_fourier_pyramid_segment( double Gx0, double Gy0, double Gz0, int idx_side, int N_side, struct supercell_geometry *supercell, int idx_shape )
{

	double alpha, x0, y0, h, Gx, Gy, Gz, G_mag;
	complex double result = 1 + 0*I;


	// Set polar angle alpha

	if( N_side == 0 ) // Special case of rectangular-based pyramid
	{

		if( ( idx_side == 0 ) || ( idx_side == 2 ) )
		{

			alpha = atan( supercell[idx_shape].L2/supercell[idx_shape].L1 );

		}
		else
		{

			alpha = atan( supercell[idx_shape].L1/supercell[idx_shape].L2 );

		}
	}
	else
	{

		alpha = M_PI/( ( double ) N_side );

	}


	// Set lengths x0, y0 and h

	if( ( N_side == 0 ) && ( ( idx_side == 0 ) || ( idx_side == 2 ) ) )
	{

		x0 = supercell[idx_shape].L1;
		y0 = supercell[idx_shape].L2;

	}
	else if( N_side == 0 && ( ( idx_side == 1 ) || ( idx_side == 3 ) ) )
	{

		x0 = supercell[idx_shape].L2;
		y0 = supercell[idx_shape].L1;

	}
	else
	{

		x0 =            supercell[idx_shape].L1;
		y0 = tan(alpha)*supercell[idx_shape].L1;

	}

	if ( N_side == 0 )
	{

		h = supercell[idx_shape].L3;

	}
	else
	{

		h = supercell[idx_shape].L2;

	}

	


	// Set wave vector

	Gx =      cos( ( 2.0*alpha*( ( double ) idx_side ) ) )*Gx0 + sin( ( 2.0*alpha*( ( double ) idx_side ) ) )*Gy0;
	Gy = -1.0*sin( ( 2.0*alpha*( ( double ) idx_side ) ) )*Gx0 + cos( ( 2.0*alpha*( ( double ) idx_side ) ) )*Gy0;
	Gz = Gz0;

	G_mag = sqrt( pow( Gx, 2.0 ) + pow( Gy, 2.0 ) + pow( Gz, 2.0 ) ); // Magnitude of wave vector


	// Evaluate Fourier transform of characteristic function


	if( fabs( G_mag ) < EPS ) // G = 0
	{

		result *= ( x0*y0*h )/3.0; // Volume

	}
	else if( ( fabs( Gx ) < EPS ) && ( fabs( Gy ) < EPS ) ) // Gx = Gy = 0, Gz != 0
	{

		result *= ( ( I*2.0*y0*x0 )/Gz )*( ( 1.0 - I*Gz*h - cexp( -1.0*I*Gz*h ) )/( Gz*Gz*h*h ) - 0.5 );

	}
	else if( ( fabs( Gx ) < EPS ) && ( fabs( Gz ) < EPS ) ) // Gx = Gz = 0, Gy != 0
	{

		result *= ( ( 2.0*h*x0 )/( Gy*Gy*y0 ) )*( 1.0 - sin( Gy*y0 )/( Gy*y0 ) );

	}
	else if( ( fabs( Gy ) < EPS ) && ( fabs( Gz ) < EPS ) ) // Gy = Gz = 0, Gx != 0
	{

		result *= ( ( 2.0*h*y0 )/( Gx*Gx*x0 ) )*( -1.0 - cexp( -1.0*I*Gx*x0 ) - ( 2.0*I )/( Gx*x0 ) + 2.0*I*cexp( -1.0*I*Gx*x0 )/( Gx*x0 ) );

	}
	else if( fabs( Gx ) < EPS ) // Gx = 0, Gy != 0, Gz != 0
	{

		if( fabs( Gy*y0 - Gz*h ) < EPS )
		{

			result *= ( ( -1.0*I*x0 )/( Gz*Gy ) )*( cexp( -1.0*I*Gz*h )*( -1.0*I + ( cexp( I*( Gy*y0 + Gz*h ) ) - 1.0 )/( Gy*y0 + Gz*h ) ) - 2.0*( cos( Gy*y0 ) - 1.0 )/( Gy*y0 ) );

		}
		else if( fabs( Gy*y0 + Gz*h ) < EPS )
		{

			result *= ( ( -1.0*I*x0 )/( Gz*Gy ) )*( cexp( -1.0*I*Gz*h )*( ( cexp( -1.0*I*( Gy*y0 - Gz*h ) ) - 1.0 )/( Gy*y0 - Gz*h ) + I ) - 2.0*( cos( Gy*y0 ) - 1.0 )/( Gy*y0 ) );

		}
		else
		{

			result *= ( ( -1.0*I*x0 )/( Gz*Gy ) )*( cexp( -1.0*I*Gz*h )*( ( cexp( -1.0*I*( Gy*y0 - Gz*h ) ) - 1.0 )/( Gy*y0 - Gz*h ) + ( cexp( I*( Gy*y0 + Gz*h ) ) - 1.0 )/( Gy*y0 + Gz*h ) ) - 2.0*( cos( Gy*y0 ) - 1.0 )/( Gy*y0 ) );

		}

	}
	else if( fabs( Gy ) < EPS ) // // Gy = 0, Gx != 0, Gz != 0
	{

		if( fabs( Gx - ( Gz*h )/x0 ) < EPS )
		{

			result *= ( ( 2.0*I*y0 )/( Gz*x0 ) )*( 0.5*cexp( -1.0*I*Gz*h )*x0*x0 - ( ( cexp( -1.0*I*Gx*x0 )*( 1.0 + I*Gx*x0 ) - 1.0 )/( Gx*Gx ) ) );

		}
		else
		{

			result *= ( ( 2.0*I*y0 )/( Gz*x0 ) )*( cexp( -1.0*I*Gz*h )*( cexp( -1.0*I*( Gx*x0 - Gz*h ) )*( 1.0 + I*( Gx*x0 - Gz*h ) ) - 1.0 )/pow( Gx - ( Gz*h )/x0, 2.0 ) - ( ( cexp( -1.0*I*Gx*x0 )*( 1.0 + I*Gx*x0 ) - 1.0 )/( Gx*Gx ) ) );

		}

	}
	else if( fabs( Gz ) < EPS ) // Gz = 0, Gx != 0, Gy != 0
	{

		if( fabs( Gx - ( Gy*y0 )/x0 ) < EPS )
		{

			result *= ( h/Gy )*(  1.0/( Gx + ( Gy*y0 )/x0 ) - ( I/( x0*pow( Gx + ( Gy*y0 )/x0, 2.0 ) ) )*( cexp( -1.0*I*x0*( Gx + ( Gy*y0 )/x0 ) ) - 1.0 ) - 0.5*I*x0 );

		}
		else if( fabs( Gx + ( Gy*y0 )/x0 ) < EPS )
		{

			result *= ( h/Gy )*( -1.0/( Gx - ( Gy*y0 )/x0 ) + ( I/( x0*pow( Gx - ( Gy*y0 )/x0, 2.0 ) ) )*( cexp( -1.0*I*x0*( Gx - ( Gy*y0 )/x0 ) ) - 1.0 ) + 0.5*I*x0 );

		}
		else
		{

			result *= ( h/Gy )*(  1.0/( Gx + ( Gy*y0 )/x0 ) - 1.0/( Gx - ( Gy*y0 )/x0 ) - ( I/( x0*pow( Gx + ( Gy*y0 )/x0, 2.0 ) ) )*( cexp( -1.0*I*x0*( Gx + ( Gy*y0 )/x0 ) ) - 1.0 )
				   		     											 + ( I/( x0*pow( Gx - ( Gy*y0 )/x0, 2.0 ) ) )*( cexp( -1.0*I*x0*( Gx - ( Gy*y0 )/x0 ) ) - 1.0 ) );

		}
	}
	else // Gx != 0, Gy != 0, Gz != 0
	{

		result *= ( -1.0/( Gy*Gz ) )*( I_e0( x0, -1.0*Gx - ( Gy*y0 )/x0 + ( Gz*h )/x0 )*cexp( -1.0*I*Gz*h )
							    - I_e0( x0, -1.0*Gx + ( Gy*y0 )/x0 + ( Gz*h )/x0 )*cexp( -1.0*I*Gz*h )
							    - I_e0( x0, -1.0*Gx - ( Gy*y0 )/x0 )
							    + I_e0( x0, -1.0*Gx + ( Gy*y0 )/x0 ) );

	}

	return result;

}

complex double characteristic_function_fourier_prism_segment( double Gx0, double Gy0, double Gz0, int idx_side, int N_side, struct supercell_geometry *supercell, int idx_shape )
{

	double alpha, x0, h, Gx, Gy, Gz, G_mag;
	complex double result = 1 + 0*I;


	alpha = M_PI/( ( double ) N_side );


	// Set lengths x0 and h

	x0 =  	supercell[idx_shape].L1;
	h = 	supercell[idx_shape].L2;


	// Set wave vector

	Gx =      cos( ( 2.0*M_PI*( ( double ) idx_side ) )/N_side )*Gx0 + sin( ( 2.0*M_PI*( ( double ) idx_side ) )/N_side )*Gy0;
	Gy = -1.0*sin( ( 2.0*M_PI*( ( double ) idx_side ) )/N_side )*Gx0 + cos( ( 2.0*M_PI*( ( double ) idx_side ) )/N_side )*Gy0;
	Gz = Gz0;

	G_mag = sqrt( pow( Gx, 2.0 ) + pow( Gy, 2.0 ) + pow( Gz, 2.0 ) ); // Magnitude of wave vector


	// Evaluate Fourier transform of characteristic function


	if( fabs( G_mag ) < EPS ) // G = 0
	{

		result *= ( x0*x0*h*tan( alpha ) ); // Volume

	}
	else if( ( fabs( Gx ) < EPS ) && ( fabs( Gy ) < EPS ) ) // Gx = Gy = 0, Gz != 0
	{

		result *= 2.0*x0*x0*tan( alpha )*sin( 0.5*Gz*h )/Gz;

	}
	else if( ( fabs( Gx ) < EPS ) && ( fabs( Gz ) < EPS ) ) // Gx = Gz = 0, Gy != 0
	{

		result *= ( 2.0*h/( Gy*Gy*tan( alpha ) ) )*( 1.0 - cos( Gy*x0*tan( alpha ) ) );

	}
	else if( ( fabs( Gy ) < EPS ) && ( fabs( Gz ) < EPS ) ) // Gy = Gz = 0, Gx != 0
	{

		result *= ( 2.0*h*tan( alpha )/Gx )*( 1.0*I*x0*cexp( -1.0*I*Gx*x0 ) + ( 1.0/Gx )*( cexp( -1.0*I*Gx*x0 ) - 1.0 ) );

	}
	else if( fabs( Gx ) < EPS ) // Gx = 0, Gy != 0, Gz != 0
	{

		result *= ( 4.0*sin( 0.5*Gz*h )/( Gy*Gy*Gz*tan( alpha ) ) )*( 1.0 - cos( Gy*x0*tan( alpha ) ) );

	}
	else if( fabs( Gy ) < EPS ) // // Gy = 0, Gx != 0, Gz != 0
	{

		result *= ( 4.0*sin( 0.5*Gz*h )*tan( alpha )/( Gz*Gx ) )*( 1.0*I*x0*cexp( -1.0*I*Gx*x0 ) + ( 1.0/Gx )*( cexp( -1.0*I*Gx*x0 ) - 1.0 ) );

	}
	else if( fabs( Gz ) < EPS ) // Gz = 0, Gx != 0, Gy != 0
	{

		if ( fabs( Gy*tan( alpha ) + Gx ) < EPS )
		{

			result *= ( h/Gy )*( 1.0*I*x0 - ( 1.0/( -Gy*tan( alpha ) + Gx ) )*( 1.0 - cexp( 1.0*I*Gy*x0*tan( alpha ) - 1.0*I*Gx*x0 ) ) );

		}
		else if ( fabs( -Gy*tan( alpha ) + Gx ) < EPS )
		{

			result *= ( h/Gy )*( ( 1.0/( Gy*tan( alpha ) + Gx ) )*( 1.0 - cexp( -1.0*I*Gy*x0*tan( alpha ) - 1.0*I*Gx*x0 ) ) - 1.0*I*x0 );

		}
		else
		{

			result *= ( h/Gy )*( ( 1.0/( Gy*tan( alpha ) + Gx ) )*( 1.0 - cexp( -1.0*I*Gy*x0*tan( alpha ) - 1.0*I*Gx*x0 ) ) - ( 1.0/( -Gy*tan( alpha ) + Gx ) )*( 1.0 - cexp( 1.0*I*Gy*x0*tan( alpha ) - 1.0*I*Gx*x0 ) ) );

		}

	}
	else // Gx != 0, Gy != 0, Gz != 0
	{

		if ( fabs( Gy*tan( alpha ) + Gx ) < EPS )
		{

			result *= ( 2.0/( Gy*Gz ) )*sin( 0.5*Gz*h )*( 1.0*I*x0 - ( 1.0/( -Gy*tan( alpha ) + Gx ) )*( 1.0 - cexp( 1.0*I*Gy*x0*tan( alpha ) - 1.0*I*Gx*x0 ) ) );

		}
		else if ( fabs( -Gy*tan( alpha ) + Gx ) < EPS )
		{

			result *= ( 2.0/( Gy*Gz ) )*sin( 0.5*Gz*h )*( ( 1.0/( Gy*tan( alpha ) + Gx ) )*( 1.0 - cexp( -1.0*I*Gy*x0*tan( alpha ) - 1.0*I*Gx*x0 ) ) - 1.0*I*x0 );

		}
		else
		{

			result *= ( 2.0/( Gy*Gz ) )*sin( 0.5*Gz*h )*( ( 1.0/( Gy*tan( alpha ) + Gx ) )*( 1.0 - cexp( -1.0*I*Gy*x0*tan( alpha ) - 1.0*I*Gx*x0 ) ) - ( 1.0/( -Gy*tan( alpha ) + Gx ) )*( 1.0 - cexp( 1.0*I*Gy*x0*tan( alpha ) - 1.0*I*Gx*x0 ) ) );

		}

	}

	return result;

}
