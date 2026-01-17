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

complex double strain_fourier_coefficient( int component, int isotropic, int idx_shape, int idx_basis, int G_prime, double Gx, double Gy, double Gz, struct supercell_geometry *supercell, struct material_params *parameters, struct characteristic_functions characteristic_function_G, int correction_order, struct convolution_char_strain *convolutions_char_strain_correction )
{

	int idx_row = -1, idx_col = -1;
	double Gi = 0.0, Gj = 0.0, G, eps0, sum;
	complex double result;

	// Assign the wave vector components Gi and Gj based on which component of the strain tensor is being computed

	G = sqrt( pow( Gx, 2.0 ) + pow( Gy, 2.0 ) + pow( Gz, 2.0 ) ); // Magnitude of wave vector

	switch( component )
	{

		case 0: // xx

			idx_row = 0;
			idx_col = 0;

			Gi = Gx;
			Gj = Gx;

			break;

		case 1: // yy

			idx_row = 1;
			idx_col = 1;

			Gi = Gy;
			Gj = Gy;

			break;

		case 2: // zz

			idx_row = 2;
			idx_col = 2;

			Gi = Gz;
			Gj = Gz;

			break;

		case 3: // yz

			idx_row = 1;
			idx_col = 2;

			Gi = Gy;
			Gj = Gz;

			break;

		case 4: // xz

			idx_row = 0;
			idx_col = 2;

			Gi = Gx;
			Gj = Gz;

			break;

		case 5: // xy

			idx_row = 0;
			idx_col = 1;

			Gi = Gx;
			Gj = Gy;

			break;

	}


	// Lattice mismatch between current shape (material) and the matrix material

	eps0 = ( parameters[supercell[idx_shape].embedded].a - parameters[idx_shape].a )/parameters[0].a;


	// Then the following computes the Fourier coefficient of the chosen component of the strain tensor for the supplied wave vector

	
	if ( correction_order == -1 ) // Only occurs for the initial contribution in the case of position-dependent elastic constants
	{

		// Initial contribution for diagonal elements in case of position-dependent elastic cosntants
		result = eps0*characteristic_function_G.chi[idx_basis]*kronecker_delta( idx_row, idx_col );

	}
	else if ( correction_order == 0 ) // Occurs for cases of both position-independent and position-dependent elastic constants
	{

		result = eps0*characteristic_function_G.chi[idx_basis];
		
		if ( include_shape_elastic_constants == 0 ) // Only option executed in case of position-independent elastic constants
		{

			if( fabs( G ) < EPS ) // |G| = 0
			{

				if( isotropic == 1 )
				{

					// Isotropic approximation to the elastic tensor

					result *= kronecker_delta( idx_row, idx_col )*( 1.0 - ( 3.0*parameters[0].C12 + 2.0*parameters[0].C44 )/( 3.0*parameters[0].C12 + 6.0*parameters[0].C44 ) );

				}
				else
				{

					// Full, anisotropic elastic tensor

					sum = 1.0/( parameters[0].C44 + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )/3.0 );


					result *=  kronecker_delta( idx_row, idx_col )*( 1.0 - ( 0.5*( parameters[0].C11 + 2.0*parameters[0].C12 )/3.0 )*( 2.0/( parameters[0].C44 + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )/3.0 ) )/( 1.0 + ( parameters[0].C12 + parameters[0].C44 )*sum ) );

				}

			}
			else // |G| != 0
			{

				if( isotropic == 1 )
				{

					// Isotropic approximation: C11 - C12 - 2*C44 = 0

					result *= kronecker_delta( idx_row, idx_col ) - ( ( ( 3.0*parameters[0].C12 + 2.0*parameters[0].C44 )*Gi*Gj )
																/( (     parameters[0].C12 + 2.0*parameters[0].C44 )*G *G  ) );

				}
				else
				{

					// Full, anisotropic elastic tensor

					sum = ( ( Gx*Gx )/( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gx*Gx ) )
						+ ( ( Gy*Gy )/( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gy*Gy ) )
						+ ( ( Gz*Gz )/( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gz*Gz ) );


					result *=  kronecker_delta( idx_row, idx_col )       - ( 0.5*( parameters[0].C11 + 2.0*parameters[0].C12 )
							*( 2.0*parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*( Gi*Gi + Gj*Gj ) )*Gi*Gj )
							/( ( 1.0 + ( parameters[0].C12 + parameters[0].C44 )*sum )*( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gi*Gi )
																		*( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gj*Gj ) );

				}

			}

		}
		else // Executed in case of position-dependent elastic constants to calculate 0th order contribution
		{

			if( fabs( G ) < EPS ) // |G| = 0
			{

				if( isotropic == 1 )
				{

					// Isotropic approximation to the elastic tensor

					result *= kronecker_delta( idx_row, idx_col )*(-1*( 3.0*parameters[idx_shape].C12 + 2.0*parameters[idx_shape].C44 )/( 3.0*parameters[0].C12 + 6.0*parameters[0].C44 ) );

				}
				else
				{

					// Full, anisotropic elastic tensor

					sum = 1.0/( parameters[0].C44 + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )/3.0 );


					result *=  kronecker_delta( idx_row, idx_col )*(( -0.5*( parameters[idx_shape].C11 + 2.0*parameters[idx_shape].C12 )/3.0 )*( 2.0/( parameters[0].C44 + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )/3.0 ) )/( 1.0 + ( parameters[0].C12 + parameters[0].C44 )*sum ) );

				}

			}
			else // |G| != 0
			{

				if( isotropic == 1 )
				{

					// Isotropic approximation: C11 - C12 - 2*C44 = 0

					result *=  -( ( 3.0*parameters[idx_shape].C12 + 2.0*parameters[idx_shape].C44 )*Gi*Gj )
																/( (     parameters[0].C12 + 2.0*parameters[0].C44 )*G *G  );

				}
				else
				{

					// Full, anisotropic elastic tensor

					sum = ( ( Gx*Gx )/( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gx*Gx ) )
						+ ( ( Gy*Gy )/( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gy*Gy ) )
						+ ( ( Gz*Gz )/( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gz*Gz ) );


					result *= -( 0.5*( parameters[idx_shape].C11 + 2.0*parameters[idx_shape].C12 )
							*( 2.0*parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*( Gi*Gi + Gj*Gj ) )*Gi*Gj )
							/( ( 1.0 + ( parameters[0].C12 + parameters[0].C44 )*sum )*( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gi*Gi )
																		*( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gj*Gj ) );

				}

			}

		}

	}
	else // Executed in case of position-dependent elastic constants for correction orders > 0
	{

		switch( component )
		{

			case 0: // xx

				result = convolutions_char_strain_correction->char_xx[idx_basis];

				break;

			case 1: // yy

				result = convolutions_char_strain_correction->char_yy[idx_basis];

				break;

			case 2: // zz

				result = convolutions_char_strain_correction->char_zz[idx_basis];

				break;

			case 3: // yz

				result = convolutions_char_strain_correction->char_yz[idx_basis];

				break;

			case 4: // xz

				result = convolutions_char_strain_correction->char_xz[idx_basis];

				break;

			case 5: // xy

				result = convolutions_char_strain_correction->char_xy[idx_basis];

				break;

		}

		if( fabs( G ) < EPS ) // |G| = 0
		{

			if( isotropic == 1 )
			{

				// Isotropic approximation to the elastic tensor

				if ( component == 0 || component == 1 || component == 2)
				{

					result *= -( 1.0/( 3.0*parameters[0].C44 ) )*( parameters[idx_shape].C11 - parameters[0].C11 + 2.0*parameters[idx_shape].C12 - 2.0*parameters[0].C12 );

				}
				else
				{

					result *=  0;

				}

			}
			else
			{

				// Full, anisotropic elastic tensor

				sum = 1.0/( parameters[0].C44 + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )/3.0 );

				if ( component == 0 || component == 1 || component == 2)
				{

					result *= -( sum/3.0 )*( parameters[idx_shape].C11 - parameters[0].C11 + 2.0*parameters[idx_shape].C12 - 2.0*parameters[0].C12 );

				}
				else
				{

					result *=  0;

				}

			}

		}
		else // |G| != 0
		{

			// Isotropic & anisotropic (isotropic may not reduce nicely)

			sum = ( ( Gx*Gx )/( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gx*Gx ) )
				+ ( ( Gy*Gy )/( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gy*Gy ) )
				+ ( ( Gz*Gz )/( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gz*Gz ) );


			switch( component )
			{

				case 0: // xx

					result *= -1.0*( Gi*( ( ( parameters[idx_shape].C11 - parameters[0].C11 )*Gx*greens_function_tensor( parameters, 0, Gx, Gy, Gz, G, sum ) ) + ( ( parameters[idx_shape].C12 - parameters[0].C12 )*Gy*greens_function_tensor( parameters, 5, Gx, Gy, Gz, G, sum ) ) + ( ( parameters[idx_shape].C12 - parameters[0].C12 )*Gz*greens_function_tensor( parameters, 4, Gx, Gy, Gz, G, sum ) ) ) );

				break;

				case 1: // yy

					result *= -1.0*( Gi*( ( ( parameters[idx_shape].C12 - parameters[0].C12 )*Gx*greens_function_tensor( parameters, 5, Gx, Gy, Gz, G, sum ) ) + ( ( parameters[idx_shape].C11 - parameters[0].C11 )*Gy*greens_function_tensor( parameters, 1, Gx, Gy, Gz, G, sum ) ) + ( ( parameters[idx_shape].C12 - parameters[0].C12 )*Gz*greens_function_tensor( parameters, 3, Gx, Gy, Gz, G, sum ) ) ) );

				break;

				case 2: // zz

					result *= -1.0*( Gi*( ( ( parameters[idx_shape].C12 - parameters[0].C12 )*Gx*greens_function_tensor( parameters, 4, Gx, Gy, Gz, G, sum ) ) + ( ( parameters[idx_shape].C12 - parameters[0].C12 )*Gy*greens_function_tensor( parameters, 3, Gx, Gy, Gz, G, sum ) ) + ( ( parameters[idx_shape].C11 - parameters[0].C11 )*Gz*greens_function_tensor( parameters, 2, Gx, Gy, Gz, G, sum ) ) ) );

				break;

				case 3: // yz

					result *= -0.5*( parameters[idx_shape].C44 - parameters[0].C44 )*( Gi*( ( Gz*greens_function_tensor( parameters, 3, Gx, Gy, Gz, G, sum ) ) + ( Gy*greens_function_tensor( parameters, 2, Gx, Gy, Gz, G, sum ) ) )
																					+ Gj*( ( Gz*greens_function_tensor( parameters, 1, Gx, Gy, Gz, G, sum ) ) + ( Gy*greens_function_tensor( parameters, 3, Gx, Gy, Gz, G, sum ) ) ) );

				break;

				case 4: // xz

					result *= -0.5*( parameters[idx_shape].C44 - parameters[0].C44 )*( Gi*( ( Gz*greens_function_tensor( parameters, 4, Gx, Gy, Gz, G, sum ) ) + ( Gx*greens_function_tensor( parameters, 2, Gx, Gy, Gz, G, sum ) ) )
																					+ Gj*( ( Gz*greens_function_tensor( parameters, 0, Gx, Gy, Gz, G, sum ) ) + ( Gx*greens_function_tensor( parameters, 4, Gx, Gy, Gz, G, sum ) ) ) );

				break;

				case 5: // xy

					result *= -0.5*( parameters[idx_shape].C44 - parameters[0].C44 )*( Gi*( ( Gy*greens_function_tensor( parameters, 5, Gx, Gy, Gz, G, sum ) ) + ( Gx*greens_function_tensor( parameters, 1, Gx, Gy, Gz, G, sum ) ) )
																					+ Gj*( ( Gy*greens_function_tensor( parameters, 0, Gx, Gy, Gz, G, sum ) ) + ( Gx*greens_function_tensor( parameters, 5, Gx, Gy, Gz, G, sum ) ) ) );

				break;

			}

		}

	}


	return result;

}


// Greeen's function tensor
double greens_function_tensor( struct material_params *parameters, int component, double Gx, double Gy, double Gz, double G, double sum )
{

	int idx_row, idx_col;
	double Gi, Gj, result;

	switch( component )
	{

		case 0: // xx

			idx_row = 0;
			idx_col = 0;

			Gi = Gx;
			Gj = Gx;

			break;

		case 1: // yy

			idx_row = 1;
			idx_col = 1;

			Gi = Gy;
			Gj = Gy;

			break;

		case 2: // zz

			idx_row = 2;
			idx_col = 2;

			Gi = Gz;
			Gj = Gz;

			break;

		case 3: // yz

			idx_row = 1;
			idx_col = 2;

			Gi = Gy;
			Gj = Gz;

			break;

		case 4: // xz

			idx_row = 0;
			idx_col = 2;

			Gi = Gx;
			Gj = Gz;

			break;

		case 5: // xy

			idx_row = 0;
			idx_col = 1;

			Gi = Gx;
			Gj = Gy;

			break;

	}

	result = ( kronecker_delta( idx_row, idx_col ) / ( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gi*Gi ) ) - ( ( parameters[0].C12 + parameters[0].C44 )*Gi*Gj / ( ( 1.0 + ( parameters[0].C12 + parameters[0].C44 )*sum )*( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gi*Gi )*( parameters[0].C44*G*G + ( parameters[0].C11 - parameters[0].C12 - 2.0*parameters[0].C44 )*Gj*Gj ) ) );

	return result;

}