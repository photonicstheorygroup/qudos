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

double interpolate_quaternary( int bowing, struct supercell_geometry supercell, struct interpolate_parameters alloy_parameters )
{

	double P_ABC, P_ABD, P_ACD, P_BCD, result = 0.0;

	// bowing = 0 ==> Linear interpolation
	// bowing = 1 ==> Interpolation with ternary and quaternary bowing parameters

	if( bowing == 0 )
	{

		switch( supercell.material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				result =         supercell.x  *supercell.y*alloy_parameters.P_AC +         supercell.x  *( 1.0 - supercell.y )*alloy_parameters.P_AD
			  	  	  + ( 1.0 - supercell.x )*supercell.y*alloy_parameters.P_BC + ( 1.0 - supercell.x )*( 1.0 - supercell.y )*alloy_parameters.P_BD;

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				result = supercell.x*alloy_parameters.P_AD + supercell.y*alloy_parameters.P_BD + ( 1.0 - supercell.x - supercell.y )*alloy_parameters.P_CD;

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				result = supercell.x*alloy_parameters.P_AB + supercell.y*alloy_parameters.P_AC + ( 1.0 - supercell.x - supercell.y )*alloy_parameters.P_AD;

				break;

		}

	}
	else if( bowing == 1 )
	{

		switch( supercell.material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				// There is a danger of dividing by zero and obtaining NaN for this alloy type ==> Deal with binary cases separately

					if( ( fabs( supercell.x - 1.0 ) < EPS ) && ( fabs( supercell.y - 1.0 ) < EPS ) ) // (x,y) = (1,1) ==> A(x)B(1-x)C(y)D(1-y) = AC
				{

					result = alloy_parameters.P_AC;

				}
				else if( ( fabs( supercell.x - 1.0 ) < EPS ) && ( fabs( supercell.y ) <       EPS ) ) // (x,y) = (1,0) ==> A(x)B(1-x)C(y)D(1-y) = AD
				{

					result = alloy_parameters.P_AD;

				}
				else if( ( fabs( supercell.x ) <       EPS ) && ( fabs( supercell.y - 1.0 ) < EPS ) ) // (x,y) = (0,1) ==> A(x)B(1-x)C(y)D(1-y) = BC
				{

					result = alloy_parameters.P_BC;

				}
				else if( ( fabs( supercell.x ) <       EPS ) && ( fabs( supercell.y ) <       EPS ) ) // (x,y) = (0,0) ==> A(x)B(1-x)C(y)D(1-y) = BD
				{

					result = alloy_parameters.P_BD;

				}
				else
				{


				 	P_ABC = supercell.x*alloy_parameters.P_AC + ( 1.0 - supercell.x )*alloy_parameters.P_BC - supercell.x*( 1.0 - supercell.x )*alloy_parameters.B_ABC;
					P_ABD = supercell.x*alloy_parameters.P_AD + ( 1.0 - supercell.x )*alloy_parameters.P_BD - supercell.x*( 1.0 - supercell.x )*alloy_parameters.B_ABD;
				 	P_ACD = supercell.y*alloy_parameters.P_AC + ( 1.0 - supercell.y )*alloy_parameters.P_AD - supercell.y*( 1.0 - supercell.y )*alloy_parameters.B_ACD;
				 	P_BCD = supercell.y*alloy_parameters.P_BC + ( 1.0 - supercell.y )*alloy_parameters.P_BD - supercell.y*( 1.0 - supercell.y )*alloy_parameters.B_BCD;

					result  = supercell.x*( 1.0 - supercell.x )*( supercell.y*P_ABC + ( 1.0 - supercell.y )*P_ABD )
					   	   + supercell.y*( 1.0 - supercell.y )*( supercell.x*P_ACD + ( 1.0 - supercell.x )*P_BCD );

					result /= supercell.x*( 1.0 - supercell.x ) + supercell.y*( 1.0 - supercell.y );

					result -= supercell.x*( 1.0 - supercell.x ) * supercell.y*( 1.0 - supercell.y )*alloy_parameters.B_ABCD;

				}

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				result = supercell.x                                                *alloy_parameters.P_AD
				   	  +             supercell.y                                    *alloy_parameters.P_BD
				   	  +                         ( 1.0 - supercell.x - supercell.y )*alloy_parameters.P_CD
				   	  - supercell.x*supercell.y                                    *alloy_parameters.B_ABD
				   	  - supercell.x            *( 1.0 - supercell.x - supercell.y )*alloy_parameters.B_ACD
				   	  -             supercell.y*( 1.0 - supercell.x - supercell.y )*alloy_parameters.B_BCD
				   	  - supercell.x*supercell.y*( 1.0 - supercell.x - supercell.y )*alloy_parameters.B_ABCD;

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				result = supercell.x                                                *alloy_parameters.P_AB
				  	  +             supercell.y                                    *alloy_parameters.P_AC
				  	  +                         ( 1.0 - supercell.x - supercell.y )*alloy_parameters.P_AD
			       	  - supercell.x*supercell.y                                    *alloy_parameters.B_ABC
				  	  - supercell.x            *( 1.0 - supercell.x - supercell.y )*alloy_parameters.B_ABD
				  	  -             supercell.y*( 1.0 - supercell.x - supercell.y )*alloy_parameters.B_ACD
				  	  - supercell.x*supercell.y*( 1.0 - supercell.x - supercell.y )*alloy_parameters.B_ABCD;

				break;

		}

	}


	return result;

}
