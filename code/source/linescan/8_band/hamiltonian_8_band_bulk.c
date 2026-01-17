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

complex double hamiltonian_8_band_bulk( int idx_row, int idx_col, double x, double y, double z, struct material_params parameters, struct strain_tensor epsilon, complex double phi_piezo, complex double phi_F )
{

	complex double result = 0.0 + I*0.0;

	switch( idx_row + 1 )
	{

		// Row 1

		case 1 :

			switch( idx_col + 1 )
			{

				// Row 1, column 1

				case 1 :

					result = parameters.E_cb + parameters.ac*( epsilon.xx[0] + epsilon.yy[0] + epsilon.zz[0] ) + phi_piezo + phi_F;

					break;


				// Row 1, column 2

				case 2 :

					result = 0.0 + I*0.0;

					break;


				// Row 1, column 3

				case 3 :

					result = 0.0 + I*0.0;

					break;


				// Row 1, column 4

				case 4 :

					result = 0.0 + I*0.0;

					break;


				// Row 1, column 5

				case 5 :

					result = 0.0 + I*0.0;

					break;


				// Row 1, column 6

				case 6 :

					result = 0.0 + I*0.0;

					break;


				// Row 1, column 7

				case 7 :

					result = 0.0 + I*0.0;

					break;


				// Row 1, column 8

				case 8 :

					result = 0.0 + I*0.0;

					break;

			}

			break;


		// Row 2

		case 2 :

			switch( idx_col + 1 )
			{

				// Row 2, column 1

				case 1 :

					result = 0.0 + I*0.0;

					break;


				// Row 2, column 2

				case 2 :

					result = parameters.E_hh + parameters.av*( epsilon.xx[0] + epsilon.yy[0] + epsilon.zz[0] ) + 0.5*parameters.b*( epsilon.xx[0] + epsilon.yy[0] - 2.0*epsilon.zz[0] ) + phi_piezo + phi_F;

					break;


				// Row 2, column 3

				case 3 :

					result = -1.0*parameters.d*( epsilon.xz[0] - I*epsilon.yz[0] );

					break;


				// Row 2, column 4

				case 4 :

					result = sqrt( 0.5 )*parameters.d*( epsilon.xz[0] - I*epsilon.yz[0] );

					break;


				// Row 2, column 5

				case 5 :

					result = 0.0 + I*0.0;

					break;


				// Row 2, column 6

				case 6 :

					result = 0.0 + I*0.0;

					break;


				// Row 2, column 7

				case 7 :

					result = 0.5*sqrt( 3.0 )*parameters.b*( epsilon.xx[0] - epsilon.yy[0] ) - I*parameters.d*epsilon.xy[0];

					break;


				// Row 2, column 8

				case 8 :

					result = 0.5*sqrt( 6.0 )*parameters.b*( epsilon.xx[0] - epsilon.yy[0] ) - I*sqrt( 2.0 )*parameters.d*epsilon.xy[0];

					break;

			}

			break;


		// Row 3

		case 3 :

			switch( idx_col + 1 )
			{

				// Row 3, column 1

				case 1 :

					result = 0.0 + I*0.0;

					break;


				// Row 3, column 2

				case 2 :

					result = -1.0*parameters.d*( epsilon.xz[0] + I*epsilon.yz[0] );

					break;


				// Row 3, column 3

				case 3 :

					result = parameters.E_lh + parameters.av*( epsilon.xx[0] + epsilon.yy[0] + epsilon.zz[0] ) - 0.5*parameters.b*( epsilon.xx[0] + epsilon.yy[0] - 2.0*epsilon.zz[0] ) + phi_piezo + phi_F;

					break;


				// Row 3, column 4

				case 4 :

					result = sqrt( 0.5 )*parameters.b*( epsilon.xx[0] + epsilon.yy[0] - 2.0*epsilon.zz[0] );

					break;


				// Row 3, column 5

				case 5 :

					result = 0.0 + I*0.0;

					break;


				// Row 3, column 6

				case 6 :

					result = -0.5*sqrt( 3.0 )*parameters.b*( epsilon.xx[0] - epsilon.yy[0] ) + I*parameters.d*epsilon.xy[0];

					break;


				// Row 3, column 7

				case 7 :

					result = 0.0 + I*0.0;

					break;


				// Row 3, column 8

				case 8 :

					result = -1.0*sqrt( 1.5 )*parameters.d*( epsilon.xz[0] - I*epsilon.yz[0] );

					break;

			}

			break;


		// Row 4

		case 4 :

			switch( idx_col + 1 )
			{

				// Row 4, column 1

				case 1 :

					result = 0.0 + I*0.0;

					break;


				// Row 4, column 2

				case 2 :

					result = sqrt( 0.5 )*parameters.d*( epsilon.xz[0] + I*epsilon.yz[0] );

					break;


				// Row 4, column 3

				case 3 :

					result = sqrt( 0.5 )*parameters.b*( epsilon.xx[0] + epsilon.yy[0] - 2.0*epsilon.zz[0] );

					break;


				// Row 4, column 4

				case 4 :

					result = parameters.E_so + parameters.av*( epsilon.xx[0] + epsilon.yy[0] + epsilon.zz[0] ) + phi_piezo + phi_F;

					break;


				// Row 4, column 5

				case 5 :

					result = 0.0 + I*0.0;

					break;


				// Row 4, column 6

				case 6 :

					result = -0.5*sqrt( 6.0 )*parameters.b*( epsilon.xx[0] - epsilon.yy[0] ) + I*sqrt( 2.0 )*parameters.d*epsilon.xy[0];

					break;


				// Row 4, column 7

				case 7 :

					result = sqrt( 1.5 )*parameters.d*( epsilon.xz[0] - I*epsilon.yz[0] );

					break;


				// Row 4, column 8

				case 8 :

					result = 0.0 + I*0.0;

					break;

			}

			break;


		// Row 5

		case 5 :

			switch( idx_col + 1 )
			{

				// Row 5, column 1

				case 1 :

					result = 0.0 + I*0.0;

					break;


				// Row 5, column 2

				case 2 :

					result = 0.0 + I*0.0;

					break;


				// Row 5, column 3

				case 3 :

					result = 0.0 + I*0.0;

					break;


				// Row 5, column 4

				case 4 :

					result = 0.0 + I*0.0;

					break;


				// Row 5, column 5

				case 5 :

					result = parameters.E_cb + parameters.ac*( epsilon.xx[0] + epsilon.yy[0] + epsilon.zz[0] ) + phi_piezo + phi_F;

					break;


				// Row 5, column 6

				case 6 :

					result = 0.0 + I*0.0;

					break;


				// Row 5, column 7

				case 7 :

					result = 0.0 + I*0.0;

					break;


				// Row 5, column 8

				case 8 :

					result = 0.0 + I*0.0;

					break;

			}

			break;


		// Row 6

		case 6 :

			switch( idx_col + 1 )
			{

				// Row 6, column 1

				case 1 :

					result = 0.0 + I*0.0;

					break;


				// Row 6, column 2

				case 2 :

					result = 0.0 + I*0.0;

					break;


				// Row 6, column 3

				case 3 :

					result = -0.5*sqrt( 3.0 )*parameters.b*( epsilon.xx[0] - epsilon.yy[0] ) - I*parameters.d*epsilon.xy[0];

					break;


				// Row 6, column 4

				case 4 :

					result = -0.5*sqrt( 6.0 )*parameters.b*( epsilon.xx[0] - epsilon.yy[0] ) - I*sqrt( 2.0 )*parameters.d*epsilon.xy[0];

					break;


				// Row 6, column 5

				case 5 :

					result = 0.0 + I*0.0;

					break;


				// Row 6, column 6

				case 6 :

					result = parameters.E_hh + parameters.av*( epsilon.xx[0] + epsilon.yy[0] + epsilon.zz[0] ) + 0.5*parameters.b*( epsilon.xx[0] + epsilon.yy[0] - 2.0*epsilon.zz[0] ) + phi_piezo + phi_F;

					break;


				// Row 6, column 7

				case 7 :

					result = -1.0*parameters.d*( epsilon.xz[0] + I*epsilon.yz[0] );

					break;


				// Row 6, column 8

				case 8 :

					result = sqrt( 0.5 )*parameters.d*( epsilon.xz[0] + I*epsilon.yz[0] );

					break;

			}

			break;


		// Row 7

		case 7 :

			switch( idx_col + 1 )
			{

				// Row 7, column 1

				case 1 :

					result = 0.0 + I*0.0;

					break;


				// Row 7, column 2

				case 2 :

					result = 0.5*sqrt( 3.0 )*parameters.b*( epsilon.xx[0] - epsilon.yy[0] ) + I*parameters.d*epsilon.xy[0];

					break;


				// Row 7, column 3

				case 3 :

					result = 0.0 + I*0.0;

					break;


				// Row 7, column 4

				case 4 :

					result = sqrt( 1.5 )*parameters.d*( epsilon.xz[0] + I*epsilon.yz[0] );

					break;


				// Row 7, column 5

				case 5 :

					result = 0.0 + I*0.0;

					break;


				// Row 7, column 6

				case 6 :

					result = -1.0*parameters.d*( epsilon.xz[0] - I*epsilon.yz[0] );

					break;


				// Row 7, column 7

				case 7 :

					result = parameters.E_lh + parameters.av*( epsilon.xx[0] + epsilon.yy[0] + epsilon.zz[0] ) - 0.5*parameters.b*( epsilon.xx[0] + epsilon.yy[0] - 2.0*epsilon.zz[0] ) + phi_piezo + phi_F;

					break;


				// Row 7, column 8

				case 8 :

					result = sqrt( 0.5 )*parameters.b*( epsilon.xx[0] + epsilon.yy[0] - 2.0*epsilon.zz[0] );

					break;

			}

			break;


		// Row 8

		case 8 :

			switch( idx_col + 1 )
			{

				// Row 8, column 1

				case 1 :

					result = 0.0 + I*0.0;

					break;


				// Row 8, column 2

				case 2 :

					result = 0.5*sqrt( 6.0 )*parameters.b*( epsilon.xx[0] - epsilon.yy[0] ) + I*sqrt( 2.0 )*parameters.d*epsilon.xy[0];

					break;


				// Row 8, column 3

				case 3 :

					result = -1.0*sqrt( 1.5 )*parameters.d*( epsilon.xz[0] + I*epsilon.yz[0] );

					break;


				// Row 8, column 4

				case 4 :

					result = 0.0 + I*0.0;

					break;


				// Row 8, column 5

				case 5 :

					result = 0.0 + I*0.0;

					break;


				// Row 8, column 6

				case 6 :

					result = 1.0*sqrt( 0.5 )*parameters.d*( epsilon.xz[0] - I*epsilon.yz[0] );

					break;


				// Row 8, column 7

				case 7 :

					result = sqrt( 0.5 )*parameters.b*( epsilon.xx[0] + epsilon.yy[0] - 2.0*epsilon.zz[0] );

					break;


				// Row 8, column 8

				case 8 :

					result = parameters.E_so + parameters.av*( epsilon.xx[0] + epsilon.yy[0] + epsilon.zz[0] ) + phi_piezo + phi_F;

					break;

			}

			break;

	}


	return result;

}
