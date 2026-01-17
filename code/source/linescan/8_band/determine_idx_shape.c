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
// Note: This function implements the real space characterstic function in the case of "hard" shape interfaces.
//       The case of Gaussian interface softening is handled in characteristic_function_real_space.c, which computes
//       the inverse Fourier transform of the reciprocal space characterstic function numerically.


#include "8_band_linescan.h"

int determine_idx_shape( int N_shape, double x0, double y0, double z0, struct supercell_geometry *supercell )
{

	int idx_shape, result = 0;

	double min_vol = 0, shape_vol = 0;

	// Loop over shapes
	for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
	{

		// Accounting for shape rotation for simplicity
		double sin_alpha = sin( M_PI*supercell[idx_shape].alpha/180 );
		double cos_alpha = cos( M_PI*supercell[idx_shape].alpha/180 );
		double sin_beta = sin( M_PI*supercell[idx_shape].beta/180 );
		double cos_beta = cos( M_PI*supercell[idx_shape].beta/180 );
		double sin_gamma = sin( M_PI*supercell[idx_shape].gamma/180 );
		double cos_gamma = cos( M_PI*supercell[idx_shape].gamma/180 );
		
		double x = ( ( cos_beta * cos_alpha * cos_gamma ) - ( sin_alpha * sin_gamma ))*x0 + ( ( cos_gamma * cos_beta * sin_alpha ) + ( sin_gamma * cos_alpha ) )*y0 + ( -sin_beta * cos_gamma )*z0;
	
		double y = ( ( cos_beta * -cos_alpha * sin_gamma ) - ( sin_alpha * cos_gamma ) )*x0 + ( ( -sin_gamma * cos_beta * sin_alpha ) + ( cos_gamma * cos_alpha ) )*y0 + ( sin_beta * sin_gamma )*z0;
	
		double z = ( sin_beta * cos_alpha )*x0 + ( sin_alpha * sin_beta )*y0 + ( cos_beta )*z0;

		switch( supercell[idx_shape].shape )
		{

			////////////////
			//// Cuboid ////
			////////////////

			case 0:
			{

				if(    ( x >= supercell[idx_shape].Ox - 0.5*supercell[idx_shape].L1 ) && ( x <= supercell[idx_shape].Ox + 0.5*supercell[idx_shape].L1 )
			 	    && ( y >= supercell[idx_shape].Oy - 0.5*supercell[idx_shape].L2 ) && ( y <= supercell[idx_shape].Oy + 0.5*supercell[idx_shape].L2 )
    			 	    && ( z >= supercell[idx_shape].Oz - 0.5*supercell[idx_shape].L3 ) && ( z <= supercell[idx_shape].Oz + 0.5*supercell[idx_shape].L3 ) )
				{

					shape_vol = supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;
			}


			////////////////
			//// Sphere ////
			////////////////

			case 1:
			{
				if( sqrt( pow( x - supercell[idx_shape].Ox, 2.0 )
					   + pow( y - supercell[idx_shape].Oy, 2.0 )
					   + pow( z - supercell[idx_shape].Oz, 2.0 ) ) <= supercell[idx_shape].L1 )
				{

					shape_vol = ( 4.0*M_PI*pow( supercell[idx_shape].L1, 3.0 ) )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			///////////////////
			//// Ellipsoid ////
			///////////////////

			case 2:
			{

				if( sqrt( pow( ( x - supercell[idx_shape].Ox )/supercell[idx_shape].L1, 2.0 )
					   + pow( ( y - supercell[idx_shape].Oy )/supercell[idx_shape].L2, 2.0 )
					   + pow( ( z - supercell[idx_shape].Oz )/supercell[idx_shape].L3, 2.0 ) ) <= 1.0 )
				{

					shape_vol = ( 4.0*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3 )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			////////////////////
			//// Hemisphere ////
			////////////////////

			case 3:
			{

				if(    ( sqrt( pow( x - supercell[idx_shape].Ox, 2.0 )
						   + pow( y - supercell[idx_shape].Oy, 2.0 )
						   + pow( z - supercell[idx_shape].Oz, 2.0 ) ) <= supercell[idx_shape].L1 )
				    && ( z >= supercell[idx_shape].Oz ) )
				{

					shape_vol = ( 2.0*M_PI*pow( supercell[idx_shape].L1, 3.0 ) )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			///////////////////////
			//// Hemiellipsoid ////
			///////////////////////

			case 4:
			{

				if(    ( sqrt( pow( ( x - supercell[idx_shape].Ox )/supercell[idx_shape].L1, 2.0 )
						   + pow( ( y - supercell[idx_shape].Oy )/supercell[idx_shape].L2, 2.0 )
						   + pow( ( z - supercell[idx_shape].Oz )/supercell[idx_shape].L3, 2.0 ) ) <= 1.0 )
				    && ( z >= supercell[idx_shape].Oz ) )
				{

					shape_vol = ( 2.0*M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3 )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			/////////////////
			//// Pyramid ////
			/////////////////

			case 5:
			{

				if(    ( x >= supercell[idx_shape].Ox - supercell[idx_shape].L1*( supercell[idx_shape].L3 - ( z - supercell[idx_shape].Oz ) )/supercell[idx_shape].L3 )
				    && ( x <= supercell[idx_shape].Ox + supercell[idx_shape].L1*( supercell[idx_shape].L3 - ( z - supercell[idx_shape].Oz ) )/supercell[idx_shape].L3 )
			 	    && ( y >= supercell[idx_shape].Oy - supercell[idx_shape].L2*( supercell[idx_shape].L3 - ( z - supercell[idx_shape].Oz ) )/supercell[idx_shape].L3 )
				    && ( y <= supercell[idx_shape].Oy + supercell[idx_shape].L2*( supercell[idx_shape].L3 - ( z - supercell[idx_shape].Oz ) )/supercell[idx_shape].L3 )
    			 	    && ( z >= supercell[idx_shape].Oz )
				    && ( z <= supercell[idx_shape].Oz + supercell[idx_shape].L3 ) )
				{

					shape_vol = ( 4.0*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3 )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			///////////////////////////
			//// Truncated pyramid ////
			///////////////////////////

			case 6:
			{

				double H = ( supercell[idx_shape].L4*supercell[idx_shape].L1 )/( supercell[idx_shape].L1 - supercell[idx_shape].L3 );

				if(    ( x >= supercell[idx_shape].Ox - supercell[idx_shape].L1*( H - ( z - supercell[idx_shape].Oz ) )/H )
				    && ( x <= supercell[idx_shape].Ox + supercell[idx_shape].L1*( H - ( z - supercell[idx_shape].Oz ) )/H )
			 	    && ( y >= supercell[idx_shape].Oy - supercell[idx_shape].L2*( H - ( z - supercell[idx_shape].Oz ) )/H )
				    && ( y <= supercell[idx_shape].Oy + supercell[idx_shape].L2*( H - ( z - supercell[idx_shape].Oz ) )/H )
    			 	    && ( z >= supercell[idx_shape].Oz )
				    && ( z <= supercell[idx_shape].Oz + supercell[idx_shape].L4 ) )
				{

					shape_vol = ( 4.0*( supercell[idx_shape].L1*supercell[idx_shape].L2*H - supercell[idx_shape].L3*supercell[idx_shape].L2*supercell[idx_shape].L3/supercell[idx_shape].L1*( H - supercell[idx_shape].L4 ) ) )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			/////////////////////////
			//// N-sided pyramid ////
			/////////////////////////

			case 7:
			{

				int segment, N = lround( 0.1*supercell[idx_shape].L3 );
				double theta;

				if( ( y - supercell[idx_shape].Oy >= 0.0 ) )
				{

					theta = atan2( ( y - supercell[idx_shape].Oy ), ( x - supercell[idx_shape].Ox ) );

				}

				else
				{

					theta = atan2( ( y - supercell[idx_shape].Oy ), ( x - supercell[idx_shape].Ox ) ) + 2*M_PI;

				}


				segment = ( int ) ( ( ( ( double ) N )*( theta + M_PI/( ( double ) N ) ) )/( 2.0*M_PI ) );
				segment = segment % N;

				double x1 =      cos( ( 2.0*M_PI*segment )/( ( double ) N ) )*x + sin( ( 2.0*M_PI*segment )/( ( double ) N ) )*y;
				double y1 = -1.0*sin( ( 2.0*M_PI*segment )/( ( double ) N ) )*x + cos( ( 2.0*M_PI*segment )/( ( double ) N ) )*y;


				if(    ( x1 >= supercell[idx_shape].Ox )
				    && ( x1 <= supercell[idx_shape].Ox + supercell[idx_shape].L1*( supercell[idx_shape].L2 - ( z - supercell[idx_shape].Oz ) ) )
			 	    && ( y1 >= supercell[idx_shape].Oy - x1*tan( M_PI/( ( double ) N ) ) )
				    && ( y1 <= supercell[idx_shape].Oy + x1*tan( M_PI/( ( double ) N ) ) )
    			 	    && ( z  >= supercell[idx_shape].Oz )
				    && ( z  <= supercell[idx_shape].Oz + supercell[idx_shape].L2 ) )
				{

					double alpha = M_PI/( ( double ) N );
					double x0    = supercell[idx_shape].L1;
					double y0    = supercell[idx_shape].L1*tan( alpha );

					shape_vol = ( ( ( double ) N )*x0*y0*supercell[idx_shape].L2 )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			///////////////////////////////////
			//// Truncated N-sided pyramid ////
			///////////////////////////////////

			case 8:
			{

				int segment, N = lround( 0.1*supercell[idx_shape].L4 ); // Assuming a 100-sided truncated pyramid (approximates a truncated circular-base cone)
				double theta, H;

				if( ( y - supercell[idx_shape].Oy >= 0.0 ) )
				{

					theta = atan2( ( y - supercell[idx_shape].Oy ), ( x - supercell[idx_shape].Ox ) );

				}

				else
				{

					theta = atan2( ( y - supercell[idx_shape].Oy ), ( x - supercell[idx_shape].Ox ) ) + 2*M_PI;

				}


				H = ( supercell[idx_shape].L3*supercell[idx_shape].L1 )/( supercell[idx_shape].L1 - supercell[idx_shape].L2 );

				segment = ( int ) ( ( ( ( double ) N )*( theta + M_PI/( ( double ) N ) ) )/( 2.0*M_PI ) );
				segment = segment % N;

				double x1 =      cos( ( 2.0*M_PI*segment )/( ( double ) N ) )*x + sin( ( 2.0*M_PI*segment )/( ( double ) N ) )*y;
				double y1 = -1.0*sin( ( 2.0*M_PI*segment )/( ( double ) N ) )*x + cos( ( 2.0*M_PI*segment )/( ( double ) N ) )*y;


				if(    ( x1 >= supercell[idx_shape].Ox )
				    && ( x1 <= supercell[idx_shape].Ox + supercell[idx_shape].L1*( H - ( z - supercell[idx_shape].Oz ) )/H )
			 	    && ( y1 >= supercell[idx_shape].Oy - x1*tan( M_PI/( ( double ) N ) ) )
				    && ( y1 <= supercell[idx_shape].Oy + x1*tan( M_PI/( ( double ) N ) ) )
    			 	    && ( z  >= supercell[idx_shape].Oz )
				    && ( z  <= supercell[idx_shape].Oz + supercell[idx_shape].L3 ) )
				{

					double alpha = M_PI/( ( double ) N );

					double x0 = supercell[idx_shape].L1;
					double x1 = supercell[idx_shape].L2;

					double y0 = supercell[idx_shape].L1*tan( alpha );
					double y1 = supercell[idx_shape].L2*tan( alpha );

					shape_vol = ( ( ( double ) N )*( x0*y0*H - x1*y1*( H - supercell[idx_shape].L3 ) ) )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			//////////////////
			//// Cylinder ////
			//////////////////

			case 9:
			{

				if(    ( sqrt( pow( x - supercell[idx_shape].Ox, 2.0 ) + pow( y - supercell[idx_shape].Oy, 2.0 ) ) <= supercell[idx_shape].L1 )
				    && ( z >= supercell[idx_shape].Oz - 0.5*supercell[idx_shape].L2 )
				    && ( z <= supercell[idx_shape].Oz + 0.5*supercell[idx_shape].L2 ) )
				{

					shape_vol = M_PI*pow( supercell[idx_shape].L1, 2.0 )*supercell[idx_shape].L2;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			//////////////////////////
			//// Cylindrical ring ////
			//////////////////////////

			case 10:
			{

				if(    ( sqrt( pow( x - supercell[idx_shape].Ox, 2.0 ) + pow( y - supercell[idx_shape].Oy, 2.0 ) ) <= supercell[idx_shape].L1 )
				    && ( sqrt( pow( x - supercell[idx_shape].Ox, 2.0 ) + pow( y - supercell[idx_shape].Oy, 2.0 ) ) >= supercell[idx_shape].L2 )
				    && ( z >= supercell[idx_shape].Oz - 0.5*supercell[idx_shape].L3 )
				    && ( z <= supercell[idx_shape].Oz + 0.5*supercell[idx_shape].L3 ) )
				{

					shape_vol = M_PI*( pow( supercell[idx_shape].L1, 2.0 ) - pow( supercell[idx_shape].L2, 2.0 ) )*supercell[idx_shape].L3;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			//////////////
			//// Cone ////
			//////////////

			case 11:
			{

				if(    ( sqrt( pow( x - supercell[idx_shape].Ox, 2.0 ) + pow( y - supercell[idx_shape].Oy, 2.0 ) ) <= ( supercell[idx_shape].L1*( supercell[idx_shape].L2 - ( z - supercell[idx_shape].Oz ) )/supercell[idx_shape].L2 ) )
				    && ( z >= supercell[idx_shape].Oz )
				    && ( z <= supercell[idx_shape].Oz + supercell[idx_shape].L2 ) )
				{

					shape_vol = ( M_PI*pow( supercell[idx_shape].L1, 2.0 )*supercell[idx_shape].L2 )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			////////////////////////
			//// Truncated cone ////
			////////////////////////

			case 12:
			{

				double H = ( supercell[idx_shape].L3*supercell[idx_shape].L1 )/( supercell[idx_shape].L1 - supercell[idx_shape].L2 );

				if(    ( sqrt( pow( x - supercell[idx_shape].Ox, 2.0 ) + pow( y - supercell[idx_shape].Oy, 2.0 ) ) <= ( supercell[idx_shape].L1*( H - ( z - supercell[idx_shape].Oz ) )/H ) )
				    && ( z >= supercell[idx_shape].Oz )
				    && ( z <= supercell[idx_shape].Oz + supercell[idx_shape].L3 ) )
				{

					shape_vol = M_PI*( pow( supercell[idx_shape].L1, 2.0 )*H - pow( supercell[idx_shape].L2, 2.0 )*( H - supercell[idx_shape].L3 ) );

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			//////////////////////////////
			//// Truncated conic ring ////
			//////////////////////////////

			case 13:
			{

				double H = ( supercell[idx_shape].L4*supercell[idx_shape].L1 )/( supercell[idx_shape].L1 - supercell[idx_shape].L2 );

				if(    ( sqrt( pow( x - supercell[idx_shape].Ox, 2.0 ) + pow( y - supercell[idx_shape].Oy, 2.0 ) ) <= ( supercell[idx_shape].L1*( H - ( z - supercell[idx_shape].Oz ) )/H ) )
				    && ( sqrt( pow( x - supercell[idx_shape].Ox, 2.0 ) + pow( y - supercell[idx_shape].Oy, 2.0 ) ) >= ( supercell[idx_shape].L3 ) )
				    && ( z >= supercell[idx_shape].Oz )
				    && ( z <= supercell[idx_shape].Oz + supercell[idx_shape].L4 ) )
				{

					shape_vol = M_PI*( pow( supercell[idx_shape].L1, 2.0 )*H - pow( supercell[idx_shape].L2, 2.0 )*( H - supercell[idx_shape].L4 ) );

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			/////////////////////////
			//// Elliptical cone ////
			/////////////////////////

			case 14:
			{
				if(    ( sqrt( pow( ( x - supercell[idx_shape].Ox )/supercell[idx_shape].L1, 2.0 ) + pow( ( y - supercell[idx_shape].Oy )/supercell[idx_shape].L2, 2.0 ) ) <= ( 1.0 - ( z - supercell[idx_shape].Oz )/supercell[idx_shape].L3 ) )
				    && ( z >= supercell[idx_shape].Oz )
				    && ( z <= supercell[idx_shape].Oz + supercell[idx_shape].L3 ) )
				{

					shape_vol = ( M_PI*supercell[idx_shape].L1*supercell[idx_shape].L2*supercell[idx_shape].L3 )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			///////////////////////////////////
			//// Truncated elliptical cone ////
			///////////////////////////////////

			case 15:
			{
				double H = ( supercell[idx_shape].L4*supercell[idx_shape].L1 )/( supercell[idx_shape].L1 - supercell[idx_shape].L3 );

				if(    ( sqrt( pow( ( x - supercell[idx_shape].Ox )/supercell[idx_shape].L1, 2.0 ) + pow( ( y - supercell[idx_shape].Oy )/supercell[idx_shape].L2, 2.0 ) ) <= ( 1.0 - ( z - supercell[idx_shape].Oz )/H ) )
				    && ( z >= supercell[idx_shape].Oz )
				    && ( z <= supercell[idx_shape].Oz + supercell[idx_shape].L4 ) )
				{

					shape_vol = ( M_PI*( supercell[idx_shape].L1*supercell[idx_shape].L2*H - supercell[idx_shape].L3*supercell[idx_shape].L2*supercell[idx_shape].L3/supercell[idx_shape].L1*( H - supercell[idx_shape].L4 ) ) )/3.0;

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;

			}


			///////////////////////
			//// N-sided prism ////
			///////////////////////

			case 16:
			{

				int segment, N = lround( 0.1*supercell[idx_shape].L3 );
				double theta;

				if( ( y - supercell[idx_shape].Oy >= 0.0 ) )
				{

					theta = atan2( ( y - supercell[idx_shape].Oy ), ( x - supercell[idx_shape].Ox ) );

				}

				else
				{

					theta = atan2( ( y - supercell[idx_shape].Oy ), ( x - supercell[idx_shape].Ox ) ) + 2*M_PI;

				}


				segment = ( int ) ( ( ( ( double ) N )*( theta + M_PI/( ( double ) N ) ) )/( 2.0*M_PI ) );
				segment = segment % N;

				double x1 =      cos( ( 2.0*M_PI*segment )/( ( double ) N ) )*x + sin( ( 2.0*M_PI*segment )/( ( double ) N ) )*y;
				double y1 = -1.0*sin( ( 2.0*M_PI*segment )/( ( double ) N ) )*x + cos( ( 2.0*M_PI*segment )/( ( double ) N ) )*y;


				if(    ( x1 >= supercell[idx_shape].Ox )
				    && ( x1 <= supercell[idx_shape].Ox + supercell[idx_shape].L1 )
			 	    && ( y1 >= supercell[idx_shape].Oy - x1*tan( M_PI/( ( double ) N ) ) )
				    && ( y1 <= supercell[idx_shape].Oy + x1*tan( M_PI/( ( double ) N ) ) )
    			 	    && ( z  >= supercell[idx_shape].Oz - supercell[idx_shape].L2/2 )
				    && ( z  <= supercell[idx_shape].Oz + supercell[idx_shape].L2/2 ) )
				{

					double alpha = M_PI/( ( double ) N );
					double x0    = supercell[idx_shape].L1;
					double y0    = supercell[idx_shape].L1*tan( alpha );

					shape_vol = ( ( ( double ) N )*x0*y0*supercell[idx_shape].L2 );

					if( ( min_vol == 0 ) || ( shape_vol < min_vol ) )
					{

						min_vol = shape_vol;
						 result = idx_shape;

					}

				}

				break;
				
			}


			//////////////////////
			//// Quantum well ////
			//////////////////////
			
			case 17:
			{

				if ( lround( 0.1*supercell[idx_shape].L1 ) == 0 )
				{

					if ( fabs( x - supercell[idx_shape].Ox ) < supercell[idx_shape].L2/2 )
					{

						shape_vol = supercell[idx_shape].L2*supercell[0].L2*supercell[0].L3;
						result = idx_shape;

					}

				}
				else if ( lround( 0.1*supercell[idx_shape].L1 ) == 1 )
				{

					if ( fabs( y - supercell[idx_shape].Oy ) < supercell[idx_shape].L2/2 )
					{

						shape_vol = supercell[idx_shape].L2*supercell[0].L1*supercell[0].L3;
						result = idx_shape;

					}

				}
				else if ( lround( 0.1*supercell[idx_shape].L1 ) == 2 )
				{

					if ( fabs( z - supercell[idx_shape].Oz ) < supercell[idx_shape].L2/2 )
					{

						shape_vol = supercell[idx_shape].L2*supercell[0].L1*supercell[0].L2;
						result = idx_shape;

					}

				}

				break;

			}


			//////////////////////
			//// Quantum wire ////
			//////////////////////
			
			case 18:
			{

				if ( lround( 0.1*supercell[idx_shape].L1 ) == 0 )
				{

					if ( fabs( y - supercell[idx_shape].Oy ) < supercell[idx_shape].L2/2 && fabs( z - supercell[idx_shape].Oz ) < supercell[idx_shape].L3/2 )
					{

						shape_vol = supercell[idx_shape].L2*supercell[idx_shape].L3*supercell[0].L1;
						result = idx_shape;

					}

				}
				else if ( lround( 0.1*supercell[idx_shape].L1 ) == 1 )
				{

					if ( fabs( x - supercell[idx_shape].Ox ) < supercell[idx_shape].L2/2 && fabs( z - supercell[idx_shape].Oz ) < supercell[idx_shape].L3/2 )
					{

						shape_vol = supercell[idx_shape].L2*supercell[idx_shape].L3*supercell[0].L2;
						result = idx_shape;

					}

				}
				else if ( lround( 0.1*supercell[idx_shape].L1 ) == 2 )
				{

					if ( fabs( x - supercell[idx_shape].Ox ) < supercell[idx_shape].L2/2 && fabs( y - supercell[idx_shape].Oy ) < supercell[idx_shape].L3/2 )
					{

						shape_vol = supercell[idx_shape].L2*supercell[idx_shape].L3*supercell[0].L3;
						result = idx_shape;

					}

				}

				break;

			}

		}

	}


	return result;

}
