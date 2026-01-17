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

void strain_region_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct strain_tensor strain_G, struct characteristic_functions *characteristic_function_G )
{

	// Isotropic approximation to elastic tensor: C11 - C12 - 2*C44 = 0

	if( isotropic != 0 )
	{

		for( int idx_shape = 0; idx_shape < N_shape; idx_shape++ )
		{

			parameters[idx_shape].C44 = 0.5*( parameters[idx_shape].C11 - parameters[idx_shape].C12 ); // C44 = ( C11 - C12 )/2

			// parameters[idx_shape].C11 = parameters[idx_shape].C12 + 2.0*parameters[idx_shape].C44; // C11 = C12 + 2*C44
			// parameters[idx_shape].C12 = parameters[idx_shape].C11 - 2.0*parameters[idx_shape].C44; // C12 = C11 - 2*C44

		}

	}


	// In case of position-independent elastic constants
	if ( include_shape_elastic_constants == 0 )
	{

		struct convolution_char_strain convolutions_char_strain_correction_dummy;

		#pragma omp parallel
		{

			// Loop over wave vectors and evaluate the Fourier coefficients of the chosen component of the strain tensor
			#pragma omp for collapse( 3 )
			for( int idx_Gx = 0; idx_Gx < basis.N_Gx_T; idx_Gx++ )
			{

				for( int idx_Gy = 0; idx_Gy < basis.N_Gy_T; idx_Gy++ )
				{

					for( int idx_Gz = 0; idx_Gz < basis.N_Gz_T; idx_Gz++ )
					{

						int idx_basis = idx_Gz + basis.N_Gz_T*idx_Gy + basis.N_Gy_T*basis.N_Gz_T*idx_Gx;

						strain_G.xx[idx_basis] = 0.0 + I*0.0;
						strain_G.yy[idx_basis] = 0.0 + I*0.0;
						strain_G.zz[idx_basis] = 0.0 + I*0.0;
						strain_G.yz[idx_basis] = 0.0 + I*0.0;
						strain_G.xz[idx_basis] = 0.0 + I*0.0;
						strain_G.xy[idx_basis] = 0.0 + I*0.0;

						for( int component = 0; component < 6; component++ )
						{

							// Loop over shapes
							for( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
							{

								switch( component )
								{

									// Note: Since we are calculating the strain for shapes (materials) embedded in a defined matrix material, the supercell matrix material
									//	    (idx_shape = 0) is not included in sums over the shapes in the supercell - i.e. idx_shape >= 1 only in the strain calculation

									case 0: // component = 0 --> xx

										strain_G.xx[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], 0, &convolutions_char_strain_correction_dummy );

										break;


									case 1: // component = 1 --> yy

										strain_G.yy[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], 0, &convolutions_char_strain_correction_dummy );

										break;


									case 2: // component = 2 --> zz

										strain_G.zz[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], 0, &convolutions_char_strain_correction_dummy );

										break;


									case 3: // component = 3 --> yz

										strain_G.yz[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], 0, &convolutions_char_strain_correction_dummy );

										break;


									case 4: // component = 4 --> xz

										strain_G.xz[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], 0, &convolutions_char_strain_correction_dummy );

										break;


									case 5: // component = 5 --> xy

										strain_G.xy[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], 0, &convolutions_char_strain_correction_dummy );

										break;


								}

							}

						}

					}

				}

			}

		}

	}
	else // In case of position-dependent elastic constants
	{

		// Will store updated strain tensor after each order correction
		struct strain_tensor strain_correction_region_G;

		strain_correction_region_G.xx = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_correction_region_G.yy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_correction_region_G.zz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_correction_region_G.yz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_correction_region_G.xz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_correction_region_G.xy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );

		// Will store convolution of shape characteristic function with corrected strain tensor after each order correction, for calculation of next order correction
		struct convolution_char_strain convolutions_char_strain_correction;

		convolutions_char_strain_correction.char_xx = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		convolutions_char_strain_correction.char_yy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		convolutions_char_strain_correction.char_zz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		convolutions_char_strain_correction.char_yz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		convolutions_char_strain_correction.char_xz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		convolutions_char_strain_correction.char_xy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );

		// Loop over shapes
		for ( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			// Loop over correction orders
			for ( int correction_order = -1; correction_order <= shape_elastic_constants_order; correction_order++ )
			{

				// Calculates convolution of previous strain tensor correction with shape characteristic function
				if ( correction_order > 0 )
				{
					for ( int idx_basis = 0; idx_basis < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T; idx_basis++ )
					{

						convolutions_char_strain_correction.char_xx[idx_basis] = 0;
						convolutions_char_strain_correction.char_yy[idx_basis] = 0;
						convolutions_char_strain_correction.char_zz[idx_basis] = 0;
						convolutions_char_strain_correction.char_yz[idx_basis] = 0;
						convolutions_char_strain_correction.char_xz[idx_basis] = 0;
						convolutions_char_strain_correction.char_xy[idx_basis] = 0;

					}

					convolutions_strain_correction( idx_shape, basis, characteristic_function_G, &strain_correction_region_G, &convolutions_char_strain_correction );

				}

				
				#pragma omp parallel
				{

					// Loop over local wave vectors
					#pragma omp for collapse( 3 )
					for( int idx_Gx = 0; idx_Gx < basis.N_Gx_T; idx_Gx++ )
					{

						for( int idx_Gy = 0; idx_Gy < basis.N_Gy_T; idx_Gy++ )
						{

							for( int idx_Gz = 0; idx_Gz < basis.N_Gz_T; idx_Gz++ )
							{

								int idx_basis = idx_Gz + basis.N_Gz_T*idx_Gy + basis.N_Gy_T*basis.N_Gz_T*idx_Gx;

								for( int component = 0; component < 6; component++ )
								{

									// Note: Sum of shape strain contributions is stored in strain_G[0]

									if (component == 0) // component = 0 --> xx
									{
									
										if ( correction_order == - 1 ) // Initial contribution
										{

											strain_G.xx[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );

										}
										else // Subsequent order corrections
										{
										
											// Calculating correction from order
											strain_correction_region_G.xx[idx_basis] = strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );
											// Adding contribution from order
											strain_G.xx[idx_basis] += strain_correction_region_G.xx[idx_basis];

										}
									
									}

									else if (component == 1) // component = 1 --> yy
									{
									
										if ( correction_order == - 1 ) // Initial contribution
										{

											strain_G.yy[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );

										}
										else // Subsequent order corrections
										{
										
											// Calculating correction from order
											strain_correction_region_G.yy[idx_basis] = strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );
											// Adding contribution from order
											strain_G.yy[idx_basis] += strain_correction_region_G.yy[idx_basis];

										}

									}


									else if (component == 2) // component = 2 --> zz
									{

										if ( correction_order == - 1 ) // Initial contribution
										{

											strain_G.zz[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );

										}
										else // Subsequent order corrections
										{
										
											// Calculating correction from order
											strain_correction_region_G.zz[idx_basis] = strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );
											// Adding contribution from order
											strain_G.zz[idx_basis] += strain_correction_region_G.zz[idx_basis];

										}
								
									}


									else if (component == 3) // component = 3 --> yz
									{
									
										if ( correction_order == - 1 ) // Initial contribution
										{

											strain_G.yz[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );

										}
										else // Subsequent order corrections
										{
										
											// Calculating correction from order
											strain_correction_region_G.yz[idx_basis] = strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );
											// Adding contribution from order
											strain_G.yz[idx_basis] += strain_correction_region_G.yz[idx_basis];

										}
										
									}


									else if (component == 4) // component = 4 --> xz
									{

										if ( correction_order == - 1 ) // Initial contribution
										{

											strain_G.xz[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );

										}
										else // Subsequent order corrections
										{
										
											strain_correction_region_G.xz[idx_basis] = strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );
											// Adding contribution from order
											strain_G.xz[idx_basis] += strain_correction_region_G.xz[idx_basis];

										}

									}


									else // component = 5 --> xy
									{
									
										if ( correction_order == - 1 ) // Initial contribution
										{

											strain_G.xy[idx_basis] += strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );

										}
										else // Subsequent order corrections
										{
										
											// Calculating correction from order
											strain_correction_region_G.xy[idx_basis] = strain_fourier_coefficient( component, isotropic, idx_shape, idx_basis, 0, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, characteristic_function_G[idx_shape], correction_order, &convolutions_char_strain_correction );
											// Adding contribution from order
											strain_G.xy[idx_basis] += strain_correction_region_G.xy[idx_basis];

										}

									}

								}

							}

						}

					}	

				}
			
			}

		}

		free( strain_correction_region_G.xx );
		free( strain_correction_region_G.yy );
		free( strain_correction_region_G.zz );
		free( strain_correction_region_G.yz );
		free( strain_correction_region_G.xz );
		free( strain_correction_region_G.xy );


		free( convolutions_char_strain_correction.char_xx );
		free( convolutions_char_strain_correction.char_yy );
		free( convolutions_char_strain_correction.char_zz );
		free( convolutions_char_strain_correction.char_yz );
		free( convolutions_char_strain_correction.char_xz );
		free( convolutions_char_strain_correction.char_xy );

	}

}
