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

void strain_region_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters )
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
	
	// Local size of arrays
	int local_size = basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs + ( ( myid < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? 1 : 0 );
	
	struct strain_tensor strain_region_G_local;
	
	strain_region_G_local.xx = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	strain_region_G_local.yy = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	strain_region_G_local.zz = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	strain_region_G_local.yz = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	strain_region_G_local.xz = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	strain_region_G_local.xy = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	
	
	// Will store number of strain tensor components calculated by each process
	int *recvcounts = ( int * ) malloc( nprocs*sizeof( int ) );
	
	// Will store displacement from first value of the "global" strain tensor component for the first value of the strain tensor component calculated by each process
	int *displs = ( int * ) malloc( nprocs*sizeof( int ) );
	
	for ( int proc = 0; proc < nprocs; proc++ )
	{
	
		recvcounts[proc] = basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs + ( ( proc < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? 1 : 0 );
		displs[proc] = ( proc == 0 ? 0 : displs[proc - 1] + recvcounts[proc - 1] ); 
	
	}

	// In case of position-independent elastic constants
	if ( include_shape_elastic_constants == 0 )
	{

		struct convolution_char_strain convolutions_char_strain_correction_dummy;

		// Loop over shapes
		for( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{
		
			// Loop over wave vectors to evaluate the Fourier coefficients of the chosen component of the strain tensor
			#pragma omp parallel for
			for( int idx = 0; idx < local_size; idx++ )
			{

				// Identifying global wave vector from local index
				int idx_basis = myid * ( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs ) + ( ( myid < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? myid : basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) + idx;
				
				int idx_Gx = idx_basis / ( basis.N_Gy_T*basis.N_Gz_T );
			
				int idx_Gy = ( idx_basis % ( basis.N_Gy_T*basis.N_Gz_T ) ) / basis.N_Gz_T;
			
				int idx_Gz = ( idx_basis % ( basis.N_Gy_T*basis.N_Gz_T ) ) % basis.N_Gz_T;

				for( int component = 0; component < 6; component++ )
				{

					// Note: Sum of shape strain contributions is stored in strain_G[0]

					if (component == 0) // component = 0 --> xx
					{
					
						strain_region_G_local.xx[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, 0, &convolutions_char_strain_correction_dummy );
					
					}

					else if (component == 1) // component = 1 --> yy
					{
					
						strain_region_G_local.yy[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, 0, &convolutions_char_strain_correction_dummy );

					}


					else if (component == 2) // component = 2 --> zz
					{
					
						strain_region_G_local.zz[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, 0, &convolutions_char_strain_correction_dummy );
				
					}


					else if (component == 3) // component = 3 --> yz
					{
					
						strain_region_G_local.yz[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, 0, &convolutions_char_strain_correction_dummy );
						
					}


					else if (component == 4) // component = 4 --> xz
					{
					
						strain_region_G_local.xz[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, 0, &convolutions_char_strain_correction_dummy );

					}


					else // component = 5 --> xy
					{
					
						strain_region_G_local.xy[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, 0, &convolutions_char_strain_correction_dummy );

					}

				}

			}
		
		}

	}
	else // In case of position-dependent elastic constants
	{

		// Will store updated global strain tensor after each order correction
		struct strain_tensor strain_correction_region_G;

		strain_correction_region_G.xx = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_correction_region_G.yy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_correction_region_G.zz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_correction_region_G.yz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_correction_region_G.xz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_correction_region_G.xy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );


		// Will store global convolution of shape characteristic function with corrected strain tensor after each order correction, for calculation of next order correction
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

				if ( correction_order > 0 ) // Calculates convolution of previous strain tensor correction with shape characteristic function
				{

					convolutions_strain_correction( idx_shape, basis, &strain_correction_region_G, &convolutions_char_strain_correction );

				}

				// Will store local strain correction
				struct strain_tensor strain_correction_region_G_local;

				if ( correction_order > -1 )
				{

					strain_correction_region_G_local.xx = ( complex double* ) calloc( local_size, sizeof( complex double ) );
					strain_correction_region_G_local.yy = ( complex double* ) calloc( local_size, sizeof( complex double ) );
					strain_correction_region_G_local.zz = ( complex double* ) calloc( local_size, sizeof( complex double ) );
					strain_correction_region_G_local.yz = ( complex double* ) calloc( local_size, sizeof( complex double ) );
					strain_correction_region_G_local.xz = ( complex double* ) calloc( local_size, sizeof( complex double ) );
					strain_correction_region_G_local.xy = ( complex double* ) calloc( local_size, sizeof( complex double ) );

				}

				
				// Loop over local wave vectors
				#pragma omp parallel for
				for( int idx = 0; idx < local_size; idx++ )
				{

					// Global wave vector index
					int idx_basis = myid * ( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs ) + ( ( myid < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? myid : basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) + idx;
					
					int idx_Gx = idx_basis / ( basis.N_Gy_T*basis.N_Gz_T );
				
					int idx_Gy = ( idx_basis % ( basis.N_Gy_T*basis.N_Gz_T ) ) / basis.N_Gz_T;
				
					int idx_Gz = ( idx_basis % ( basis.N_Gy_T*basis.N_Gz_T ) ) % basis.N_Gz_T;

					for( int component = 0; component < 6; component++ )
					{

						if (component == 0) // component = 0 --> xx
						{
						
							if ( correction_order == - 1 ) // Initial contribution
							{

								strain_region_G_local.xx[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );

							}
							else // Subsequent order corrections
							{
							
								// Calculating correction from order
								strain_correction_region_G_local.xx[idx] = strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );
								// Adding contribution from order
								strain_region_G_local.xx[idx] += strain_correction_region_G_local.xx[idx];

							}
						
						}

						else if (component == 1) // component = 1 --> yy
						{
						
							if ( correction_order == - 1 ) // Initial contribution
							{

								strain_region_G_local.yy[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );

							}
							else // Subsequent order corrections
							{
							
								// Calculating correction from order
								strain_correction_region_G_local.yy[idx] = strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );
								// Adding contribution from order
								strain_region_G_local.yy[idx] += strain_correction_region_G_local.yy[idx];

							}

						}


						else if (component == 2) // component = 2 --> zz
						{

							if ( correction_order == - 1 ) // Initial contribution
							{

								strain_region_G_local.zz[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );

							}
							else // Subsequent order corrections
							{
							
								// Calculating correction from order
								strain_correction_region_G_local.zz[idx] = strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );
								// Adding contribution from order
								strain_region_G_local.zz[idx] += strain_correction_region_G_local.zz[idx];

							}
					
						}


						else if (component == 3) // component = 3 --> yz
						{
						
							if ( correction_order == - 1 ) // Initial contribution
							{

								strain_region_G_local.yz[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );

							}
							else // Subsequent order corrections
							{
							
								// Calculating correction from order
								strain_correction_region_G_local.yz[idx] = strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );
								// Adding contribution from order
								strain_region_G_local.yz[idx] += strain_correction_region_G_local.yz[idx];

							}
							
						}


						else if (component == 4) // component = 4 --> xz
						{

							if ( correction_order == - 1 ) // Initial contribution
							{

								strain_region_G_local.xz[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );

							}
							else // Subsequent order corrections
							{
							
								// Calculating correction from order
								strain_correction_region_G_local.xz[idx] = strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );
								// Adding contribution from order
								strain_region_G_local.xz[idx] += strain_correction_region_G_local.xz[idx];

							}

						}


						else // component = 5 --> xy
						{
						
							if ( correction_order == - 1 ) // Initial contribution
							{

								strain_region_G_local.xy[idx] += strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );

							}
							else // Subsequent order corrections
							{
							
								// Calculating correction from order
								strain_correction_region_G_local.xy[idx] = strain_fourier_coefficient( component, idx_shape, idx_basis, basis.Gx[idx_Gx], basis.Gy[idx_Gy], basis.Gz[idx_Gz], supercell, parameters, correction_order, &convolutions_char_strain_correction );
								// Adding contribution from order
								strain_region_G_local.xy[idx] += strain_correction_region_G_local.xy[idx];

							}

						}

					}

				}

				if ( correction_order > -1 ) // Gathering local parts of strain correction order to next calculate convolution required for next order correction
				{

					MPI_Allgatherv( strain_correction_region_G_local.xx, local_size, MPI_DOUBLE_COMPLEX, strain_correction_region_G.xx, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
					MPI_Allgatherv( strain_correction_region_G_local.yy, local_size, MPI_DOUBLE_COMPLEX, strain_correction_region_G.yy, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
					MPI_Allgatherv( strain_correction_region_G_local.zz, local_size, MPI_DOUBLE_COMPLEX, strain_correction_region_G.zz, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
					MPI_Allgatherv( strain_correction_region_G_local.yz, local_size, MPI_DOUBLE_COMPLEX, strain_correction_region_G.yz, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
					MPI_Allgatherv( strain_correction_region_G_local.xz, local_size, MPI_DOUBLE_COMPLEX, strain_correction_region_G.xz, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
					MPI_Allgatherv( strain_correction_region_G_local.xy, local_size, MPI_DOUBLE_COMPLEX, strain_correction_region_G.xy, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );

					free( strain_correction_region_G_local.xx );
					free( strain_correction_region_G_local.yy );
					free( strain_correction_region_G_local.zz );
					free( strain_correction_region_G_local.yz );
					free( strain_correction_region_G_local.xz );
					free( strain_correction_region_G_local.xy );

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

	// Gathering local parts of strain tensor and duplicating on all processes
	MPI_Allgatherv( strain_region_G_local.xx, local_size, MPI_DOUBLE_COMPLEX, strain_region_G.xx, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	MPI_Allgatherv( strain_region_G_local.yy, local_size, MPI_DOUBLE_COMPLEX, strain_region_G.yy, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	MPI_Allgatherv( strain_region_G_local.zz, local_size, MPI_DOUBLE_COMPLEX, strain_region_G.zz, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	MPI_Allgatherv( strain_region_G_local.yz, local_size, MPI_DOUBLE_COMPLEX, strain_region_G.yz, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	MPI_Allgatherv( strain_region_G_local.xz, local_size, MPI_DOUBLE_COMPLEX, strain_region_G.xz, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	MPI_Allgatherv( strain_region_G_local.xy, local_size, MPI_DOUBLE_COMPLEX, strain_region_G.xy, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	
	free( strain_region_G_local.xx );
	free( strain_region_G_local.yy );
	free( strain_region_G_local.zz );
	free( strain_region_G_local.yz );
	free( strain_region_G_local.xz );
	free( strain_region_G_local.xy );
	
	free( recvcounts );
	free( displs );

}
