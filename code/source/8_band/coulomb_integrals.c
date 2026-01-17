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

void coulomb_integrals( struct supercell_geometry *supercell, struct material_params *parameters, int *N_states, struct convolution_state *convolution_states, struct convolution_state *convolution_states_reversed, int local_size_conv, struct coulomb_corrections CME_corrections, struct coulomb_matrix_element CME, int type )
{

	// Note that Vijkl is defined according to Eq. 20 in N. Vukmirovic and S. Tomic, J. Appl. Phys. 103, 103718 (2008), doi: 10.1063/1.2936318
	// This definition corresponds to Viklj in Sec. 2.8 in T. Murphy, PhD Thesis, University College Cork (2025)
	
	// Numbers of Vijkl contributions calculated locally in each block before reducing between processes (MPI_Reduce is called for each i iteration) 
	int N_CME_ee = 2*N_states[1]*2*N_states[1]*2*N_states[1]; // Vijkl involving electron states only
	int N_CME_hh = 2*N_states[0]*2*N_states[0]*2*N_states[0]; // Vijkl involving hole states only
	int N_CME_eh = 2*N_states[1]*2*N_states[0]*2*N_states[1]; // Vijkl involving electron and hole states

	// Will store local contributions to Vijkl calculations
	complex double *CME_local;

	if ( type == 0 ) // Vijkl involving hole states only
	{
	
		CME_local = calloc( N_CME_hh, sizeof( complex double ) );

	}
	else if ( type == 1 ) // Vijkl involving electron states only
	{

		CME_local = calloc( N_CME_ee, sizeof( complex double ) );

	}
	else // Vijkl involving electron and hole states
	{

		CME_local = calloc( N_CME_eh, sizeof( complex double ) );

	}
	

	if ( type == 0 ) // Vijkl involving hole states only
	{

		for ( int i = 0; i < 2*N_states[0]; i++ )
		{

			for ( int j = i; j < 2*N_states[0]; j++ ) // Start index modified to account for symmetry properties
			{

				for ( int k = i; k < 2*N_states[0]; k++ ) // Start index modified to account for symmetry properties
				{

					for ( int l = ( j == i ? k : ( k == i ? j : ( k < j ? i + 1 : i ) ) ); l < 2*N_states[0]; l++ ) // Start index modified to account for symmetry properties
					{

						// Calculating dot product
						if ( l >= k )
						{

							cblas_zdotu_sub( local_size_conv, convolution_states_reversed[(l * (l + 1) / 2) + k].conv_state, 1, convolution_states[(j * (j + 1) / 2) + i].conv_state, 1, &CME_local[j*2*N_states[0]*2*N_states[0] + k*2*N_states[0] + l] );

						}
						else
						{

							cblas_zdotc_sub( local_size_conv, convolution_states_reversed[(k * (k + 1) / 2) + l].conv_state, 1, convolution_states[(j * (j + 1) / 2) + i].conv_state, 1, &CME_local[j*2*N_states[0]*2*N_states[0] + k*2*N_states[0] + l] );

						}

						// If monopole correction is included, add local contribution
						if ( include_CME >= 2 ) 
						{

							if ( i == j && k == l )
							{

								CME_local[j*2*N_states[0]*2*N_states[0] + k*2*N_states[0] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*CME_corrections.a_mad;

							}

						}

						// If dipole correction is included, add local contribution
						if ( include_CME >= 3 && myid == 0 )
						{

							CME_local[j*2*N_states[0]*2*N_states[0] + k*2*N_states[0] + l] += ( ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 4*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 3 ) ) )*( CME_corrections.dx[j * (2*N_states[0] + 2*N_states[1]) + i]*CME_corrections.dx[l * (2*N_states[0] + 2*N_states[1]) + k] + CME_corrections.dy[j * (2*N_states[0] + 2*N_states[1]) + i]*CME_corrections.dy[l * (2*N_states[0] + 2*N_states[1]) + k] + CME_corrections.dz[j * (2*N_states[0] + 2*N_states[1]) + i]*CME_corrections.dz[l * (2*N_states[0] + 2*N_states[1]) + k] );

						}

						// If quadrupole correction is included, add local contribution
						if ( include_CME >= 4 )
						{

							if ( i == j && l >= k )
							{

								CME_local[j*2*N_states[0]*2*N_states[0] + k*2*N_states[0] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 2*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 2 ) ) )*( CME_corrections.Q[(l * (l + 1) / 2) + k] );

							}
							else if ( i == j && l < k )
							{

								CME_local[j*2*N_states[0]*2*N_states[0] + k*2*N_states[0] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 2*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 2 ) ) )*( conj( CME_corrections.Q[(k * (k + 1) / 2) + l] ) );

							}
							
							if ( k == l )
							{

								CME_local[j*2*N_states[0]*2*N_states[0] + k*2*N_states[0] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 2*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 2 ) ) )*( CME_corrections.Q[(j * (j + 1) / 2) + i] );

							}

						}

					}

				}

			}

			// Sum contributions from each process for each Vjkl block 
			MPI_Reduce( CME_local, &CME.hh[i*2*N_states[0]*2*N_states[0]*2*N_states[0]], N_CME_hh, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );

		}

	}
	else if ( type == 1 ) // Vijkl involving electron states only
	{

		for ( int i = 0; i < 2*N_states[1]; i++ )
		{

			for ( int j = i; j < 2*N_states[1]; j++ ) // Start index modified to account for symmetry properties
			{

				for ( int k = i; k < 2*N_states[1]; k++ ) // Start index modified to account for symmetry properties
				{

					for ( int l = ( j == i ? k : ( k == i ? j : ( k < j ? i + 1 : i ) ) ); l < 2*N_states[1]; l++ ) // Start index modified to account for symmetry properties
					{

						if ( l >= k ) // Calculating dot product
						{

							cblas_zdotu_sub( local_size_conv, convolution_states_reversed[(( 2*N_states[0] + l ) * ( 2*N_states[0] + l + 1 ) / 2) + 2*N_states[0] + k].conv_state, 1, convolution_states[(( 2*N_states[0] + j ) * ( 2*N_states[0] + j + 1 ) / 2) + 2*N_states[0] + i].conv_state, 1, &CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] );

						}
						else
						{

							cblas_zdotc_sub( local_size_conv, convolution_states_reversed[(( 2*N_states[0] + k ) * ( 2*N_states[0] + k + 1 ) / 2) + 2*N_states[0] + l].conv_state, 1, convolution_states[(( 2*N_states[0] + j ) * ( 2*N_states[0] + j + 1 ) / 2) + 2*N_states[0] + i].conv_state, 1, &CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] );

						}

						if ( include_CME >= 2 ) // If monopole correction is included, add local contribution
						{

							if ( i == j && k == l )
							{

								CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*CME_corrections.a_mad;

							}

						}

						if ( include_CME >= 3 && myid == 0 ) // If dipole correction is included, add contribution which is entirely stored on process 0
						{

							CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] += ( ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 4*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 3 ) ) )*( CME_corrections.dx[( 2*N_states[0] + j ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + i]*CME_corrections.dx[( 2*N_states[0] + l ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + k] + CME_corrections.dy[( 2*N_states[0] + j ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + i]*CME_corrections.dy[( 2*N_states[0] + l ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + k] + CME_corrections.dz[( 2*N_states[0] + j ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + i]*CME_corrections.dz[( 2*N_states[0] + l ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + k] );

						}

						if ( include_CME >= 4 ) // If quadrupole correction is included, add local contribution
						{

							if ( i == j && l >= k )
							{

								CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 2*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 2 ) ) )*( CME_corrections.Q[(( 2*N_states[0] + l ) * ( 2*N_states[0] + l + 1) / 2) + 2*N_states[0] + k] );

							}
							else if ( i == j && l < k )
							{

								CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 2*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 2) ) )*( conj( CME_corrections.Q[(( 2*N_states[0] + k ) * ( 2*N_states[0] + k + 1) / 2) + 2*N_states[0] + l] ) );

							}
							
							if ( k == l )
							{

								CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 2*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 2 ) ) )*( CME_corrections.Q[(( 2*N_states[0] + j ) * ( 2*N_states[0] + j + 1) / 2) + 2*N_states[0] + i] );

							}

						}

					}

				}

			}

			// Sum contributions from each process for each Vjkl block
			MPI_Reduce( CME_local, &CME.ee[i*2*N_states[1]*2*N_states[1]*2*N_states[1]], N_CME_ee, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );

		}

	}
	else if ( type == 2 ) // Direct Vijkl involving hole and electron states
	{

		for ( int i = 0; i < 2*N_states[0]; i++ )
		{

			for ( int j = i; j < 2*N_states[0]; j++ ) // Start index modified to account for symmetry properties
			{

				for ( int k = 0; k < 2*N_states[1]; k++ )
				{

					for ( int l = ( j == i ? k : 0 ); l < 2*N_states[1]; l++ ) // Start index modified to account for symmetry properties
					{

						if ( l >= k ) // Calculating dot product
						{

							cblas_zdotu_sub( local_size_conv, convolution_states_reversed[(( 2*N_states[0] + l ) * ( 2*N_states[0] + l + 1 ) / 2) + 2*N_states[0] + k].conv_state, 1, convolution_states[(j * (j + 1) / 2) + i].conv_state, 1, &CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] );

						}
						else
						{

							cblas_zdotc_sub( local_size_conv, convolution_states_reversed[(( 2*N_states[0] + k ) * ( 2*N_states[0] + k + 1 ) / 2) + 2*N_states[0] + l].conv_state, 1, convolution_states[(j * (j + 1) / 2) + i].conv_state, 1, &CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] );

						}

						if ( include_CME >= 2 ) // If monopole correction is included, add local contribution
						{

							if ( i == j && k == l )
							{

								CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*CME_corrections.a_mad;

							}

						}

						if ( include_CME >= 3 && myid == 0 ) // If dipole correction is included, add contribution which is entirely stored on process 0
						{

							CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] += ( ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 4*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 3 ) ) )*( CME_corrections.dx[j * (2*N_states[0] + 2*N_states[1]) + i]*CME_corrections.dx[( 2*N_states[0] + l ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + k] + CME_corrections.dy[j * (2*N_states[0] + 2*N_states[1]) + i]*CME_corrections.dy[( 2*N_states[0] + l ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + k] + CME_corrections.dz[j * (2*N_states[0] + 2*N_states[1]) + i]*CME_corrections.dz[( 2*N_states[0] + l ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + k] );

						}

						if ( include_CME >= 4 ) // If quadrupole correction is included, add local contribution
						{

							if ( i == j && l >= k )
							{

								CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 2*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 2 ) ) )*( CME_corrections.Q[(( 2*N_states[0] + l ) * ( 2*N_states[0] + l + 1 ) / 2) + 2*N_states[0] + k] );

							}
							else if ( i == j && l < k )
							{

								CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 2*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 2 ) ) )*( conj( CME_corrections.Q[(( 2*N_states[0] + k ) * ( 2*N_states[0] + k + 1 ) / 2) + 2*N_states[0] + l] ) );

							}
							
							if ( k == l )
							{

								CME_local[j*2*N_states[1]*2*N_states[1] + k*2*N_states[1] + l] += ( -ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 2*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 2 ) ) )*( CME_corrections.Q[(j * (j + 1) / 2) + i] );

							}

						}

					}

				}

			}

			// Sum contributions from each process for each Vjkl block
			MPI_Reduce( CME_local, &CME.eh_direct[i*2*N_states[0]*2*N_states[1]*2*N_states[1]], N_CME_eh, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );

		}

	}
	else if ( type == 3) // Exchange Vijkl involving hole and electron states
	{

		for ( int i = 0; i < 2*N_states[0]; i++ )
		{

			for ( int j = 0; j < 2*N_states[1]; j++ )
			{

				for ( int k = j; k < 2*N_states[1]; k++ ) // Start index modified to account for symmetry properties
				{

					for ( int l = ( k == j ? i : 0 ); l < 2*N_states[0]; l++ ) // Start index modified to account for symmetry properties
					{

						// Calculating dot product
						cblas_zdotc_sub( local_size_conv, convolution_states_reversed[(( 2*N_states[0] + k ) * ( 2*N_states[0] + k + 1 ) / 2) + l].conv_state, 1, convolution_states[(( 2*N_states[0] + j ) * ( 2*N_states[0] + j + 1 ) / 2) + i].conv_state, 1, &CME_local[j*2*N_states[1]*2*N_states[0] + k*2*N_states[0] + l] );
						

						if ( include_CME >= 3 && myid == 0 ) // If dipole correction is included, add contribution which is entirely stored on process 0
						{

							CME_local[j*2*N_states[1]*2*N_states[0] + k*2*N_states[0] + l] += ( ESQEPS / ( 4*M_PI*pow( parameters[0].nr, 2.0 ) ) )*( 4*M_PI / ( 3*pow( supercell[0].CL1*supercell[0].CL2*supercell[0].CL3, 3 ) ) )*( CME_corrections.dx[( 2*N_states[0] + j ) * (2*N_states[0] + 2*N_states[1]) + i]*CME_corrections.dx[( l ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + k] + CME_corrections.dy[( 2*N_states[0] + j ) * (2*N_states[0] + 2*N_states[1]) + i]*CME_corrections.dy[( l ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + k] + CME_corrections.dz[( 2*N_states[0] + j ) * (2*N_states[0] + 2*N_states[1]) + i]*CME_corrections.dz[( l ) * (2*N_states[0] + 2*N_states[1]) + 2*N_states[0] + k] );

						}

					}

				}

			}

			// Sum contributions from each process for each Vjkl block
			MPI_Reduce( CME_local, &CME.eh_exchange[i*2*N_states[1]*2*N_states[1]*2*N_states[0]], N_CME_eh, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD );

		}

	}

	free( CME_local );

}
