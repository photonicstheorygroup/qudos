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

void strain_reciprocal_space( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell )
{
	
	double strain_volume_inverse = 1.0 / ( supercell[0].L1*supercell[0].L2*supercell[0].L3 );
	
	// Local size of involved arrays
	int local_size = basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs + ( ( myid < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? 1 : 0 );
	
	// Will store local components of strain tensor
	struct strain_tensor strain_G_local;
	
	strain_G_local.xx = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	strain_G_local.yy = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	strain_G_local.zz = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	strain_G_local.yz = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	strain_G_local.xz = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	strain_G_local.xy = ( complex double* ) calloc( local_size, sizeof( complex double ) );
	
	
	// Will store number of strain tensor components calculated by each process
	int *recvcounts = ( int * ) malloc( nprocs*sizeof( int ) );
	
	// Will store displacement from first value of the "global" strain tensor components for the first value of the strain tensor calculated by each process
	int *displs = ( int * ) malloc( nprocs*sizeof( int ) );
	
	for ( int proc = 0; proc < nprocs; proc++ )
	{
	
		recvcounts[proc] = basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs + ( ( proc < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? 1 : 0 );
		displs[proc] = ( proc == 0 ? 0 : displs[proc - 1] + recvcounts[proc - 1] ); 
	
	}
	
	// Loop over local wave vectors and evaluate the Fourier coefficients of the chosen component of the strain tensor
	#pragma omp parallel for
	for( int idx = 0; idx < local_size; idx++ )
	{

		// Global wave vector index for electronic structure supercell
		int idx_basis = myid * ( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T / nprocs ) + ( ( myid < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) ? myid : basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T % nprocs ) + idx;
		
		int idx_Gx = idx_basis / ( basis.N_Gy_T*basis.N_Gz_T );
	
		int idx_Gy = ( idx_basis % ( basis.N_Gy_T*basis.N_Gz_T ) ) / basis.N_Gz_T;
	
		int idx_Gz = ( idx_basis % ( basis.N_Gy_T*basis.N_Gz_T ) ) % basis.N_Gz_T;

		strain_G_local.xx[idx] = 0.0 + I*0.0;
		strain_G_local.yy[idx] = 0.0 + I*0.0;
		strain_G_local.zz[idx] = 0.0 + I*0.0;
		strain_G_local.yz[idx] = 0.0 + I*0.0;
		strain_G_local.xz[idx] = 0.0 + I*0.0;
		strain_G_local.xy[idx] = 0.0 + I*0.0;
		
		for( int idx_GSx = 0; idx_GSx < strain_basis.N_Gx_T; idx_GSx++ )
		{
		
			if ( basis.N_Gx == strain_basis.N_Gx && supercell[0].L1 == supercell[0].SL1 ) // Copying in the case where the strain supercell is the same size as the electronic structure supercell
			{
			
				idx_GSx = idx_Gx;
			
			}

			for( int idx_GSy = 0; idx_GSy < strain_basis.N_Gy_T; idx_GSy++ )
			{
			
				if ( basis.N_Gy == strain_basis.N_Gy && supercell[0].L2 == supercell[0].SL2 ) // Copying in the case where the strain supercell is the same size as the electronic structure supercell
				{
				
					idx_GSy = idx_Gy;
				
				}

				for( int idx_GSz = 0; idx_GSz < strain_basis.N_Gz_T; idx_GSz++ )
				{
				
					if ( basis.N_Gz == strain_basis.N_Gz && supercell[0].L3 == supercell[0].SL3 ) // Copying in the case where the strain supercell is the same size as the electronic structure supercell
					{
					
						idx_GSz = idx_Gz;
					
					}
				
					// Global wave vector index for strain supercell
					int idx_basis_strain = idx_GSz + strain_basis.N_Gz_T*idx_GSy + strain_basis.N_Gy_T*strain_basis.N_Gz_T*idx_GSx;
					
					double Gx = strain_basis.Gx[idx_GSx] - basis.Gx[idx_Gx];
					double Gy = strain_basis.Gy[idx_GSy] - basis.Gy[idx_Gy];
					double Gz = strain_basis.Gz[idx_GSz] - basis.Gz[idx_Gz];
					
					double I0x;
					double I0y;
					double I0z;
					
					if ( fabs( Gx ) < EPS ) // Gx  = 0
					{
					
						I0x = supercell[0].L1;
					
					}
					else
					{
					
						I0x = 2*sin( Gx*supercell[0].L1/2 ) / Gx;
					
					}
					
					if ( fabs( Gy ) < EPS ) // Gy  = 0
					{
					
						I0y = supercell[0].L2;
					
					}
					else
					{
					
						I0y = 2*sin( Gy*supercell[0].L2/2 ) / Gy;
					
					}
					
					if ( fabs( Gz ) < EPS ) // Gz  = 0
					{
					
						I0z = supercell[0].L3;
					
					}
					else
					{
					
						I0z = 2*sin( Gz*supercell[0].L3/2 ) / Gz;
					
					}

					// component = 0 --> xx

					strain_G_local.xx[idx] += strain_volume_inverse*strain_region_G.xx[idx_basis_strain]*I0x*I0y*I0z;


					// component = 1 --> yy

					strain_G_local.yy[idx] += strain_volume_inverse*strain_region_G.yy[idx_basis_strain]*I0x*I0y*I0z;


					// component = 2 --> zz

					strain_G_local.zz[idx] += strain_volume_inverse*strain_region_G.zz[idx_basis_strain]*I0x*I0y*I0z;

					// component = 3 --> yz

					strain_G_local.yz[idx] += strain_volume_inverse*strain_region_G.yz[idx_basis_strain]*I0x*I0y*I0z;

					// component = 4 --> xz

					strain_G_local.xz[idx] += strain_volume_inverse*strain_region_G.xz[idx_basis_strain]*I0x*I0y*I0z;

					// component = 5 --> xy

					strain_G_local.xy[idx] += strain_volume_inverse*strain_region_G.xy[idx_basis_strain]*I0x*I0y*I0z;
					
					if ( basis.N_Gz == strain_basis.N_Gz && supercell[0].L3 == supercell[0].SL3 )
					{
					
						break;
					
					}
					
				}
				
				if ( basis.N_Gy == strain_basis.N_Gy && supercell[0].L2 == supercell[0].SL2 )
				{
				
					break;
				
				}
				
			}
			
			if ( basis.N_Gx == strain_basis.N_Gx && supercell[0].L1 == supercell[0].SL1 )
			{
			
				break;
			
			}
			
		}
		
	}
	
	// Gathering strain tensor components calculated for electronic structure supercell
	MPI_Allgatherv( strain_G_local.xx, local_size, MPI_DOUBLE_COMPLEX, strain_G.xx, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	MPI_Allgatherv( strain_G_local.yy, local_size, MPI_DOUBLE_COMPLEX, strain_G.yy, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	MPI_Allgatherv( strain_G_local.zz, local_size, MPI_DOUBLE_COMPLEX, strain_G.zz, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	MPI_Allgatherv( strain_G_local.yz, local_size, MPI_DOUBLE_COMPLEX, strain_G.yz, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	MPI_Allgatherv( strain_G_local.xz, local_size, MPI_DOUBLE_COMPLEX, strain_G.xz, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	MPI_Allgatherv( strain_G_local.xy, local_size, MPI_DOUBLE_COMPLEX, strain_G.xy, recvcounts, displs, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );
	
	free( strain_G_local.xx );
	free( strain_G_local.yy );
	free( strain_G_local.zz );
	free( strain_G_local.yz );
	free( strain_G_local.xz );
	free( strain_G_local.xy );
	
	free( recvcounts );
	free( displs );

}
