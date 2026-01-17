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

void B_coulomb_supercell( struct wave_vectors basis, struct wave_vectors coulomb_basis, struct supercell_geometry *supercell, int local_size_coulomb, struct convolution_state convolution_states_interim, struct convolution_state convolution_states_interim_local, int *recvcounts_coulomb, int *displs_coulomb, struct convolution_state convolution_states_global )
{
	
	double coulomb_volume_inverse = 1.0 / ( supercell[0].L1*supercell[0].L2*supercell[0].L3 );
	
	// Loop over local wave vectors and evaluate the Fourier coefficients of the B function for the Coulomb supercell
	#pragma omp parallel for
	for( int idx = 0; idx < local_size_coulomb; idx++ )
	{

		// Global wave vector index for Coulomb supercell
		int idx_basis_coulomb = myid * ( coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T / nprocs ) + ( ( myid < coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T % nprocs ) ? myid : coulomb_basis.N_Gx_T*coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T % nprocs ) + idx;
		
		int idx_GCx = idx_basis_coulomb / ( coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T );
	
		int idx_GCy = ( idx_basis_coulomb % ( coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T ) ) / coulomb_basis.N_Gz_T;
	
		int idx_GCz = ( idx_basis_coulomb % ( coulomb_basis.N_Gy_T*coulomb_basis.N_Gz_T ) ) % coulomb_basis.N_Gz_T;

		convolution_states_interim_local.conv_state[idx] = 0.0 + I*0.0;
		
		for( int idx_Gx = 0; idx_Gx < basis.N_Gx_T; idx_Gx++ )
		{
		
			if ( coulomb_basis.N_Gx == basis.N_Gx && supercell[0].L1 == supercell[0].CL1 ) // Copying in the case where the Coulomb supercell is the same size as the electronic structure supercell
			{
			
				idx_Gx = idx_GCx;
			
			}

			for( int idx_Gy = 0; idx_Gy < basis.N_Gy_T; idx_Gy++ )
			{
			
				if ( coulomb_basis.N_Gy == basis.N_Gy && supercell[0].L2 == supercell[0].CL2 ) // Copying in the case where the strain supercell is the same size as the electronic structure supercell
				{
				
					idx_Gy = idx_GCy;
				
				}

				for( int idx_Gz = 0; idx_Gz < basis.N_Gz_T; idx_Gz++ )
				{
				
					if ( coulomb_basis.N_Gz == basis.N_Gz && supercell[0].L3 == supercell[0].CL3 ) // Copying in the case where the strain supercell is the same size as the electronic structure supercell
					{
					
						idx_Gz = idx_GCz;
					
					}
				
					int idx_basis = idx_Gz + basis.N_Gz_T*idx_Gy + basis.N_Gy_T*basis.N_Gz_T*idx_Gx; // Global wave vector index for electronic structure supercell
					
					double Gx = coulomb_basis.Gx[idx_GCx] - basis.Gx[idx_Gx];
					double Gy = coulomb_basis.Gy[idx_GCy] - basis.Gy[idx_Gy];
					double Gz = coulomb_basis.Gz[idx_GCz] - basis.Gz[idx_Gz];
					
					double I0x;
					double I0y;
					double I0z;
					
					if ( fabs( Gx ) < EPS ) // Gx = 0
					{
					
						I0x = supercell[0].L1;
					
					}
					else if ( fabs( Gx ) > coulomb_basis.Gx[coulomb_basis.N_Gx_T - 1] + EPS ) // Upsampling wave-vector cut-off to be consistent with FFT method
					{
					
						I0x = 0;
					
					}
					else
					{
					
						I0x = 2*sin( Gx*supercell[0].L1/2 ) / Gx;
					
					}
					
					if ( fabs( Gy ) < EPS ) // Gy = 0
					{
					
						I0y = supercell[0].L2;
					
					}
					else if ( fabs( Gy ) > coulomb_basis.Gy[coulomb_basis.N_Gy_T - 1] + EPS ) // Upsampling wave-vector cut-off to be consistent with FFT method
					{
					
						I0y = 0;
					
					}
					else
					{
					
						I0y = 2*sin( Gy*supercell[0].L2/2 ) / Gy;
					
					}
					
					if ( fabs( Gz ) < EPS )  // Gz = 0
					{
					
						I0z = supercell[0].L3;
					
					}
					else if ( fabs( Gz ) > coulomb_basis.Gz[coulomb_basis.N_Gz_T - 1] + EPS ) // Upsampling wave-vector cut-off to be consistent with FFT method
					{
					
						I0z = 0;
					
					}
					else
					{
					
						I0z = 2*sin( Gz*supercell[0].L3/2 ) / Gz;
					
					}

					
					convolution_states_interim_local.conv_state[idx] += coulomb_volume_inverse*convolution_states_interim.conv_state[idx_basis]*I0x*I0y*I0z;

					
					if ( coulomb_basis.N_Gz == basis.N_Gz && supercell[0].L3 == supercell[0].CL3 )
					{
					
						break;
					
					}
					
				}
				
				if ( coulomb_basis.N_Gy == basis.N_Gy && supercell[0].L2 == supercell[0].CL2 )
				{
				
					break;
				
				}
				
			}
			
			if ( coulomb_basis.N_Gx == basis.N_Gx && supercell[0].L1 == supercell[0].CL1 )
			{
			
				break;
			
			}
			
		}
		
	}
	
	// Gathering B components calculated for Coulomb supercell
	MPI_Allgatherv( convolution_states_interim_local.conv_state, local_size_coulomb, MPI_DOUBLE_COMPLEX, convolution_states_global.conv_state, recvcounts_coulomb, displs_coulomb, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD );

}
