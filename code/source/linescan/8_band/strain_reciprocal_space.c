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

void strain_reciprocal_space( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell, struct strain_tensor strain_G, struct strain_tensor strain_G_expanded )
{
	
	double strain_volume_inverse = 1.0 / ( supercell[0].L1*supercell[0].L2*supercell[0].L3 );

	#pragma omp parallel
	{

		// Loop over electronic structure supercell basis set
		#pragma omp for collapse( 3 )
		for( int idx_Gx = 0; idx_Gx < basis.N_Gx_T; idx_Gx++ )
		{

			for( int idx_Gy = 0; idx_Gy < basis.N_Gy_T; idx_Gy++ )
			{

				for( int idx_Gz = 0; idx_Gz < basis.N_Gz_T; idx_Gz++ )
				{

					// Wave vector index of electronic structure supercell strain arrays
					int idx_basis = idx_Gz + basis.N_Gz_T*idx_Gy + basis.N_Gy_T*basis.N_Gz_T*idx_Gx;

					strain_G.xx[idx_basis] = 0.0 + I*0.0;
					strain_G.yy[idx_basis] = 0.0 + I*0.0;
					strain_G.zz[idx_basis] = 0.0 + I*0.0;
					strain_G.yz[idx_basis] = 0.0 + I*0.0;
					strain_G.xz[idx_basis] = 0.0 + I*0.0;
					strain_G.xy[idx_basis] = 0.0 + I*0.0;
					
					// Loop over strain supercell basis set
					for( int idx_GSx = 0; idx_GSx < strain_basis.N_Gx_T; idx_GSx++ )
					{
					
						if ( basis.N_Gx == strain_basis.N_Gx && supercell[0].L1 == supercell[0].SL1 ) // If x dimension lengths are the same
						{
						
							idx_GSx = idx_Gx;
						
						}

						for( int idx_GSy = 0; idx_GSy < strain_basis.N_Gy_T; idx_GSy++ ) 
						{
						
							if ( basis.N_Gy == strain_basis.N_Gy && supercell[0].L2 == supercell[0].SL2 ) // If y dimension lengths are the same
							{
							
								idx_GSy = idx_Gy;
							
							}

							for( int idx_GSz = 0; idx_GSz < strain_basis.N_Gz_T; idx_GSz++ )
							{
							
								if ( basis.N_Gz == strain_basis.N_Gz && supercell[0].L3 == supercell[0].SL3 ) // If z dimension lengths are the same
								{
								
									idx_GSz = idx_Gz;
								
								}
							
								// Wave vector index of strain supercell strain arrays
								int idx_basis_strain = idx_GSz + strain_basis.N_Gz_T*idx_GSy + strain_basis.N_Gy_T*strain_basis.N_Gz_T*idx_GSx;
								
								// Wave vector components
								double Gx = strain_basis.Gx[idx_GSx] - basis.Gx[idx_Gx];
								double Gy = strain_basis.Gy[idx_GSy] - basis.Gy[idx_Gy];
								double Gz = strain_basis.Gz[idx_GSz] - basis.Gz[idx_Gz];
								
								double I0x;
								double I0y;
								double I0z;
								
								if ( fabs( Gx ) < EPS )
								{
								
									I0x = supercell[0].L1;
								
								}
								else
								{
								
									I0x = 2*sin( Gx*supercell[0].L1/2 ) / Gx;
								
								}
								
								if ( fabs( Gy ) < EPS )
								{
								
									I0y = supercell[0].L2;
								
								}
								else
								{
								
									I0y = 2*sin( Gy*supercell[0].L2/2 ) / Gy;
								
								}
								
								if ( fabs( Gz ) < EPS )
								{
								
									I0z = supercell[0].L3;
								
								}
								else
								{
								
									I0z = 2*sin( Gz*supercell[0].L3/2 ) / Gz;
								
								}

								// component = 0 --> xx

								strain_G.xx[idx_basis] += strain_volume_inverse*strain_G_expanded.xx[idx_basis_strain]*I0x*I0y*I0z;


								// component = 1 --> yy

								strain_G.yy[idx_basis] += strain_volume_inverse*strain_G_expanded.yy[idx_basis_strain]*I0x*I0y*I0z;


								// component = 2 --> zz

								strain_G.zz[idx_basis] += strain_volume_inverse*strain_G_expanded.zz[idx_basis_strain]*I0x*I0y*I0z;

								// component = 3 --> yz

								strain_G.yz[idx_basis] += strain_volume_inverse*strain_G_expanded.yz[idx_basis_strain]*I0x*I0y*I0z;

								// component = 4 --> xz

								strain_G.xz[idx_basis] += strain_volume_inverse*strain_G_expanded.xz[idx_basis_strain]*I0x*I0y*I0z;

								// component = 5 --> xy

								strain_G.xy[idx_basis] += strain_volume_inverse*strain_G_expanded.xy[idx_basis_strain]*I0x*I0y*I0z;
								
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

			}

		}
	}

}
