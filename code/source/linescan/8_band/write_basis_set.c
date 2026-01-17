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

void write_basis_set( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell )
{

	FILE *fp1;

	fp1 = fopen( "basis_set.dat", "w" );

	fprintf( fp1, "Electronic structure supercell size:\n\n" );

	fprintf( fp1, "\t(Lx, Ly, Lz) = (%.4e, %.4e, %.4e) [nm]\n", 0.1*supercell[0].L1, 0.1*supercell[0].L2, 0.1*supercell[0].L3 );
	fprintf( fp1, "\tVolume       =  %.4e [nm^3]           \n", 0.1*supercell[0].L1 *0.1*supercell[0].L2 *0.1*supercell[0].L3 );
	fprintf( fp1, "\n\n" );
	
	fprintf( fp1, "Strain supercell size:\n\n" );

	fprintf( fp1, "\t(LSx, LSy, LSz) = (%.4e, %.4e, %.4e) [nm]\n", 0.1*supercell[0].SL1, 0.1*supercell[0].SL2, 0.1*supercell[0].SL3 );
	fprintf( fp1, "\tVolume       =  %.4e [nm^3]           \n", 0.1*supercell[0].SL1 *0.1*supercell[0].SL2 *0.1*supercell[0].SL3 );
	fprintf( fp1, "\n\n" );

	fprintf( fp1, "Plane wave vectors:\n\n" );

	fprintf( fp1, "\t(M_Gx, M_Gy, M_Gz) = (%d, %d, %d)\n", (int) ( basis.N_Gx - 1 )/2, (int) ( basis.N_Gy - 1 )/2, (int) ( basis.N_Gz - 1 )/2 );
	fprintf( fp1, "\t(N_Gx, N_Gy, N_Gz) = (%d, %d, %d) = (2*M_Gx + 1, 2*M_Gy + 1, 2*M_Gz + 1)\n", basis.N_Gx, basis.N_Gy, basis.N_Gz );
	fprintf( fp1, "\n" );
	fprintf( fp1, "\tmax(|Gx|, |Gy|, |Gz|) = (%.4e, %.4e, %.4e) [nm^-1]\n", 10.0*basis.Gx[((basis.N_Gx_T - 1) / 2) + (basis.N_Gx/2)],
															  10.0*basis.Gy[((basis.N_Gy_T - 1) / 2) + (basis.N_Gy/2)],
															  10.0*basis.Gz[((basis.N_Gz_T - 1) / 2) + (basis.N_Gz/2)] );
	fprintf( fp1, "\n" );
	fprintf( fp1, "\tCut-off energy = %.4e [eV]\n", 0.5*HBSQM*( basis.Gx[((basis.N_Gx_T - 1) / 2) + (basis.N_Gx/2)]*basis.Gx[((basis.N_Gx_T - 1) / 2) + (basis.N_Gx/2)]
												   + basis.Gy[((basis.N_Gy_T - 1) / 2) + (basis.N_Gy/2)]*basis.Gy[((basis.N_Gy_T - 1) / 2) + (basis.N_Gy/2)]
												   + basis.Gz[((basis.N_Gz_T - 1) / 2) + (basis.N_Gz/2)]*basis.Gz[((basis.N_Gz_T - 1) / 2) + (basis.N_Gz/2)] ) );
												   
	fprintf( fp1, "Strain supercell plane wave vectors:\n\n" );

	fprintf( fp1, "\t(M_S_Gx, M_S_Gy, M_S_Gz) = (%d, %d, %d)\n", (int) ( strain_basis.N_Gx - 1 )/2, (int) ( strain_basis.N_Gy - 1 )/2, (int) ( strain_basis.N_Gz - 1 )/2 );
	fprintf( fp1, "\t(N_S_Gx, N_S_Gy, N_S_Gz) = (%d, %d, %d) = (2*M_S_Gx + 1, 2*M_S_Gy + 1, 2*M_S_Gz + 1)\n", strain_basis.N_Gx, strain_basis.N_Gy, strain_basis.N_Gz );
	fprintf( fp1, "\n" );
	fprintf( fp1, "\tmax(|GSx|, |GSy|, |GSz|) = (%.4e, %.4e, %.4e) [nm^-1]\n", 10.0*strain_basis.Gx[((strain_basis.N_Gx_T - 1) / 2) + (strain_basis.N_Gx/2)],
															  10.0*strain_basis.Gy[((strain_basis.N_Gy_T - 1) / 2) + (strain_basis.N_Gy/2)],
															  10.0*strain_basis.Gz[((strain_basis.N_Gz_T - 1) / 2) + (strain_basis.N_Gz/2)] );
	fprintf( fp1, "\n" );

	fclose( fp1 );


}
