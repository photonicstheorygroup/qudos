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

void write_momentum_matrix_elements( int Num_CB_States, int Num_VB_States, int *polarisations, complex double *momentum_matrix_elements )
{

    int size = 80;
	char *str1;
	str1 = malloc( sizeof(char)*size );
    FILE *fp1;


    strcpy( str1, "momentum_matrix_elements.dat" );
          
	fp1 = fopen( str1, "w" );

	for( int idx_state_cb = 0; idx_state_cb < 2*Num_CB_States; idx_state_cb++ )
	{

		for ( int idx_state_vb = 0; idx_state_vb < 2*Num_VB_States; idx_state_vb++ )
		{
		
			// Column 1  - CB state index 
			// Column 2  - VB state index      			     	                       
			// Column 3  - Momentum matrix element for [100] [eV^(1/2)] - Real part
			// Column 4  - Momentum matrix element for [100] [eV^(1/2)] - Imag part
			// Column 5  - Momentum matrix element for [100] [eV^(1/2)] - Absolute value
			// Column 6  - Momentum matrix element for [010] [eV^(1/2)] - Real part
			// Column 7  - Momentum matrix element for [010] [eV^(1/2)] - Imag part
			// Column 8  - Momentum matrix element for [010] [eV^(1/2)] - Absolute value
			// Column 9  - Momentum matrix element for [001] [eV^(1/2)] - Real part
			// Column 10  - Momentum matrix element for [001] [eV^(1/2)] - Imag part
			// Column 11  - Momentum matrix element for [001] [eV^(1/2)] - Absolute value
			// Column 12  - Momentum matrix element for [110] [eV^(1/2)] - Real part
			// Column 13  - Momentum matrix element for [110] [eV^(1/2)] - Imag part
			// Column 14  - Momentum matrix element for [110] [eV^(1/2)] - Absolute value
			// Column 15  - Momentum matrix element for [-110] [eV^(1/2)] - Real part
			// Column 16  - Momentum matrix element for [-110] [eV^(1/2)] - Imag part
			// Column 17  - Momentum matrix element for [-110] [eV^(1/2)] - Absolute value
		
			// The momentum matrix elements are given in units of eV^(1/2) so that the squared element is given in terms of eV, akin to the Kane energy.
			// As the result of the required vector-matrix-vector product has units of eV Angstrom (as a result of the units the code utilises for electronic structure calculations), the units of eV^(1/2) are acquired by dividing the vector-matrix-vector product by sqrt( HBSQM/2 ), where HBQSM is equal to ( Reduced Planck constant )^2 )/( Free electron mass  ) in units of [eV (Angstrom)^2]. 
		
			fprintf( fp1, "%d\t%d\t", idx_state_cb + 1, idx_state_vb + 1);
		
			for (int polarisation = 0; polarisation < Num_Possible_Polarisations; polarisation++ ) // Loop over polarisations
			{
			
				if ( polarisations[polarisation] == 1 )
				{
				
					fprintf( fp1, "%17.10e\t",   creal( momentum_matrix_elements[polarisation*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] / sqrt( HBSQM/2 ) ) );
					fprintf( fp1, "%17.10e\t",   cimag( momentum_matrix_elements[polarisation*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] / sqrt( HBSQM/2 ) ) );
					fprintf( fp1, "%17.10e\t",   cabs( momentum_matrix_elements[polarisation*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] / sqrt( HBSQM/2 ) ) );
					
				}
				else
				{
				
					fprintf( fp1, "N/A\tN/A\tN/A\t");
				
				}
			
			}

			if ( polarisations[0] == 1 && polarisations[1] == 1 )
			{

				fprintf( fp1, "%17.10e\t",       creal( ( 1.0 / sqrt( 2.0 ) ) * ( momentum_matrix_elements[0*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] + momentum_matrix_elements[1*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] ) / sqrt( HBSQM/2 ) ) );
				fprintf( fp1, "%17.10e\t",       cimag( ( 1.0 / sqrt( 2.0 ) ) * ( momentum_matrix_elements[0*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] + momentum_matrix_elements[1*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] ) / sqrt( HBSQM/2 ) ) );
				fprintf( fp1, "%17.10e\t",       cabs( ( 1.0 / sqrt( 2.0 ) ) * ( momentum_matrix_elements[0*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] + momentum_matrix_elements[1*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] ) / sqrt( HBSQM/2 ) ) );

				fprintf( fp1, "%17.10e\t",       creal( ( 1.0 / sqrt( 2.0 ) ) * ( -momentum_matrix_elements[0*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] + momentum_matrix_elements[1*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] ) / sqrt( HBSQM/2 ) ) );
				fprintf( fp1, "%17.10e\t",       cimag( ( 1.0 / sqrt( 2.0 ) ) * ( -momentum_matrix_elements[0*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] + momentum_matrix_elements[1*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] ) / sqrt( HBSQM/2 ) ) );
				fprintf( fp1, "%17.10e\t",       cabs( ( 1.0 / sqrt( 2.0 ) ) * ( -momentum_matrix_elements[0*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] + momentum_matrix_elements[1*2*Num_CB_States*2*Num_VB_States + idx_state_cb*2*Num_VB_States + idx_state_vb] ) / sqrt( HBSQM/2 ) ) );

			}
			else
			{
				
				fprintf( fp1, "N/A\tN/A\tN/A\tN/A\tN/A\tN/A\t");
				
			}
			
			fprintf( fp1, "\n");
			
		}
				
		fprintf( fp1, "\n");

	}

	fclose( fp1 );
	
	free( str1 );

}
