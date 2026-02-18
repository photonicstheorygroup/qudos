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

void write_coulomb_matrix_elements( int Num_VB_States, int Num_CB_States, struct coulomb_matrix_element CME )
{

    int size = 80;
	char *str1;
	str1 = malloc( sizeof(char)*size );
    FILE *fp1;


	// Output Coulomb matrix elements involving hole states only
    strcpy( str1, "coulomb_matrix_elements_hh.dat" );
          
	fp1 = fopen( str1, "w" );

	// Count of uniquely calculated states, total should equal N_CME_hh_unique from 8_band_main.c
	//int state_count = 0;

	int i_s = 0;
	int j_s = 0;
	int k_s = 0;
	int l_s = 0;

	// The Coulomb matrix elements have been calculated and stored according to the defintion of Vijkl in Eq. 20 in 'Plane wave methodology for single quantum dot electronic structure calculations', Vukmirovic, 2008.
	// However, the Coulomb matrix elements are written to file according to the definition of Vijkl in Eq. 2.86 in 'Electronic and optical properties of III-Nitride nanostructures', PhD Thesis, S. Patra, University College Cork (2019).
	// E.g. For the Coulomb matrix elements involving only the hole states, Vijkl from Patra is equal to Viklj from Vukmirovic.

	// Loop each index over hole states
	for( int i = Num_VB_States - 1; i >= 0; i-- )
	{

		for( int k = Num_VB_States - 1; k >= 0; k-- )
		{

			for( int l = Num_VB_States - 1; l >= 0; l-- )
			{

				for( int j = Num_VB_States - 1; j >= 0; j-- )
				{

					// Column 1  - Hole state index i 
					// Column 2  - Hole state index k
					// Column 3  - Hole state index l
					// Column 4  - Hole state index j      			     	                       
					// Column 5  - Real part of Coulomb matrix element [eV]
					// Column 6  - Imag part of Coulomb matrix element [eV]
					// Column 7  - Magnitude of Coulomb matrix element [eV]
				
				
					fprintf( fp1, "%d\t%d\t%d\t%d\t", Num_VB_States - i, Num_VB_States - k, Num_VB_States - l, Num_VB_States - j );
					
					if ( j >= i && k >= i && l >= ( j == i ? k : ( k == i ? j : ( k < j ? i + 1 : i ) ) ) ) // Indicates uniquely calculated element, see index limits due to symmetries in coulomb_integrals.c
					{
					
						fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.hh[i*Num_VB_States*Num_VB_States*Num_VB_States + j*Num_VB_States*Num_VB_States + k*Num_VB_States + l] ), cimag( CME.hh[i*Num_VB_States*Num_VB_States*Num_VB_States + j*Num_VB_States*Num_VB_States + k*Num_VB_States + l] ), cabs( CME.hh[i*Num_VB_States*Num_VB_States*Num_VB_States + j*Num_VB_States*Num_VB_States + k*Num_VB_States + l] ) );

						//state_count++;
						
					}
					else // Check for symmetry property Vkijl=Viklj
					{

						i_s = k;
						j_s = l;
						k_s = i;
						l_s = j;

						if ( j_s >= i_s && k_s >= i_s && l_s >= ( j_s == i_s ? k_s : ( k_s == i_s ? j_s : ( k_s < j_s ? i_s + 1 : i_s ) ) ) )
						{

							fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.hh[i_s*Num_VB_States*Num_VB_States*Num_VB_States + j_s*Num_VB_States*Num_VB_States + k_s*Num_VB_States + l_s] ), cimag( CME.hh[i_s*Num_VB_States*Num_VB_States*Num_VB_States + j_s*Num_VB_States*Num_VB_States + k_s*Num_VB_States + l_s] ), cabs( CME.hh[i_s*Num_VB_States*Num_VB_States*Num_VB_States + j_s*Num_VB_States*Num_VB_States + k_s*Num_VB_States + l_s] ) );

						}
						else // Check for symmetry property Vjlki=Viklj*
						{

							i_s = j;
							j_s = i;
							k_s = l;
							l_s = k;

							if ( j_s >= i_s && k_s >= i_s && l_s >= ( j_s == i_s ? k_s : ( k_s == i_s ? j_s : ( k_s < j_s ? i_s + 1 : i_s ) ) ) )
							{

								fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.hh[i_s*Num_VB_States*Num_VB_States*Num_VB_States + j_s*Num_VB_States*Num_VB_States + k_s*Num_VB_States + l_s] ), -cimag( CME.hh[i_s*Num_VB_States*Num_VB_States*Num_VB_States + j_s*Num_VB_States*Num_VB_States + k_s*Num_VB_States + l_s] ), cabs( CME.hh[i_s*Num_VB_States*Num_VB_States*Num_VB_States + j_s*Num_VB_States*Num_VB_States + k_s*Num_VB_States + l_s] ) );

							}
							else // Symmetry property Vjlki=Viklj*
							{

								i_s = l;
								j_s = k;
								k_s = j;
								l_s = i;

								if ( j_s >= i_s && k_s >= i_s && l_s >= ( j_s == i_s ? k_s : ( k_s == i_s ? j_s : ( k_s < j_s ? i_s + 1 : i_s ) ) ) )
								{

									fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.hh[i_s*Num_VB_States*Num_VB_States*Num_VB_States + j_s*Num_VB_States*Num_VB_States + k_s*Num_VB_States + l_s] ), -cimag( CME.hh[i_s*Num_VB_States*Num_VB_States*Num_VB_States + j_s*Num_VB_States*Num_VB_States + k_s*Num_VB_States + l_s] ), cabs( CME.hh[i_s*Num_VB_States*Num_VB_States*Num_VB_States + j_s*Num_VB_States*Num_VB_States + k_s*Num_VB_States + l_s] ) );

								}

							}

						}
					
					}
					
					fprintf( fp1, "\n");
					
				}
						
				fprintf( fp1, "\n");

			}

		}
		
	}

	fclose( fp1 );


	// Output Coulomb matrix elements involving electron states only
	strcpy( str1, "coulomb_matrix_elements_ee.dat" );

	fp1 = fopen( str1, "w" );

	// Count of uniquely calculated states, total should equal N_CME_ee_unique from 8_band_main.c
	//state_count = 0;

	// Loop each index over electron states
	for( int i = 0; i < Num_CB_States; i++ )
	{

		for( int k = 0; k < Num_CB_States; k++ )
		{

			for( int l = 0; l < Num_CB_States; l++ )
			{

				for( int j = 0; j < Num_CB_States; j++ )
				{

					// Column 1  - Electron state index i 
					// Column 2  - Electron state index k
					// Column 3  - Electron state index l
					// Column 4  - Electron state index j      			     	                       
					// Column 5  - Real part of Coulomb matrix element [eV]
					// Column 6  - Imag part of Coulomb matrix element [eV]
					// Column 7  - Magnitude of Coulomb matrix element [eV]
				
				
					fprintf( fp1, "%d\t%d\t%d\t%d\t", i + 1, k + 1, l + 1, j + 1);
					
					if ( j >= i && k >= i && l >= ( j == i ? k : ( k == i ? j : ( k < j ? i + 1 : i ) ) ) ) // Indicates uniquely calculated element, see index limits due to symmetries in coulomb_integrals.c
					{
					
						fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.ee[i*Num_CB_States*Num_CB_States*Num_CB_States + j*Num_CB_States*Num_CB_States + k*Num_CB_States + l] ), cimag( CME.ee[i*Num_CB_States*Num_CB_States*Num_CB_States + j*Num_CB_States*Num_CB_States + k*Num_CB_States + l] ), cabs( CME.ee[i*Num_CB_States*Num_CB_States*Num_CB_States + j*Num_CB_States*Num_CB_States + k*Num_CB_States + l] ) );

						//state_count++;
						
					}
					else // Check for symmetry property Vkijl=Viklj
					{

						i_s = k;
						j_s = l;
						k_s = i;
						l_s = j;

						if ( j_s >= i_s && k_s >= i_s && l_s >= ( j_s == i_s ? k_s : ( k_s == i_s ? j_s : ( k_s < j_s ? i_s + 1 : i_s ) ) ) )
						{

							fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.ee[i_s*Num_CB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ), cimag( CME.ee[i_s*Num_CB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ), cabs( CME.ee[i_s*Num_CB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ) );

						}
						else // Check for symmetry property Vjlki=Viklj*
						{

							i_s = j;
							j_s = i;
							k_s = l;
							l_s = k;

							if ( j_s >= i_s && k_s >= i_s && l_s >= ( j_s == i_s ? k_s : ( k_s == i_s ? j_s : ( k_s < j_s ? i_s + 1 : i_s ) ) ) )
							{

								fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.ee[i_s*Num_CB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ), -cimag( CME.ee[i_s*Num_CB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ), cabs( CME.ee[i_s*Num_CB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ) );

							}
							else // Symmetry property Vjlki=Viklj*
							{

								i_s = l;
								j_s = k;
								k_s = j;
								l_s = i;

								if ( j_s >= i_s && k_s >= i_s && l_s >= ( j_s == i_s ? k_s : ( k_s == i_s ? j_s : ( k_s < j_s ? i_s + 1 : i_s ) ) ) )
								{

									fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.ee[i_s*Num_CB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ), -cimag( CME.ee[i_s*Num_CB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ), cabs( CME.ee[i_s*Num_CB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ) );

								}

							}

						}
					
					}
					
					fprintf( fp1, "\n");
					
				}
						
				fprintf( fp1, "\n");

			}

		}
		
	}

	fclose( fp1 );


	// Output direct Coulomb matrix elements involving electron and hole states
	strcpy( str1, "coulomb_matrix_elements_eh_direct.dat" );
          
	fp1 = fopen( str1, "w" );

	// Count of uniquely calculated states, total should equal N_CME_eh_unique from 8_band_main.c
	//state_count = 0;

	// Loop each index appropriately over electron or hole states
	for( int k = 0; k < Num_CB_States; k++ )
	{

		for( int i = Num_VB_States - 1; i >= 0; i-- )
		{

			for( int j = Num_VB_States - 1; j >= 0; j-- )
			{

				for( int l = 0; l < Num_CB_States; l++ )
				{

					// Column 1  - Electron state index k 
					// Column 2  - Hole state index i
					// Column 3  - Hole state index j
					// Column 4  - Electron state index l      			     	                       
					// Column 5  - Real part of Coulomb matrix element [eV]
					// Column 6  - Imag part of Coulomb matrix element [eV]
					// Column 7  - Magnitude of Coulomb matrix element [eV]
				
				
					fprintf( fp1, "%d\t%d\t%d\t%d\t", k + 1, Num_VB_States - i, Num_VB_States - j, l + 1);
					
					if ( j >= i && l >= ( j == i ? k : 0 ) ) // Indicates uniquely calculated element, see index limits due to symmetries in coulomb_integrals.c
					{
					
						fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.eh_direct[i*Num_VB_States*Num_CB_States*Num_CB_States + j*Num_CB_States*Num_CB_States + k*Num_CB_States + l] ), cimag( CME.eh_direct[i*Num_VB_States*Num_CB_States*Num_CB_States + j*Num_CB_States*Num_CB_States + k*Num_CB_States + l] ), cabs( CME.eh_direct[i*Num_VB_States*Num_CB_States*Num_CB_States + j*Num_CB_States*Num_CB_States + k*Num_CB_States + l] ) );

						//state_count++;
						
					}
					else // Symmetry property Vjlki=Viklj*
					{

						i_s = j;
						j_s = i;
						k_s = l;
						l_s = k;

						if ( j_s >= i_s && l_s >= ( j_s == i_s ? k_s : 0 ) )
						{

							fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.eh_direct[i_s*Num_VB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ), -cimag( CME.eh_direct[i_s*Num_VB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ), cabs( CME.eh_direct[i_s*Num_VB_States*Num_CB_States*Num_CB_States + j_s*Num_CB_States*Num_CB_States + k_s*Num_CB_States + l_s] ) );

						}
					
					}
					
					fprintf( fp1, "\n");
					
				}
						
				fprintf( fp1, "\n");

			}

		}
		
	}

	fclose( fp1 );


	// Output exchange Coulomb matrix elements involving electron and hole states
	strcpy( str1, "coulomb_matrix_elements_eh_exchange.dat" );
          
	fp1 = fopen( str1, "w" );

	// Count of uniquely calculated states, total should equal N_CME_eh_unique from 8_band_main.c
	//state_count = 0;

	// Loop each index appropriately over electron or hole states
	for( int k = 0; k < Num_CB_States; k++ )
	{

		for( int i = Num_VB_States - 1; i >= 0; i-- )
		{

			for( int j = 0; j < Num_CB_States; j++ )
			{

				for( int l = Num_VB_States - 1; l >= 0; l-- )
				{

					// Column 1  - Electron state index k 
					// Column 2  - Hole state index i
					// Column 3  - Electron state index j
					// Column 4  - Hole state index l      			     	                       
					// Column 5  - Real part of Coulomb matrix element [eV]
					// Column 6  - Imag part of Coulomb matrix element [eV]
					// Column 7  - Magnitude of Coulomb matrix element [eV]
				
				
					fprintf( fp1, "%d\t%d\t%d\t%d\t", k + 1, Num_VB_States - i, j + 1, Num_VB_States - l );
					
					if ( k >= j && l >= ( k == j ? i : 0 ) ) // Indicates uniquely calculated element, see index limits due to symmetries in coulomb_integrals.c
					{
					
						fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.eh_exchange[i*Num_CB_States*Num_CB_States*Num_VB_States + j*Num_CB_States*Num_VB_States + k*Num_VB_States + l] ), cimag( CME.eh_exchange[i*Num_CB_States*Num_CB_States*Num_VB_States + j*Num_CB_States*Num_VB_States + k*Num_VB_States + l] ), cabs( CME.eh_exchange[i*Num_CB_States*Num_CB_States*Num_VB_States + j*Num_CB_States*Num_VB_States + k*Num_VB_States + l] ) );

						//state_count++;
						
					}
					else // Symmetry property Vjlki=Viklj*
					{

						i_s = l;
						j_s = k;
						k_s = j;
						l_s = i;

						if ( k_s >= j_s && l_s >= ( k_s == j_s ? i_s : 0 ) )
						{

							fprintf( fp1, "%17.10e\t%17.10e\t%17.10e\t",   creal( CME.eh_exchange[i_s*Num_CB_States*Num_CB_States*Num_VB_States + j_s*Num_CB_States*Num_VB_States + k_s*Num_VB_States + l_s] ), -cimag( CME.eh_exchange[i_s*Num_CB_States*Num_CB_States*Num_VB_States + j_s*Num_CB_States*Num_VB_States + k_s*Num_VB_States + l_s] ), cabs( CME.eh_exchange[i_s*Num_CB_States*Num_CB_States*Num_VB_States + j_s*Num_CB_States*Num_VB_States + k_s*Num_VB_States + l_s] ) );

						}
					
					}
					
					fprintf( fp1, "\n");
					
				}
						
				fprintf( fp1, "\n");

			}

		}
		
	}

	fclose( fp1 );

	
	free( str1 );

}
