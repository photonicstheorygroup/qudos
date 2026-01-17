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

void write_calculation_times()
{

     int size = 80;
	char *str1;
	str1 = malloc( sizeof(char)*size );
     FILE *fp1;

     int *recvcounts = ( int * ) malloc( nprocs*sizeof( int ) );
	int *displs = ( int * ) malloc( nprocs*sizeof( int ) );

	for ( int proc = 0; proc < nprocs; proc++ )
	{

		recvcounts[proc] = 1;

		displs[proc] = ( proc == 0 ? 0 : displs[proc - 1] + recvcounts[proc - 1] );
	
	}

     double *calculation_times = ( double * ) malloc( nprocs*sizeof( double ) );

     // Write calculation times

     if ( myid == 0 )
     {
          
          strcpy( str1, "calculation_times.dat" );
          fp1 = fopen( str1, "w" );

          MPI_Allgatherv( &component_run_times.char_s_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.char_es_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.strain_s_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.piezo_s_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.strain_es_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.piezo_es_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.conv, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.ham_elements, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.ham_diag, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.eig, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.bloch, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.linescan, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.cross_section, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.real_3D, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );
          
          MPI_Allgatherv( &component_run_times.overlaps, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.car_loc, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.mme, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.cme_B, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.cme_I, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.cme_total, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          MPI_Allgatherv( &component_run_times.total, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", calculation_times[idx_proc] );

          }
          fprintf( fp1, "%17.10e\n", calculation_times[nprocs - 1] );

          fclose( fp1 );


     }
     else
     {

          MPI_Allgatherv( &component_run_times.char_s_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
     
          MPI_Allgatherv( &component_run_times.char_es_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.strain_s_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
        
          MPI_Allgatherv( &component_run_times.piezo_s_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.strain_es_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.piezo_es_grid, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.conv, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.ham_elements, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.ham_diag, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.eig, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.bloch, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.linescan, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.cross_section, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.real_3D, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.overlaps, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.car_loc, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          
          MPI_Allgatherv( &component_run_times.mme, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.cme_B, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.cme_I, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.cme_total, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_run_times.total, 1, MPI_DOUBLE, calculation_times, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

     }

     

     // Row 1  - Characteristic functions in reciprocal space for strain region                     [s]
     // Row 2  - Characteristic functions in reciprocal space for electronic structure region       [s]
     // Row 3  - Strain tensor components in reciprocal space for strain region                     [s]
     // Row 4  - Piezoelectric potential in reciprocal space for strain region                      [s]
     // Row 5  - Strain tensor components in reciprocal space for electronic structure region       [s]
     // Row 6  - Piezoelectric potential in reciprocal space for electronic structure region        [s]
     // Row 7  - Convolutions between strain tensor components and characteristic functions         [s]
     // Row 8  - Hamiltonian matrix elements                                                        [s]
     // Row 9  - Diagonalise Hamiltonian                                                            [s]
     // Row 10 - Total eigenstates calculation                                                      [s]
     // Row 11 - Bloch character of states                                                          [s]
     // Row 12 - Linescans                                                                          [s]
     // Row 13 - Cross-sections                                                                     [s]
     // Row 14 - 3D real space output                                                               [s]
     // Row 15 - Overlaps between states                                                            [s]
     // Row 16 - Carrier localisation in shapes                                                     [s]
     // Row 17 - Momentum matrix elements                                                           [s]
     // Row 18 - Coulomb matrix element B convolutions                                              [s]
     // Row 19 - Coulomb matrix element i, j, k, l integrals                                        [s]
     // Row 20 - Total coulomb matrix element time                                                  [s]
     // Row 21 - Total computation time                                                             [s]


     free( recvcounts );
     free( displs );
     free( calculation_times );
     free( str1 );

}
