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

void write_memory_usage()
{

     int idx_proc, size = 80;
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

     long long *memory_usage = ( long long * ) malloc( nprocs*sizeof( long long ) );
     long long *peak_memory_usage = ( long long * ) malloc( nprocs*sizeof( long long ) );

     double GB = ( 1024.0*1024.0*1024.0 );

     // Write calculation memory usage

     if ( myid == 0 )
     {
          
          strcpy( str1, "memory_usage.dat" );
          fp1 = fopen( str1, "w" );

          MPI_Allgatherv( &component_memory_usage.char_mem, 1, MPI_LONG_LONG, memory_usage, recvcounts, displs, MPI_LONG_LONG, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", memory_usage[idx_proc] / GB );
               peak_memory_usage[idx_proc] = memory_usage[idx_proc];

          }
          fprintf( fp1, "%17.10e\n", memory_usage[nprocs - 1] / GB );

          MPI_Allgatherv( &component_memory_usage.strain_mem, 1, MPI_LONG_LONG, memory_usage, recvcounts, displs, MPI_LONG_LONG, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", memory_usage[idx_proc] / GB );
               if ( memory_usage[idx_proc] > peak_memory_usage[idx_proc] )
               {

                    peak_memory_usage[idx_proc] = memory_usage[idx_proc];

               }

          }
          fprintf( fp1, "%17.10e\n", memory_usage[nprocs - 1] / GB );
          if ( memory_usage[nprocs - 1] > peak_memory_usage[nprocs - 1] )
          {

               peak_memory_usage[nprocs - 1] = memory_usage[nprocs - 1];

          }

          MPI_Allgatherv( &component_memory_usage.piezo_mem, 1, MPI_LONG_LONG, memory_usage, recvcounts, displs, MPI_LONG_LONG, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", memory_usage[idx_proc] / GB );
               if ( memory_usage[idx_proc] > peak_memory_usage[idx_proc] )
               {

                    peak_memory_usage[idx_proc] = memory_usage[idx_proc];

               }

          }
          fprintf( fp1, "%17.10e\n", memory_usage[nprocs - 1] / GB );
          if ( memory_usage[nprocs - 1] > peak_memory_usage[nprocs - 1] )
          {

               peak_memory_usage[nprocs - 1] = memory_usage[nprocs - 1];

          }

          MPI_Allgatherv( &component_memory_usage.strain_downsample_mem, 1, MPI_LONG_LONG, memory_usage, recvcounts, displs, MPI_LONG_LONG, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", memory_usage[idx_proc] / GB );
               if ( memory_usage[idx_proc] > peak_memory_usage[idx_proc] )
               {

                    peak_memory_usage[idx_proc] = memory_usage[idx_proc];

               }

          }
          fprintf( fp1, "%17.10e\n", memory_usage[nprocs - 1] / GB );
          if ( memory_usage[nprocs - 1] > peak_memory_usage[nprocs - 1] )
          {

               peak_memory_usage[nprocs - 1] = memory_usage[nprocs - 1];

          }

          MPI_Allgatherv( &component_memory_usage.convolutions_mem, 1, MPI_LONG_LONG, memory_usage, recvcounts, displs, MPI_LONG_LONG, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", memory_usage[idx_proc] / GB );
               if ( memory_usage[idx_proc] > peak_memory_usage[idx_proc] )
               {

                    peak_memory_usage[idx_proc] = memory_usage[idx_proc];

               }

          }
          fprintf( fp1, "%17.10e\n", memory_usage[nprocs - 1] / GB );
          if ( memory_usage[nprocs - 1] > peak_memory_usage[nprocs - 1] )
          {

               peak_memory_usage[nprocs - 1] = memory_usage[nprocs - 1];

          }

          MPI_Allgatherv( &component_memory_usage.eig_work_mem, 1, MPI_LONG_LONG, memory_usage, recvcounts, displs, MPI_LONG_LONG, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", memory_usage[idx_proc] / GB );
               if ( memory_usage[idx_proc] > peak_memory_usage[idx_proc] )
               {

                    peak_memory_usage[idx_proc] = memory_usage[idx_proc];

               }

          }
          fprintf( fp1, "%17.10e\n", memory_usage[nprocs - 1] / GB );
          if ( memory_usage[nprocs - 1] > peak_memory_usage[nprocs - 1] )
          {

               peak_memory_usage[nprocs - 1] = memory_usage[nprocs - 1];

          }

          MPI_Allgatherv( &component_memory_usage.elec_struc_mem, 1, MPI_LONG_LONG, memory_usage, recvcounts, displs, MPI_LONG_LONG, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", memory_usage[idx_proc] / GB );
               if ( memory_usage[idx_proc] > peak_memory_usage[idx_proc] )
               {

                    peak_memory_usage[idx_proc] = memory_usage[idx_proc];

               }

          }
          fprintf( fp1, "%17.10e\n", memory_usage[nprocs - 1] / GB );
          if ( memory_usage[nprocs - 1] > peak_memory_usage[nprocs - 1] )
          {

               peak_memory_usage[nprocs - 1] = memory_usage[nprocs - 1];

          }

          MPI_Allgatherv( &component_memory_usage.mme_mem, 1, MPI_LONG_LONG, memory_usage, recvcounts, displs, MPI_LONG_LONG, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", memory_usage[idx_proc] / GB );
               if ( memory_usage[idx_proc] > peak_memory_usage[idx_proc] )
               {

                    peak_memory_usage[idx_proc] = memory_usage[idx_proc];

               }

          }
          fprintf( fp1, "%17.10e\n", memory_usage[nprocs - 1] / GB );
          if ( memory_usage[nprocs - 1] > peak_memory_usage[nprocs - 1] )
          {

               peak_memory_usage[nprocs - 1] = memory_usage[nprocs - 1];

          }

          MPI_Allgatherv( &component_memory_usage.real_space_mem, 1, MPI_LONG_LONG, memory_usage, recvcounts, displs, MPI_LONG_LONG, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", memory_usage[idx_proc] / GB );
               if ( memory_usage[idx_proc] > peak_memory_usage[idx_proc] )
               {

                    peak_memory_usage[idx_proc] = memory_usage[idx_proc];

               }

          }
          fprintf( fp1, "%17.10e\n", memory_usage[nprocs - 1] / GB );
          if ( memory_usage[nprocs - 1] > peak_memory_usage[nprocs - 1] )
          {

               peak_memory_usage[nprocs - 1] = memory_usage[nprocs - 1];

          }

          MPI_Allgatherv( &component_memory_usage.CME_mem, 1, MPI_LONG_LONG, memory_usage, recvcounts, displs, MPI_LONG_LONG, MPI_COMM_WORLD );
          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", memory_usage[idx_proc] / GB );
               if ( memory_usage[idx_proc] > peak_memory_usage[idx_proc] )
               {

                    peak_memory_usage[idx_proc] = memory_usage[idx_proc];

               }

          }
          fprintf( fp1, "%17.10e\n", memory_usage[nprocs - 1] / GB );
          if ( memory_usage[nprocs - 1] > peak_memory_usage[nprocs - 1] )
          {

               peak_memory_usage[nprocs - 1] = memory_usage[nprocs - 1];

          }

          for ( int idx_proc = 0; idx_proc < nprocs - 1; idx_proc++ )
          {

               fprintf( fp1, "%17.10e\t", peak_memory_usage[idx_proc] / GB );

          }
          fprintf( fp1, "%17.10e\n", peak_memory_usage[nprocs - 1] / GB );


          fclose( fp1 );

     }
     else
     {

          MPI_Allgatherv( &component_memory_usage.char_mem, 1, MPI_DOUBLE, memory_usage, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
     
          MPI_Allgatherv( &component_memory_usage.strain_mem, 1, MPI_DOUBLE, memory_usage, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_memory_usage.piezo_mem, 1, MPI_DOUBLE, memory_usage, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
        
          MPI_Allgatherv( &component_memory_usage.strain_downsample_mem, 1, MPI_DOUBLE, memory_usage, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_memory_usage.convolutions_mem, 1, MPI_DOUBLE, memory_usage, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_memory_usage.eig_work_mem, 1, MPI_DOUBLE, memory_usage, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
          
          MPI_Allgatherv( &component_memory_usage.elec_struc_mem, 1, MPI_DOUBLE, memory_usage, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_memory_usage.mme_mem, 1, MPI_DOUBLE, memory_usage, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_memory_usage.real_space_mem, 1, MPI_DOUBLE, memory_usage, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

          MPI_Allgatherv( &component_memory_usage.CME_mem, 1, MPI_DOUBLE, memory_usage, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

     }

     

     // Row 1  - Peak usage for characteristic functions calculation                                       [GB]
     // Row 2  - Peak usage for strain calculation                                                         [GB]
     // Row 3  - Peak usage for piezoelectric potential calculation                                        [GB]
     // Row 4  - Peak usage for strain and piezoelectric potential downsampling calculation (if enabled)   [GB]
     // Row 5  - Peak usage for convolutions calculation                                                   [GB]
     // Row 6  - Usage for eigensolver workspace                                                           [GB]
     // Row 7  - Peak usage for electronic structure calculation                                           [GB]
     // Row 8  - Peak usage for momentum matrix elements calculation                                       [GB]
     // Row 9  - Peak usage for real space calculation                                                     [GB]
     // Row 10  - Peak usage for Coulomb matrix element calculation                                        [GB]
     // Row 11  - Overall peak memory usage                                                                [GB]


     free( recvcounts );
     free( displs );
     free( memory_usage );
     free( peak_memory_usage );
     free( str1 );

}
