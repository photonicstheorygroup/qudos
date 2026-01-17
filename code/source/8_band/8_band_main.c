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


//////////////////////////
//// Global variables ////
//////////////////////////

// Descriptions given in 8_band.h

char *endptr, graph_folder[80];

int N_shape, isotropic, softening, include_strain, include_piezo, include_shape_elastic_constants, shape_elastic_constants_order, optimise, include_overlaps, include_localisation, include_CME, include_MME, PACKED, Num_Possible_Polarisations, Polarisation_Indicator;

double sigma_softening;

MKL_INT N_H;

struct characteristic_functions *characteristic_function_G;
struct characteristic_functions *characteristic_function_strain_G;
struct strain_tensor strain_G;
struct strain_tensor strain_region_G;
struct convolution_char_strain *convolutions_char_strain;
struct piezo_potential phi_piezo_G;
struct piezo_potential phi_piezo_region_G;
struct timing component_run_times;

complex double *dH_dk, *yV;



///////////////////////////////
//// MPI & BLACS variables ////
///////////////////////////////

// Descriptions given in 8_band.h

int myid, nprocs;
MKL_INT na_cols, na_rows, na_colsz, na_rowsz;
MKL_INT np_cols, np_rows;
MKL_INT my_prow, my_pcol;

MKL_INT my_blacs_ctxt, sc_desc[9], sc_descz[9], info;

MKL_INT block_size;


int main( int argc, char *argv[] )
{
	
	double start, start_1, finish, elapsed; // Variables to measure run time of the different calculation stages


	
	////////////////////////////////////
	//// MPI & BLACS initialisation ////
	////////////////////////////////////
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	double start_total = MPI_Wtime();

	// Optimally setting numbers of BLACS process rows and columns
	for (np_cols = (int) sqrt((double) nprocs); np_cols > 1; np_cols--) 
	{

			if (nprocs % np_cols == 0)
			{
	       		
				break;
	     	
			}

   	}
	np_rows = nprocs/np_cols;
	
	// Setting up BLACS process grid
	my_blacs_ctxt = MPI_Comm_c2f(MPI_COMM_WORLD);
	BLACS_GRIDINIT(&my_blacs_ctxt, "C", &np_rows, &np_cols);
	BLACS_GRIDINFO(&my_blacs_ctxt, &np_rows, &np_cols, &my_prow, &my_pcol);



	///////////////////////////////////////////////////////////////////////////////////////
	//// Read parameters.in file, command line parameters, and supercell.in input file ////
	///////////////////////////////////////////////////////////////////////////////////////

	if ( myid == 0 )
	{

		printf( "Reading supercell parameters and setting up plane wave basis set...\n");

	}

	// Parameters to be read from parameters.in / command line
	int N_Gx, N_Gy, N_Gz; 		// No. of plane waves for envelope functions
	int strain_region_type; 	// Free strain supercell lengths or commensurate strain supercell lengths
	int N_S_Gx, N_S_Gy, N_S_Gz; // No. of plane waves for strain supercell (corresponding to above for envelope functions of electronic structure supercell)
	int coulomb_region_type; 	// Free Coulomb supercell lengths or commensurate Coulomb supercell lengths
	int N_C_Gx, N_C_Gy, N_C_Gz; // No. of plane waves for Coulomb supercell (corresponding to above for envelope functions of electronic structure supercell)
	int N_vb, N_cb; 			// No. of VB states and CB states to be computed (excluding factor of 2 for spin degenaracy)
	int N_grid; 				// No. of real-space points for linescans
	int cross_sections; 		// Enable calculation of real-space cross-sections
	int N_grid_cross_sections; 	// No. of real-sapce points for cross-sections (same for both dimensions)
	int enable_3D; 				// Enable calculation of real-space 3D outputs
	int N_grid_3D; 				// No. of real-sapce points for 3D outputs (same for all three dimensions)
	double T; 					// Temperature
	double Fx, Fy, Fz; 			// Electric field component for each dimension
	double Kx, Ky, Kz;			// Supercell wave vector component for each direction (for miniband dispersion)
	
	PACKED = 0;
	
	int *directions = ( int* ) malloc( N_DIRECT*sizeof( int ) ); // Array of bools to indicate whether to compute linescans along each direction
	
	int *directions_cross_sections = ( int* ) malloc( N_DIRECT_CS*sizeof( int ) ); // Array of bools to indicate whether to compute cross-sections with normals along each direction

	int *input_check = ( int* ) calloc( N_PARAM, sizeof( int ) ); // Array to check that all input parameters have been provided
	
	// Read calculation parameters
	read_parameters( &N_shape, &N_Gx, &N_Gy, &N_Gz, &strain_region_type, &N_S_Gx, &N_S_Gy, &N_S_Gz, &N_vb, &N_cb, &T, &Fx, &Fy, &Fz, &Kx, &Ky, &Kz, &N_grid, &include_strain, &isotropic, &include_shape_elastic_constants, &shape_elastic_constants_order, &include_piezo, &optimise, &include_overlaps, &include_localisation, &include_MME, &include_CME, &coulomb_region_type, &N_C_Gx, &N_C_Gy, &N_C_Gz, &softening, &sigma_softening, directions, &cross_sections, &N_grid_cross_sections, directions_cross_sections, &enable_3D, &N_grid_3D, input_check );
	
	// Read command line inputs
	read_command_line_input( argc, argv, &N_shape, &N_Gx, &N_Gy, &N_Gz, &strain_region_type, &N_S_Gx, &N_S_Gy, &N_S_Gz, &N_vb, &N_cb, &T, &Fx, &Fy, &Fz, &Kx, &Ky, &Kz, &N_grid, &include_strain, &isotropic, &include_shape_elastic_constants, &shape_elastic_constants_order, &include_piezo, &optimise, &include_overlaps, &include_localisation, &include_MME, &include_CME, &coulomb_region_type, &N_C_Gx, &N_C_Gy, &N_C_Gz, &softening, &sigma_softening, directions, &cross_sections, &N_grid_cross_sections, directions_cross_sections, &enable_3D, &N_grid_3D, input_check );

	free( input_check );


	struct vector F; F.x = Fx; F.y = Fy; F.z = Fz; // Electric field
	
	struct supercell_geometry *supercell = ( _supercell_geometry* ) malloc( N_shape*sizeof( _supercell_geometry ) ); // To hold supercell and shape parameters

	// Read supercell parameters
	read_supercell( N_shape, strain_region_type, coulomb_region_type, supercell );


	// In the case of use of relative strain supercell lengths, check if the number of plane waves used for the strain supercell is above the minimum allowed.
	if ( include_strain == 1 && strain_region_type == 1 && myid == 0 )
	{

		int strain_plane_wave_error = 0;
		if ( ( N_S_Gx/2 ) < ( N_Gx/2 )*supercell[0].RL1/2 + ( ( N_Gx/2 ) % 2 )*( supercell[0].RL1 % 2 ) )
		{

			printf( "Error. M_S_Gx is too small. M_S_Gx must be greater than or equal to ( M_Gx times the relative length of the strain supercell -- strain supercell length divided by electronic structure supercell length -- along x  ) divided by 2. This minimum should be rounded up if not integer. It is recommended to let M_S_Gx equal to M_Gx times the relative length of the strain supercell along x.\n\n" );
			strain_plane_wave_error = 1;

		}
		if ( ( N_S_Gy/2 ) < ( N_Gy/2 )*supercell[0].RL2/2 + ( ( N_Gy/2 ) % 2 )*( supercell[0].RL2 % 2 ) )
		{

			printf( "Error. M_S_Gy is too small. M_S_Gy must be greater than or equal to ( M_Gy times the relative length of the strain supercell -- strain supercell length divided by electronic structure supercell length -- along y ) divided by 2. This minimum should be rounded up if not integer. It is recommended to let M_S_Gy equal to M_Gy times the relative length of the strain supercell along y.\n\n" );
			strain_plane_wave_error = 1;

		}
		if ( ( N_S_Gz/2 ) < ( N_Gz/2 )*supercell[0].RL3/2 + ( ( N_Gz/2 ) % 2 )*( supercell[0].RL3 % 2 ) )
		{

			printf( "Error. M_S_Gz is too small. M_S_Gz must be greater than or equal to ( M_Gz times the relative length of the strain supercell -- strain supercell length divided by electronic structure supercell length -- along z ) divided by 2. This minimum should be rounded up if not integer. It is recommended to let M_S_Gz equal to M_Gz times the relative length of the strain supercell along z.\n\n" );
			strain_plane_wave_error = 1;

		}

		if ( strain_plane_wave_error == 1 )
		{

			MPI_Abort( MPI_COMM_WORLD, 1 );

		}

	}

	// Check if strain supercell and electronic structure supercell are equally sized to allow to later skip downsampling step
	if ( strain_region_type == 1 && supercell[0].RL1 == 1 && supercell[0].RL2 == 1 && supercell[0].RL3 == 1 && N_Gx == N_S_Gx && N_Gy == N_S_Gy && N_Gz == N_S_Gz )
	{

		strain_region_type = 0;

	}

	

	/////////////////////////////////////
	//// Set up plane wave basis set ////
	/////////////////////////////////////
	
	// Electronic structure supercell basis set
	struct wave_vectors basis;

	basis.N_Gx = N_Gx;
	basis.N_Gy = N_Gy;
	basis.N_Gz = N_Gz;
	
	basis.N_Gx_T = 2*N_Gx - 1;
	basis.N_Gy_T = 2*N_Gy - 1;
	basis.N_Gz_T = 2*N_Gz - 1;
	
	basis.Gx = ( double* ) malloc( basis.N_Gx_T*sizeof( double ) );
	basis.Gy = ( double* ) malloc( basis.N_Gy_T*sizeof( double ) );
	basis.Gz = ( double* ) malloc( basis.N_Gz_T*sizeof( double ) );

	basis.Kx = Kx;
	basis.Ky = Ky;
	basis.Kz = Kz;

	// Compute basis set for electronic structure supercell
	basis_set( basis, supercell, 0, 0 );
	
	
	// Strain supercell basis set
	struct wave_vectors strain_basis;
	
	if ( include_strain == 1 )
	{
	
		strain_basis.N_Gx = N_S_Gx;
		strain_basis.N_Gy = N_S_Gy;
		strain_basis.N_Gz = N_S_Gz;

		strain_basis.N_Gx_T = 2*N_S_Gx - 1;
		strain_basis.N_Gy_T = 2*N_S_Gy - 1;
		strain_basis.N_Gz_T = 2*N_S_Gz - 1;

		strain_basis.Gx = ( double* ) malloc( strain_basis.N_Gx_T*sizeof( double ) );
		strain_basis.Gy = ( double* ) malloc( strain_basis.N_Gy_T*sizeof( double ) );
		strain_basis.Gz = ( double* ) malloc( strain_basis.N_Gz_T*sizeof( double ) );

		// Compute basis set for strain supercell
		basis_set( strain_basis, supercell, 1, 0 );
		
	}
	      
	if ( myid == 0 )
	{
	
		// Write electronic structure supercell and strain supercell basis set to file
		write_basis_set( basis, strain_basis, supercell );
	
	}



	////////////////////////////////////////////////////////////////////////////////
	//// Set number of local rows & columns for Hamiltonian matrix distribution ////
	////////////////////////////////////////////////////////////////////////////////
	
	N_H = (MKL_INT)N_BAND*basis.N_Gx*basis.N_Gy*basis.N_Gz; // Order of plane-wave Hamiltonian matrix
	
	MKL_INT f = 0; // NUMROC parameter
	MKL_INT m = 2*N_cb + 2*N_vb; // Total number of eigenstates to compute (including factor of 2 for spin degenaracy)

	if ( nprocs == 1 )
	{

		block_size = N_H; // Block size is matrix order for single-process OpenMP only case

	}
	else
	{

		block_size = 64; // Optimal block size is machine dependent. Can be modified here.

	}
	
	na_rows = NUMROC(&N_H, &block_size, &my_prow, &f, &np_rows); // Setting number of process-local Hamiltonian matrix rows
	na_cols = NUMROC(&N_H, &block_size, &my_pcol, &f, &np_cols); // Setting number of process-local Hamiltonian matrix columns
	na_rowsz = NUMROC(&N_H, &block_size, &my_prow, &f, &np_rows); // Setting number of process-local eigenvector matrix rows
	na_colsz = NUMROC(&m, &block_size, &my_pcol, &f, &np_cols); // Setting number of process-local eigenvector matrix columns
	
	if ( nprocs == 1 ) // In case of single process, OpenMP only configuration
	{
	
		na_rows = N_H;
		na_cols = N_H;
		na_rowsz = N_H;
		na_colsz = m;
	
	}
	
	DESCINIT(sc_desc, &N_H, &N_H, &block_size, &block_size, &f, &f, &my_blacs_ctxt, &na_rows, &info); // Initialising array descriptor for distributed Hamiltonian matrix
	DESCINIT(sc_descz, &N_H, &N_H, &block_size, &block_size, &f, &f, &my_blacs_ctxt, &na_rows, &info); // Initialising array descriptor for distributed eigenvector matrix
	


	////////////////////////////////////////////////////
	//// Interpolate material parameters for alloys ////
	////////////////////////////////////////////////////

	struct material_params *parameters = ( _material_params* ) malloc( N_shape*sizeof( _material_params ) ); // Holds material parameters of supercell and shapes

	// Computes material parameters of each shape in supercell
	material_parameters( N_shape, T,             supercell, parameters );

	// Writes supercell and shape parameters to file
	write_supercell    ( N_shape, T, Fx, Fy, Fz, supercell, parameters );
	


	// Create folder to store output files
	if ( myid == 0 )
	{
	
		time_t t = time( NULL );
		struct tm current_time;
		current_time = *localtime(&t);
		sprintf( graph_folder, "Output_files_%d_%02d_%02d_%02d_%02d_%02d", current_time.tm_year+1900, current_time.tm_mon+1, current_time.tm_mday, current_time.tm_hour, current_time.tm_min, current_time.tm_sec ); // Output files folder will be named using the current time
		char create_graph_folder[80] = "mkdir -p ";
		strcat( create_graph_folder, graph_folder );
		system( create_graph_folder );
		
	}


	if ( myid == 0 )
	{

		printf( "Done.\n\n");

	}


	
	//////////////////////////////////////////////////////////////////////
	//// Calculate shape characteristic functions in reciprocal space ////
	//////////////////////////////////////////////////////////////////////

	// Strain supercell 
	if ( include_strain == 1 )
	{

		if ( myid == 0 )
		{

			// If strain supercell is same size as electronic structure supercell (with same basis set size as should be the case)
			if ( basis.N_Gx == strain_basis.N_Gx && supercell[0].L1 == supercell[0].SL1 && basis.N_Gy == strain_basis.N_Gy && supercell[0].L2 == supercell[0].SL2 && basis.N_Gz == strain_basis.N_Gz && supercell[0].L3 == supercell[0].SL3 )
			{

				printf( "Calculating shape characteristic function Fourier transforms...\n");

			}
			else
			{

				printf( "Calculating shape characteristic function Fourier transforms for strain supercell...\n");

			}

		}
	
		// Stores shape characteristic functions for the strain supercell basis set
		characteristic_function_strain_G = ( _characteristic_functions* ) malloc( N_shape*sizeof( _characteristic_functions ) );

		int shape_index_start; // Indicates whether to calculate Fourier tranform of electronic structure supercell characteristic function inside the strain supercell, if using efficient FFT-based method (commensurate supercell lengths) to downsample strain
		if ( strain_region_type == 0 ) // Free strain supercell lengths
		{

			shape_index_start = 1;

		}
		else // Commensurate strain supercell lengths with electronic structure supercell lengths
		{

			shape_index_start = 0;

		}

		for( int idx_shape = shape_index_start; idx_shape < N_shape; idx_shape++ )
		{

			characteristic_function_strain_G[idx_shape].chi = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );

		}

		start = MPI_Wtime();

		// Compute shape charcteristic functions in reciprocal space for strain supercell
		characteristic_functions_reciprocal_space( N_shape, strain_basis, supercell, strain_region_type + 1 );

		finish = MPI_Wtime();
		elapsed = finish - start;
		component_run_times.char_s_grid = elapsed;
		
		
		if ( myid == 0 )
		{

			printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );

		}
	
	}
	
	// Electronic structure supercell
	characteristic_function_G = ( _characteristic_functions* ) malloc( N_shape*sizeof( _characteristic_functions ) );

	for( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
	{
		
		characteristic_function_G[idx_shape].chi = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		
	}

	// Simple copy if strain supercell size is equal to electronic structure supercell size
	if ( include_strain == 1 && basis.N_Gx == strain_basis.N_Gx && supercell[0].L1 == supercell[0].SL1 && basis.N_Gy == strain_basis.N_Gy && supercell[0].L2 == supercell[0].SL2 && basis.N_Gz == strain_basis.N_Gz && supercell[0].L3 == supercell[0].SL3 )
	{
	
		for( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			for ( int idx_basis = 0; idx_basis < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T; idx_basis++ )
			{
			
				characteristic_function_G[idx_shape].chi[idx_basis] = characteristic_function_strain_G[idx_shape].chi[idx_basis];
			
			}

		}
	
	}
	else
	{

		if ( myid == 0 )
		{

			printf( "Calculating shape characteristic function Fourier transforms for electronic structure supercell...\n");

		}

		start = MPI_Wtime();

		// Compute shape charcteristic functions in reciprocal space for electronic structure supercell
		characteristic_functions_reciprocal_space( N_shape, basis, supercell, 0 );

		finish = MPI_Wtime();
		elapsed = finish - start;
		
		if ( myid == 0 )
		{
		
			printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );
			
		}

	}

	component_run_times.char_es_grid = elapsed;
	
	

	//////////////////////////////////////////////////////////////////////////
	//// Calculate strain tensor in reciprocal space for strain supercell ////
	//////////////////////////////////////////////////////////////////////////
	
	if ( include_strain == 1)
	{

		if ( myid == 0 )
		{

			printf( "Calculating strain in reciprocal space for strain supercell...\n");

		}

		strain_region_G.xx = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
		strain_region_G.yy = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
		strain_region_G.zz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
		strain_region_G.yz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
		strain_region_G.xz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
		strain_region_G.xy = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );

		start = MPI_Wtime();

		// Compute strain tensor components in reciprocal space for strain supercell
		strain_region_reciprocal_space( N_shape, strain_basis, supercell, parameters );

		finish = MPI_Wtime();
		elapsed = finish - start;
		component_run_times.strain_s_grid = elapsed;

		if ( myid == 0 )
		{

			printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );

		}

	}



	////////////////////////////////////////////////////////////////////////////////////
	//// Calculate piezoelectric potential in reciprocal space for strain supercell ////
	////////////////////////////////////////////////////////////////////////////////////

	if ( include_strain == 1)
	{

		if( include_piezo != 0 )
		{

			if ( myid == 0 )
			{

				printf( "Calculating piezoelectric potential in reciprocal space for strain supercell...\n");

			}

			phi_piezo_region_G.pz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );

			start = MPI_Wtime();

			// Compute piezoelectric potential in reciprocal space for strain supercell
			piezo_region_reciprocal_space( N_shape, strain_basis, supercell, parameters );

			finish = MPI_Wtime();
			elapsed = finish - start;
			component_run_times.piezo_s_grid = elapsed;

			if ( myid == 0 )
			{

				printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );

			}

		}

		// Freeing memory used for shape characteristic functions in strain supercell as no longer needed
		for( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			free( characteristic_function_strain_G[idx_shape].chi );

		}

	}
		


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Calculate strain tensor and piezoelectric potential in reciprocal space for electronic structure supercell ////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if ( include_strain == 1)
	{
		
		// Downsampling strain tensor components for electronic structure supercell

		if ( myid == 0 )
		{

			printf( "Downsampling strain for electronic structure supercell...\n");

		}
		
		strain_G.xx = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_G.yy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_G.zz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_G.yz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_G.xz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		strain_G.xy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );

		start = MPI_Wtime();

		if ( strain_region_type == 0 ) // If using free strain supercell lengths
		{

			strain_reciprocal_space( basis, strain_basis, supercell );

		}
		else // If using commensurate strain and electronic structure supercell lengths
		{

			strain_reciprocal_space_fft( basis, strain_basis, supercell );

		}
		
		finish = MPI_Wtime();
    	elapsed = finish - start;
		component_run_times.strain_es_grid = elapsed;
		
		if ( myid == 0 )
		{

	    	printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );
	    	
		}

		free( strain_region_G.xx );
		free( strain_region_G.yy );
		free( strain_region_G.zz );
		free( strain_region_G.yz );
		free( strain_region_G.xz );
		free( strain_region_G.xy );


		// Downsampling piezoelectric potential for electronic structure supercell
		if( include_piezo != 0 )
		{

			if ( myid == 0 )
			{

				printf( "Downsampling piezoelectric potential for electronic structure supercell...\n");

			}

			phi_piezo_G.pz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );

			start = MPI_Wtime();

			if ( strain_region_type == 0 ) // If using free strain supercell lengths
			{

				piezo_reciprocal_space( basis, strain_basis, supercell );

			}
			else // If using commensurate strain and electronic structure supercell lengths
			{

				piezo_reciprocal_space_fft( basis, strain_basis, supercell );

			}

			finish = MPI_Wtime();
    		elapsed = finish - start;
			component_run_times.piezo_es_grid = elapsed;

			if ( myid == 0)
			{

				printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );

			}

			free( phi_piezo_region_G.pz );

		}

		if ( strain_region_type == 1 ) // Freeing Fourier transform of electronic structure supercell characteristic function in strain supercell
		{

			free( characteristic_function_strain_G[0].chi );

		}
		free( characteristic_function_strain_G );
		
	}
	

	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Calculate convolution of shape characteristic functions with strain tensor compoenents for Hamiltonian matrix element calculation ////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	if ( include_strain == 1)
	{

		if ( myid == 0 )
		{

			printf( "Calculating convolutions between strain tensor components and shape characteristic functions...\n");

		}
		
		// Will store convolutions between shape characteristic functions and strain tensor components
		convolutions_char_strain = ( _convolution_char_strain* ) malloc( N_shape*sizeof( _convolution_char_strain ) );

		for ( int idx_shape = 0; idx_shape < N_shape; idx_shape++ )
		{
		
			convolutions_char_strain[idx_shape].char_xx = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			convolutions_char_strain[idx_shape].char_yy = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			convolutions_char_strain[idx_shape].char_zz = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			convolutions_char_strain[idx_shape].char_yz = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			convolutions_char_strain[idx_shape].char_xz = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			convolutions_char_strain[idx_shape].char_xy = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			convolutions_char_strain[idx_shape].char_hydrostatic = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			convolutions_char_strain[idx_shape].char_biaxial = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			
		}
		
		start = MPI_Wtime();

		// Compute convolutions
		convolutions( N_shape, basis, supercell );
		
		finish = MPI_Wtime();
    	elapsed = finish - start;
		component_run_times.conv = elapsed;
		
		if ( myid == 0 )
		{

	    	printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );
	    	
		}
		
	}
	

	
	////////////////////////////////////////
	//// Calculate electronic structure ////
	////////////////////////////////////////


	// Allocate memory to store eigenstates
	
	int N_states[2]; N_states[0] = N_vb; N_states[1] = N_cb;
	
	struct eigenstate *eigenstate_vb = malloc( 2*N_states[0]*sizeof( _eigenstate ) ); // Factor of 2 to account for spin degeneracy
	struct eigenstate *eigenstate_cb = malloc( 2*N_states[1]*sizeof( _eigenstate ) ); // Factor of 2 to account for spin degeneracy
	
	for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
	{

			eigenstate_vb[idx_state].fourier       = calloc( na_rows     , sizeof( complex double ) ); // Fourier coefficients
			eigenstate_vb[idx_state].bloch         = calloc( N_BAND  , sizeof(         double ) ); // Bloch characters

	}

	for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
	{

		eigenstate_cb[idx_state].fourier       = calloc( na_rows     , sizeof( complex double ) ); // Fourier coefficients
		eigenstate_cb[idx_state].bloch         = calloc( N_BAND  , sizeof(         double ) ); // Bloch characters
		
	}

	if ( myid == 0 )
	{

		printf( "Calculating electronic structure...\n");

	}
	
	start = MPI_Wtime();
	
	// Calculate eigenstates
	eigenstates( N_shape, N_states, F, basis, supercell, parameters, eigenstate_vb, eigenstate_cb );
	
	finish = MPI_Wtime();
    elapsed = finish - start;
	component_run_times.eig = elapsed;
	
	if ( myid == 0 )
	{
	
		printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );
		
	}
	
	// Compute eigenstate Bloch character
	if ( myid == 0 )
	{

		printf( "Calculating eigenstate Bloch character...\n");

	}

	start = MPI_Wtime();
	
	bloch_character( 0, N_states, basis, eigenstate_vb );
	bloch_character( 1, N_states, basis, eigenstate_cb );

	finish = MPI_Wtime();
    elapsed = finish - start;
	component_run_times.bloch = elapsed;

	if ( myid == 0 )
	{
	
		printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );
		
	}
	
	// Write computed state energies to file
	if ( myid == 0)
	{

	    write_energies(N_states, eigenstate_vb, eigenstate_cb);

	}



	////////////////////////////////////////////////////////////////////////////////////////////
	//// Calculate and write real space envelope function and probability density linescans ////
	////////////////////////////////////////////////////////////////////////////////////////////

	if ( myid == 0 )
	{

		printf( "Calculating real space linescans...\n");

	}

	start_1 = MPI_Wtime();
	
	// Direction indicator
	int direction;

	// Real space grid
	double *x = ( double* ) malloc( N_grid*sizeof( double ) );
	double *y = ( double* ) malloc( N_grid*sizeof( double ) );
	double *z = ( double* ) malloc( N_grid*sizeof( double ) );
	
	if ( myid == 0 )
	{
	
		for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
		{

			eigenstate_vb[idx_state].probability_1 = calloc( N_grid, sizeof(         double ) );
			eigenstate_vb[idx_state].probability_2 = calloc( N_grid, sizeof(         double ) );
			eigenstate_vb[idx_state].probability_3 = calloc( N_grid, sizeof(         double ) );
			eigenstate_vb[idx_state].probability_4 = calloc( N_grid, sizeof(         double ) );
			eigenstate_vb[idx_state].probability_5 = calloc( N_grid, sizeof(         double ) );
			eigenstate_vb[idx_state].probability_6 = calloc( N_grid, sizeof(         double ) );
			eigenstate_vb[idx_state].probability_7 = calloc( N_grid, sizeof(         double ) );
			eigenstate_vb[idx_state].probability_8 = calloc( N_grid, sizeof(         double ) );
			eigenstate_vb[idx_state].F_1           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_vb[idx_state].F_2           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_vb[idx_state].F_3           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_vb[idx_state].F_4           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_vb[idx_state].F_5           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_vb[idx_state].F_6           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_vb[idx_state].F_7           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_vb[idx_state].F_8           = calloc( N_grid, sizeof( complex double ) );

		}

		for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
		{

			eigenstate_cb[idx_state].probability_1 = calloc( N_grid, sizeof(         double ) );
			eigenstate_cb[idx_state].probability_2 = calloc( N_grid, sizeof(         double ) );
			eigenstate_cb[idx_state].probability_3 = calloc( N_grid, sizeof(         double ) );
			eigenstate_cb[idx_state].probability_4 = calloc( N_grid, sizeof(         double ) );
			eigenstate_cb[idx_state].probability_5 = calloc( N_grid, sizeof(         double ) );
			eigenstate_cb[idx_state].probability_6 = calloc( N_grid, sizeof(         double ) );
			eigenstate_cb[idx_state].probability_7 = calloc( N_grid, sizeof(         double ) );
			eigenstate_cb[idx_state].probability_8 = calloc( N_grid, sizeof(         double ) );
			eigenstate_cb[idx_state].F_1           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_cb[idx_state].F_2           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_cb[idx_state].F_3           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_cb[idx_state].F_4           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_cb[idx_state].F_5           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_cb[idx_state].F_6           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_cb[idx_state].F_7           = calloc( N_grid, sizeof( complex double ) );
			eigenstate_cb[idx_state].F_8           = calloc( N_grid, sizeof( complex double ) );

		}
	
	}

	for( direction = 0; direction < N_DIRECT; direction++ ) // Iterating directions
	{

		if( directions[direction] == 1 ) // Compute real space linescans along current direction
		{

			start = MPI_Wtime();

			// Set real space grid for current linescan direction
			real_space_grid_linescan( direction, N_grid, supercell, x, y, z );

			// Envelope function linescans			
			envelope_functions ( supercell, N_states, 0, N_grid, x, y, z, basis, eigenstate_vb, direction ); // Valence    band envelope functions
			envelope_functions ( supercell, N_states, 1, N_grid, x, y, z, basis, eigenstate_cb, direction ); // Conduction band envelope functions

			finish = MPI_Wtime();
    		elapsed = finish - start;
			
			if (myid == 0)
			{
			
				write_envelope_functions ( N_states, N_grid, x, y, z, basis, eigenstate_vb, eigenstate_cb, direction );
				write_probability_density( N_states, N_grid, x, y, z, basis, eigenstate_vb, eigenstate_cb, direction );
			
			}
			
			if ( myid == 0 )
			{
			
				printf( "\tTime to compute for direction %d = %9e seconds.\n", direction, elapsed );
				
			}

		}

	}

	elapsed = finish - start_1;
	component_run_times.linescan = elapsed;

	if ( myid == 0 )
	{
	
		printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );
		
	}
	
	if ( myid == 0 )
	{
	
		for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
		{

			free( eigenstate_vb[idx_state].probability_1 );
			free( eigenstate_vb[idx_state].probability_2 );
			free( eigenstate_vb[idx_state].probability_3 );
			free( eigenstate_vb[idx_state].probability_4 );
			free( eigenstate_vb[idx_state].probability_5 );
			free( eigenstate_vb[idx_state].probability_6 );
			free( eigenstate_vb[idx_state].probability_7 );
			free( eigenstate_vb[idx_state].probability_8 );
			free( eigenstate_vb[idx_state].F_1 );
			free( eigenstate_vb[idx_state].F_2 );
			free( eigenstate_vb[idx_state].F_3 );
			free( eigenstate_vb[idx_state].F_4 );
			free( eigenstate_vb[idx_state].F_5 );
			free( eigenstate_vb[idx_state].F_6 );
			free( eigenstate_vb[idx_state].F_7 );
			free( eigenstate_vb[idx_state].F_8 );

		}

		for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
		{

			free( eigenstate_cb[idx_state].probability_1 );
			free( eigenstate_cb[idx_state].probability_2 );
			free( eigenstate_cb[idx_state].probability_3 );
			free( eigenstate_cb[idx_state].probability_4 );
			free( eigenstate_cb[idx_state].probability_5 );
			free( eigenstate_cb[idx_state].probability_6 );
			free( eigenstate_cb[idx_state].probability_7 );
			free( eigenstate_cb[idx_state].probability_8 );
			free( eigenstate_cb[idx_state].F_1 );
			free( eigenstate_cb[idx_state].F_2 );
			free( eigenstate_cb[idx_state].F_3 );
			free( eigenstate_cb[idx_state].F_4 );
			free( eigenstate_cb[idx_state].F_5 );
			free( eigenstate_cb[idx_state].F_6 );
			free( eigenstate_cb[idx_state].F_7 );
			free( eigenstate_cb[idx_state].F_8 );

		}
        
    }

	free( x );
	free( y );
	free( z );
        

        
    /////////////////////////////////////////////////////////////////////////////////////////////////
	//// Calculate and write real space envelope function and probability density cross-sections ////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	if ( cross_sections )
	{

		if ( myid == 0 )
		{

			printf( "Calculating real space cross-sections...\n");

		}

		start_1 = MPI_Wtime();
	
		int direction; // Direction of normal
		
		int N_grid_cross_sections_squared = N_grid_cross_sections * N_grid_cross_sections; // Number of 2D real-space points

		double *x = ( double* ) malloc( N_grid_cross_sections_squared*sizeof( double ) );
		double *y = ( double* ) malloc( N_grid_cross_sections_squared*sizeof( double ) );
		double *z = ( double* ) malloc( N_grid_cross_sections_squared*sizeof( double ) );
		
		if ( myid == 0 )
		{
		
			for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
			{

				eigenstate_vb[idx_state].probability_1 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_2 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_3 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_4 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_5 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_6 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_7 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_8 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_vb[idx_state].F_1           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_2           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_3           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_4           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_5           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_6           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_7           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_8           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );

			}

			for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
			{

		      	eigenstate_cb[idx_state].probability_1 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_2 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_3 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_4 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_5 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_6 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_7 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_8 = calloc( N_grid_cross_sections_squared, sizeof(         double ) );
				eigenstate_cb[idx_state].F_1           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_2           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_3           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_4           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_5           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_6           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_7           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_8           = calloc( N_grid_cross_sections_squared, sizeof( complex double ) );

			}
		
		}
		
		for( direction = 0; direction < N_DIRECT_CS; direction++ )
		{

			if( directions_cross_sections[direction] == 1 ) // Compute real space cross-sections for normal along current direction
			{

				start = MPI_Wtime();

				// Set real space grid for current cross-section normal direction
				real_space_grid_cross_section( direction, N_grid_cross_sections, supercell, x, y, z );


				// Envelope function cross-sections
				envelope_functions ( supercell, N_states, 0, N_grid_cross_sections_squared, x, y, z, basis, eigenstate_vb, direction ); // Valence    band envelope functions
		   		envelope_functions ( supercell, N_states, 1, N_grid_cross_sections_squared, x, y, z, basis, eigenstate_cb, direction ); // Conduction band envelope functions

				finish = MPI_Wtime();
    			elapsed = finish - start;
		   		
		   		if (myid == 0)
		   		{
		   		
		 			write_envelope_functions_cross_section ( N_states, N_grid_cross_sections, x, y, z, basis, eigenstate_vb, eigenstate_cb, direction );
					write_probability_density_cross_section ( N_states, N_grid_cross_sections, x, y, z, basis, eigenstate_vb, eigenstate_cb, direction );
		 			
		 		}
				
				if ( myid == 0 )
				{
				
					printf( "\tTime to compute for normal along direction %d = %9e seconds.\n", direction, elapsed );
					
				}

			}


		}
		
		free( x );
		free( y );
		free( z );

		elapsed = finish - start_1;
		component_run_times.cross_section = elapsed;

		if ( myid == 0 )
		{
		
			printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );
			
		}
		
		if( myid == 0 )
		{
		
			for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
			{

			       free( eigenstate_vb[idx_state].probability_1 );
			       free( eigenstate_vb[idx_state].probability_2 );
			       free( eigenstate_vb[idx_state].probability_3 );
			       free( eigenstate_vb[idx_state].probability_4 );
			       free( eigenstate_vb[idx_state].probability_5 );
			       free( eigenstate_vb[idx_state].probability_6 );
			       free( eigenstate_vb[idx_state].probability_7 );
			       free( eigenstate_vb[idx_state].probability_8 );
			       free( eigenstate_vb[idx_state].F_1 );
			       free( eigenstate_vb[idx_state].F_2 );
			       free( eigenstate_vb[idx_state].F_3 );
			       free( eigenstate_vb[idx_state].F_4 );
			       free( eigenstate_vb[idx_state].F_5 );
			       free( eigenstate_vb[idx_state].F_6 );
			       free( eigenstate_vb[idx_state].F_7 );
			       free( eigenstate_vb[idx_state].F_8 );

			}

			for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
			{

			       free( eigenstate_cb[idx_state].probability_1 );
			       free( eigenstate_cb[idx_state].probability_2 );
			       free( eigenstate_cb[idx_state].probability_3 );
			       free( eigenstate_cb[idx_state].probability_4 );
			       free( eigenstate_cb[idx_state].probability_5 );
			       free( eigenstate_cb[idx_state].probability_6 );
			       free( eigenstate_cb[idx_state].probability_7 );
			       free( eigenstate_cb[idx_state].probability_8 );
			       free( eigenstate_cb[idx_state].F_1 );
			       free( eigenstate_cb[idx_state].F_2 );
			       free( eigenstate_cb[idx_state].F_3 );
			       free( eigenstate_cb[idx_state].F_4 );
			       free( eigenstate_cb[idx_state].F_5 );
			       free( eigenstate_cb[idx_state].F_6 );
			       free( eigenstate_cb[idx_state].F_7 );
			       free( eigenstate_cb[idx_state].F_8 );

			}
		
		}
        
    }



	////////////////////////////////////////////////////////////////////////////////////////
	//// Calculate and write 3D real space envelope functions and probability densities ////
	////////////////////////////////////////////////////////////////////////////////////////
	
	if ( enable_3D )
	{
		
		int N_grid_3D_cubed = N_grid_3D * N_grid_3D * N_grid_3D;

		double *x = ( double* ) malloc( N_grid_3D_cubed*sizeof( double ) );
		double *y = ( double* ) malloc( N_grid_3D_cubed*sizeof( double ) );
		double *z = ( double* ) malloc( N_grid_3D_cubed*sizeof( double ) );

		if ( myid == 0 )
		{

			printf( "Calculating real space 3D outputs...\n");

		}
		
		if ( myid == 0 )
		{
		
			for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
			{

				eigenstate_vb[idx_state].probability_1 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_2 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_3 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_4 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_5 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_6 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_7 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_vb[idx_state].probability_8 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_vb[idx_state].F_1           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_2           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_3           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_4           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_5           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_6           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_7           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_vb[idx_state].F_8           = calloc( N_grid_3D_cubed, sizeof( complex double ) );

			}

			for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
			{

		      	eigenstate_cb[idx_state].probability_1 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_2 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_3 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_4 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_5 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_6 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_7 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_cb[idx_state].probability_8 = calloc( N_grid_3D_cubed, sizeof(         double ) );
				eigenstate_cb[idx_state].F_1           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_2           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_3           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_4           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_5           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_6           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_7           = calloc( N_grid_3D_cubed, sizeof( complex double ) );
				eigenstate_cb[idx_state].F_8           = calloc( N_grid_3D_cubed, sizeof( complex double ) );

			}
		
		}

		start = MPI_Wtime();

		// Set 3D real space grid
		real_space_grid_3D( N_grid_3D, supercell, x, y, z );


		// 3D envelope function
		envelope_functions ( supercell, N_states, 0, N_grid_3D_cubed, x, y, z, basis, eigenstate_vb, direction ); // Valence    band envelope functions
		envelope_functions ( supercell, N_states, 1, N_grid_3D_cubed, x, y, z, basis, eigenstate_cb, direction ); // Conduction band envelope functions

		finish = MPI_Wtime();
		elapsed = finish - start;
		component_run_times.real_3D = elapsed;
		
		// Write envelope functions and probability densities
		if (myid == 0)
		{
		
			write_envelope_functions_3D( N_states, N_grid_3D, x, y, z, basis, eigenstate_vb, eigenstate_cb );
			write_probability_density_3D( N_states, N_grid_3D, x, y, z, basis, eigenstate_vb, eigenstate_cb );
			
		}
		
		if ( myid == 0 )
		{
		
			printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );
			
		}
		
		free( x );
		free( y );
		free( z );
		
		if( myid == 0 )
		{
		
			for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
			{

			       free( eigenstate_vb[idx_state].probability_1 );
			       free( eigenstate_vb[idx_state].probability_2 );
			       free( eigenstate_vb[idx_state].probability_3 );
			       free( eigenstate_vb[idx_state].probability_4 );
			       free( eigenstate_vb[idx_state].probability_5 );
			       free( eigenstate_vb[idx_state].probability_6 );
			       free( eigenstate_vb[idx_state].probability_7 );
			       free( eigenstate_vb[idx_state].probability_8 );
			       free( eigenstate_vb[idx_state].F_1 );
			       free( eigenstate_vb[idx_state].F_2 );
			       free( eigenstate_vb[idx_state].F_3 );
			       free( eigenstate_vb[idx_state].F_4 );
			       free( eigenstate_vb[idx_state].F_5 );
			       free( eigenstate_vb[idx_state].F_6 );
			       free( eigenstate_vb[idx_state].F_7 );
			       free( eigenstate_vb[idx_state].F_8 );

			}

			for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
			{

			       free( eigenstate_cb[idx_state].probability_1 );
			       free( eigenstate_cb[idx_state].probability_2 );
			       free( eigenstate_cb[idx_state].probability_3 );
			       free( eigenstate_cb[idx_state].probability_4 );
			       free( eigenstate_cb[idx_state].probability_5 );
			       free( eigenstate_cb[idx_state].probability_6 );
			       free( eigenstate_cb[idx_state].probability_7 );
			       free( eigenstate_cb[idx_state].probability_8 );
			       free( eigenstate_cb[idx_state].F_1 );
			       free( eigenstate_cb[idx_state].F_2 );
			       free( eigenstate_cb[idx_state].F_3 );
			       free( eigenstate_cb[idx_state].F_4 );
			       free( eigenstate_cb[idx_state].F_5 );
			       free( eigenstate_cb[idx_state].F_6 );
			       free( eigenstate_cb[idx_state].F_7 );
			       free( eigenstate_cb[idx_state].F_8 );

			}
		
		}
        
    }


	// Move real-space envelope function and probability density outputs to output folder

	for( direction = 0; direction < N_DIRECT; direction++ ) // Iterating directions
	{

		if( directions[direction] == 1 ) // Compute real space linescans along current direction
		{

			
			if (myid == 0)
			{

				copy_output( 0, direction, 2*N_states[0], 2*N_states[1] );
			
			}

		}

	}

	if ( cross_sections )
	{

		for( direction = 0; direction < N_DIRECT_CS; direction++ )
		{

			if( directions_cross_sections[direction] == 1 ) // Compute real space cross-sections for normal along current direction
			{
				
				if (myid == 0)
				{

					copy_output( 1, direction, 2*N_states[0], 2*N_states[1] );
					
				}

			}

		}

	}

	if ( enable_3D )
	{

		if (myid == 0)
		{

			copy_output( 2, 0, 2*N_states[0], 2*N_states[1] );
			
		}

	}



	//////////////////////////////////////////////////////////////
	//// Calculate conduction and valence band state overlaps ////
	//////////////////////////////////////////////////////////////

	if ( include_overlaps )
	{

		if ( myid == 0 )
		{

			printf( "Calculating conduction and valence band state overlaps...\n");

		}

		start = MPI_Wtime();

		complex double *overlaps = ( complex double* ) malloc( 4*N_states[0]*N_states[1]*sizeof(complex double) ); // Will store overlaps

		for( int idx_state_cb = 0; idx_state_cb < 2*N_states[1]; idx_state_cb++ ) // Loop over CB states
		{

			for( int idx_state_vb = 0; idx_state_vb < 2*N_states[0]; idx_state_vb++ ) // Loop over VB states
			{

				overlaps[idx_state_cb * 2*N_states[0] + idx_state_vb] = spatial_overlap( basis, eigenstate_cb[idx_state_cb], 2*N_states[0] + idx_state_cb, eigenstate_vb[2*N_states[0] - 1 - idx_state_vb], 2*N_states[0] - 1 - idx_state_vb );
				
			}               

		}

		if (myid == 0) // Write overlaps to file
		{

			write_overlaps( N_states, overlaps );
			
		}

		finish = MPI_Wtime();
		elapsed = finish - start;
		component_run_times.overlaps = elapsed;

		if ( myid == 0 )
		{

			printf( "Done. Time to compute = %13e seconds.\n\n", elapsed );
			
		}
	
	}
	
	

	//////////////////////////////
	//// Carrier localisation ////
	//////////////////////////////

	if ( include_localisation )
	{

		if ( myid == 0 )
		{

			printf( "Calculating carrier localisation...\n");

		}
		
		start = MPI_Wtime();

		int idx_shape;
		complex double sum;


		// Will store local components of FFTs for shape characteristic functions
		struct characteristic_functions *characteristic_function_G_local = ( _characteristic_functions* ) malloc( N_shape*sizeof( _characteristic_functions ) );
		// Shape characteristic function FFTs will only be stored once, the following will indicate whether it has already been computed for each shape
		int *shape_indicator = ( int* ) calloc( N_shape, sizeof( int ) );
		


		// Hole states
		complex double *localisation_vb = calloc( 2*N_states[0]*N_shape, sizeof( complex double ) );

		for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
		{

			sum = 0.0 + I*0.0;

			// Shapes
			for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
			{

				localisation_vb[idx_shape + N_shape*idx_state] = carrier_localisation( idx_shape, basis, supercell, N_states, characteristic_function_G_local, shape_indicator, 0, eigenstate_vb, idx_state );
				// To later calculate proportion in rest of supercell (outside of any shape)
				if ( supercell[idx_shape].embedded == 0 )
				{

					sum += localisation_vb[idx_shape + N_shape*idx_state];

				}

			}

			// To later calculate proportion in rest of supercell (outside of any shape)
			idx_shape = 0;
			localisation_vb[idx_shape + N_shape*idx_state] = sum;

		}


		// Electron states
		complex double *localisation_cb = calloc( 2*N_states[1]*N_shape, sizeof( complex double ) );

		for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ ) // Valence states
		{

			sum = 0.0 + I*0.0;

			// Shapes
			for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
			{

				localisation_cb[idx_shape + N_shape*idx_state] = carrier_localisation( idx_shape, basis, supercell, N_states, characteristic_function_G_local, shape_indicator, 1, eigenstate_cb, idx_state );
				// To later calculate proportion in rest of supercell (outside of any shape)
				if ( supercell[idx_shape].embedded == 0 )
				{

					sum += localisation_cb[idx_shape + N_shape*idx_state];

				}

			}

			// To later calculate proportion in rest of supercell (outside of any shape)
			idx_shape = 0;
			localisation_cb[idx_shape + N_shape*idx_state] = sum;

		}


		for ( int shape_idx = 1; shape_idx < N_shape; shape_idx++ )
		{

			if ( shape_indicator[shape_idx] == 1 )
			{

				free( characteristic_function_G_local[shape_idx].chi );

			}

		}

		free( characteristic_function_G_local );
		free( shape_indicator );


		// Summing process local contributions to the global value
		MPI_Allreduce( MPI_IN_PLACE, localisation_vb, 2*N_states[0]*N_shape, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD );
		MPI_Allreduce( MPI_IN_PLACE, localisation_cb, 2*N_states[1]*N_shape, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD );


		// Calculating proportions in rest of supercell (outside of any shape)
		for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ ) // Valence states
		{

			localisation_vb[N_shape*idx_state] = 1.0 - localisation_vb[N_shape*idx_state];

		}
		for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ ) // Conduction band states
		{

			localisation_cb[N_shape*idx_state] = 1.0 - localisation_cb[N_shape*idx_state];

		}


		if ( myid == 0 )
		{
		
			write_carrier_localisation( N_shape, N_states, localisation_vb, localisation_cb );
			
		}
		
		finish = MPI_Wtime();
		elapsed = finish - start;
		component_run_times.car_loc = elapsed;
		
		if ( myid == 0 )
		{
		
			printf( "Done. Time to compute = %13e seconds.\n\n", elapsed );
				
		}

		free( localisation_vb );
		free( localisation_cb );

	}
	
	
	
	//////////////////////////////////////////////////////
	//// Calculate interband momentum matrix elements ////
	//////////////////////////////////////////////////////

	if ( include_MME )
	{
		
		if ( myid == 0 )
		{

			printf( "Calculating momentum matrix elements between states...\n");

		}
		
		start = MPI_Wtime();
		
		int Num_CB_States = N_states[1]; // Not including degeneracies
		int Num_VB_States = N_states[0]; // Not including degeneracies
		
		Num_Possible_Polarisations = 3; // Calculating separately for [100], [010], and [001]
		int polarisations[Num_Possible_Polarisations];
		polarisations[0] = 1; // Enabling [100]
		polarisations[1] = 1; // Enabling [010]
		polarisations[2] = 1; // Enabling [001]
		int Num_Polarisations = 0;
		for ( int i = 0; i < Num_Possible_Polarisations; i++ )
		{
		
			if ( polarisations[i] == 1 )
			{

				Num_Polarisations++;
				
			}
			
		}
		
		// Will store momentum matrix elements
		complex double *momentum_matrix_elements = ( complex double* ) malloc( Num_Polarisations*2*Num_CB_States*2*Num_VB_States*sizeof(complex double) );
		
		// Reserving memory to store derivative of Hamiltonian matrix
		if (PACKED && nprocs == 1) // In case of 1 MPI process and packed matrix storage
		{
		
			dH_dk = ( complex double* ) calloc( N_H*(N_H + 1)/2, sizeof( complex double ) );
		
		}
		else // Otherwise
		{
		
			dH_dk = ( complex double* ) malloc( na_rows*na_cols*sizeof( complex double ) );
			
		}
		
		// Will store matrix-vector product part of vector-matrix-vector product to avoid unnecessary recomputation when computing momentum matrix elements involving the same state
		yV = (complex double *) calloc( na_rowsz, sizeof( complex double ) );
		
		// Indicates whether the matrix representation of the derivative of the Hamiltonian matrix has been computed for the current polarisation
		Polarisation_Indicator = -1;
		
		
		for ( int polarisation = 0; polarisation < Num_Possible_Polarisations; polarisation++ ) // Iterates over polarisations
		{
			
			if (polarisations[polarisation] == 1)
			{
			
				for( int idx_state_cb = 0; idx_state_cb < 2*Num_CB_States; idx_state_cb++ ) // Iterates over CB states
				{

					for( int idx_state_vb = 0; idx_state_vb < 2*Num_VB_States; idx_state_vb++ ) // Iterates over VB states
					{

						momentum_matrix_elements[(polarisation*2*Num_CB_States*2*Num_VB_States) + (idx_state_cb*2*Num_VB_States) + idx_state_vb] = momentum_matrix_element( polarisation, N_shape, basis, supercell, parameters, eigenstate_cb[idx_state_cb], 2*N_states[0] + idx_state_cb, eigenstate_vb[2*N_states[0] - 1 - idx_state_vb], 2*N_states[0] - 1 - idx_state_vb, N_states );
						
					}               

				}
				
			}
			
		}

		finish = MPI_Wtime();
		elapsed = finish - start;
		component_run_times.mme = elapsed;
		
		// Write momentum matrix elements
		if ( myid == 0 )
		{
		
			write_momentum_matrix_elements( Num_CB_States, Num_VB_States, polarisations, momentum_matrix_elements );
		
		}

		if ( myid == 0 )
		{
		
			printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );
			
		}
		
		free( dH_dk );
		free( yV );
		free( momentum_matrix_elements );

	}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Free dynamically allocated memory that is no longer required prior to Coulomb matrix elements calculation ////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	for( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
	{

		free( characteristic_function_G[idx_shape].chi );

	}

	free( characteristic_function_G );
	
	if (include_strain == 1)
	{
	
		free( strain_basis.Gx   );
		free( strain_basis.Gy   );
		free( strain_basis.Gz   );
	
		free( strain_G.xx );
		free( strain_G.yy );
		free( strain_G.zz );
		free( strain_G.yz );
		free( strain_G.xz );
		free( strain_G.xy );


		for( int idx_shape = 0; idx_shape < N_shape; idx_shape++ )
		{
		
			free( convolutions_char_strain[idx_shape].char_xx );
			free( convolutions_char_strain[idx_shape].char_yy );
			free( convolutions_char_strain[idx_shape].char_zz );
			free( convolutions_char_strain[idx_shape].char_yz );
			free( convolutions_char_strain[idx_shape].char_xz );
			free( convolutions_char_strain[idx_shape].char_xy );
			free( convolutions_char_strain[idx_shape].char_hydrostatic );
			free( convolutions_char_strain[idx_shape].char_biaxial );
		
		}

		if ( include_piezo != 0 )
		{

			free( phi_piezo_G.pz );

		}

		free( convolutions_char_strain );
	
	}



	////////////////////////////////////////////////
	//// Calculate Coulomb matrix elements /////////
	////////////////////////////////////////////////

	if ( include_CME >= 1 )
	{

		if ( myid == 0 )
		{

			printf( "Calculating Coulomb matrix elements...\n");

		}

		int N_CME_ee = 2*N_states[1]*2*N_states[1]*2*N_states[1]*2*N_states[1]; // Number of Coulomb matrix elements involving electron states only (including symmetries)
		int N_CME_hh = 2*N_states[0]*2*N_states[0]*2*N_states[0]*2*N_states[0]; // Number of Coulomb matrix elements involving hole states only (including symmetries)
		int N_CME_eh = 2*N_states[0]*2*N_states[1]*2*N_states[0]*2*N_states[1]; // Number of direct and exchange Coulomb matrix elements involving electron and hole states (including symmetries)

		// Unused variables that give corresponding numbers of unique Coulomb matrix elements when taking symmetries into account
		//int N_CME_ee_unique = ( ( ( 2*N_states[1] )*( 2*N_states[1] + 1 ) ) / 2 ) * ( ( ( 2*N_states[1] )*( 2*N_states[1] + 1 ) ) / 2 ) - 2*N_states[1]*2*N_states[1]*( 2*N_states[1] - 1 ) / 2;
		//int N_CME_hh_unique = ( ( ( 2*N_states[0] )*( 2*N_states[0] + 1 ) ) / 2 ) * ( ( ( 2*N_states[0] )*( 2*N_states[0] + 1 ) ) / 2 ) - 2*N_states[0]*2*N_states[0]*( 2*N_states[0] - 1 ) / 2;
		//int N_CME_eh_unique = ( 2*N_states[0]*2*N_states[1] )*( 2*N_states[0]*2*N_states[1] + 1 ) / 2;

		struct coulomb_matrix_element CME;

		if ( myid == 0 )
		{
		
			CME.ee = calloc( N_CME_ee , sizeof( complex double ) );
			CME.hh = calloc( N_CME_hh , sizeof( complex double ) );
			CME.eh_direct = calloc( N_CME_eh , sizeof( complex double ) );
			CME.eh_exchange = calloc( N_CME_eh , sizeof( complex double ) );

		}
		else
		{
		
			CME.ee = NULL;
			CME.hh = NULL;
			CME.eh_direct = NULL;
			CME.eh_exchange = NULL;

		}


		// Coulomb supercell basis set
		struct wave_vectors coulomb_basis;

		coulomb_basis.N_Gx = N_C_Gx;
		coulomb_basis.N_Gy = N_C_Gy;
		coulomb_basis.N_Gz = N_C_Gz;

		coulomb_basis.N_Gx_T = 2*N_C_Gx - 1;
		coulomb_basis.N_Gy_T = 2*N_C_Gy - 1;
		coulomb_basis.N_Gz_T = 2*N_C_Gz - 1;

		coulomb_basis.Gx = ( double* ) malloc( coulomb_basis.N_Gx_T*sizeof( double ) );
		coulomb_basis.Gy = ( double* ) malloc( coulomb_basis.N_Gy_T*sizeof( double ) );
		coulomb_basis.Gz = ( double* ) malloc( coulomb_basis.N_Gz_T*sizeof( double ) );

		// Set Coulomb supercell basis set
		basis_set( coulomb_basis, supercell, 0, 1 );
			
		start = MPI_Wtime();

		// Calculate Coulomb matrix elements
		coulomb_matrix_elements( basis, coulomb_basis, supercell, parameters, N_states, eigenstate_vb, eigenstate_cb, CME, coulomb_region_type );

		finish = MPI_Wtime();
    	elapsed = finish - start;
		component_run_times.cme_total = elapsed;

		if ( myid == 0 )
		{
		
			printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );
			
		}

		if ( myid == 0 )
		{

			// Write Coulomb matrix elements to file
			write_coulomb_matrix_elements( 2*N_states[0], 2*N_states[1], CME );

			free ( CME.ee );
			free ( CME.hh );
			free ( CME.eh_direct );
			free ( CME.eh_exchange );

		}

		free( coulomb_basis.Gx );
		free( coulomb_basis.Gy );
		free( coulomb_basis.Gz );

	}



	////////////////////////////////////////////////
	//// Write code component calculation times ////
	////////////////////////////////////////////////

	double end_total = MPI_Wtime();
	component_run_times.total = end_total - start_total;

	write_calculation_times();
        
        
	// Move output files to output files folder
	
	if ( myid == 0 )
	{
	
		char copy_supercell[80] = "mv supercell.dat "; // Moves supercell.dat to output files folder 
		strcat( copy_supercell, graph_folder );
		system( copy_supercell );

		char copy_basis_set[80] = "mv basis_set.dat "; // Moves basis_set.dat to output files folder 
		strcat( copy_basis_set, graph_folder );
		system( copy_basis_set );
		
		char copy_cb_energies[80] = "mv energies_cb.dat "; // Moves energies_cb.dat to output files folder 
		strcat( copy_cb_energies, graph_folder );
		system( copy_cb_energies );
		
		char copy_vb_energies[80] = "mv energies_vb.dat "; // Moves energies_vb.dat to output files folder 
		strcat( copy_vb_energies, graph_folder );
		system( copy_vb_energies );

		if ( include_overlaps )
		{

			char copy_overlaps[80] = "mv overlaps.dat "; // Moves overlaps.dat to output files folder 
			strcat( copy_overlaps, graph_folder );
			system( copy_overlaps );

		}
		
		if ( include_localisation )
		{

			char copy_carrier_localisation_cb[80] = "mv carrier_localisation_cb.dat "; // Moves carrier_localisation_cb.dat to output files folder 
			strcat( copy_carrier_localisation_cb, graph_folder );
			system( copy_carrier_localisation_cb );
			
			char copy_carrier_localisation_vb[80] = "mv carrier_localisation_vb.dat "; // Moves carrier_localisation_vb.dat to output files folder 
			strcat( copy_carrier_localisation_vb, graph_folder );
			system( copy_carrier_localisation_vb );

		}
		
		if ( include_MME )
		{
		
			char copy_momentum_matrix_elements[80] = "mv momentum_matrix_elements.dat "; // Moves momentum_matrix_elements.dat to output files folder 
			strcat( copy_momentum_matrix_elements, graph_folder );
			system( copy_momentum_matrix_elements );

		}

		char copy_calculation_times[80] = "mv calculation_times.dat "; // Moves calculation_times.dat to output files folder 
		strcat( copy_calculation_times, graph_folder );
		system( copy_calculation_times );

		if ( include_CME >= 1 )
		{

			char copy_cme_hh[80] = "mv coulomb_matrix_elements_hh.dat "; // Moves calculation_times.dat to output files folder 
			strcat( copy_cme_hh, graph_folder );
			system( copy_cme_hh );

			char copy_cme_ee[80] = "mv coulomb_matrix_elements_ee.dat "; // Moves calculation_times.dat to output files folder 
			strcat( copy_cme_ee, graph_folder );
			system( copy_cme_ee );

			char copy_cme_eh_direct[80] = "mv coulomb_matrix_elements_eh_direct.dat "; // Moves calculation_times.dat to output files folder 
			strcat( copy_cme_eh_direct, graph_folder );
			system( copy_cme_eh_direct );

			char copy_cme_eh_exchange[80] = "mv coulomb_matrix_elements_eh_exchange.dat "; // Moves calculation_times.dat to output files folder 
			strcat( copy_cme_eh_exchange, graph_folder );
			system( copy_cme_eh_exchange );

		}
	
	}
	


	///////////////////////////////////////////
	//// Free dynamically allocated memory ////
	///////////////////////////////////////////

	free( directions );

	free( basis.Gx   );
	free( basis.Gy   );
	free( basis.Gz   );
	free( supercell  );
	free( parameters );
	
	
	for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
	{

		free( eigenstate_vb[idx_state].fourier );
		free( eigenstate_vb[idx_state].bloch );

	}

	for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
	{

		free( eigenstate_cb[idx_state].fourier );
		free( eigenstate_cb[idx_state].bloch );

	}
	
	free( eigenstate_cb );
	free( eigenstate_vb );
	
	//free( overlaps );
	
	BLACS_GRIDEXIT(&my_blacs_ctxt);
	MPI_Finalize();

	if ( myid == 0 )
	{

		printf( "Run complete.\n");

	}

	return 0;

}
