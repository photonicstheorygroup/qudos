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


//////////////////////////
//// Global variables ////
//////////////////////////

// Descriptions given in 8_band.h

int N_shape, piezo_order, isotropic, include_shape_elastic_constants, shape_elastic_constants_order, softening;
double sigma_softening;
char graph_folder[80], *endptr;
struct convolution_char_strain *convolutions_char_strain;


int main( int argc, char *argv[] )
{

	struct timeval start, finish; // Variables to measure run time of the different calculation stages


	///////////////////////////////////////////////////////////////////////////////////////
	//// Read parameters.in file, command line parameters, and supercell.in input file ////
	///////////////////////////////////////////////////////////////////////////////////////

	printf( "Reading supercell parameters and setting up plane wave basis set...\n");

	// Parameters to be read from parameters.in / command line
	int idx_shape;					// variable to index shape for later loops over supercell shapes
	int N_Gx, N_Gy, N_Gz;			// No. of plane waves for envelope functions
	int strain_region_type;			// Free strain supercell lengths or commensurate strain supercell lengths
	int N_S_Gx, N_S_Gy, N_S_Gz;		// No. of plane waves for strain supercell (corresponding to above for envelope functions of electronic structure supercell)
	int N_grid;						// No. of real-space points for linescans
	int cross_sections;				// Enable calculation of real-space cross-sections
	int N_grid_cross_sections;		// No. of real-sapce points for cross-sections (same for both dimensions)
	double T;						// Temperature
	double Fx, Fy, Fz;				// Electric field component for each dimension

	int *directions = ( int* ) malloc( N_DIRECT*sizeof( int ) );	// Array of bools to indicate whether to compute linescans along each direction
	
	int *directions_cross_sections = ( int* ) malloc( N_DIRECT_CS*sizeof( int ) );	// Array of bools to indicate whether to compute cross-sections with normals along each direction


	int *input_check = ( int* ) calloc( N_PARAM, sizeof( int ) ); // Array to check that all input parameters have been provided

	// Read calculation parameters
	read_parameters( &N_shape, &N_Gx, &N_Gy, &N_Gz, &strain_region_type, &N_S_Gx, &N_S_Gy, &N_S_Gz, &T, &piezo_order, &Fx, &Fy, &Fz, &N_grid, &isotropic, &include_shape_elastic_constants, &shape_elastic_constants_order, &softening, &sigma_softening, directions, &cross_sections, &N_grid_cross_sections, directions_cross_sections, input_check );

	// Read command line inputs
	read_command_line_input( argc, argv, &N_shape, &N_Gx, &N_Gy, &N_Gz, &strain_region_type, &N_S_Gx, &N_S_Gy, &N_S_Gz, &T, &piezo_order, &Fx, &Fy, &Fz, &N_grid, &isotropic, &include_shape_elastic_constants, &shape_elastic_constants_order, &softening, &sigma_softening, directions, &cross_sections, &N_grid_cross_sections, directions_cross_sections, input_check );

	free( input_check );
	
	
	time_t t = time( NULL );
	struct tm current_time;
	current_time = *localtime(&t);
	sprintf( graph_folder, "Output_files_linescan_%d_%02d_%02d_%02d_%02d_%02d", current_time.tm_year+1900, current_time.tm_mon+1, current_time.tm_mday, current_time.tm_hour, current_time.tm_min, current_time.tm_sec ); // Graphs folder will be named using the current time


	struct supercell_geometry *supercell = ( _supercell_geometry* ) malloc( N_shape*sizeof( _supercell_geometry ) );	// To hold supercell and shape parameters

	// Read supercell parameters
	read_supercell( N_shape, strain_region_type, supercell );


	// In the case of use of relative strain supercell lengths, check if the number of plane waves used for the strain supercell is above the minimum allowed.
	if ( strain_region_type == 1 )
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

			exit( 1 );

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
	
	basis.N_Gx_T = N_Gx;
	basis.N_Gy_T = N_Gy;
	basis.N_Gz_T = N_Gz;

	basis.Gx = ( double* ) malloc( basis.N_Gx_T*sizeof( double ) );
	basis.Gy = ( double* ) malloc( basis.N_Gy_T*sizeof( double ) );
	basis.Gz = ( double* ) malloc( basis.N_Gz_T*sizeof( double ) );

	// Compute basis set for electronic structure supercell
	basis_set( basis, supercell, 0 );
	      
	struct wave_vectors strain_basis;

	strain_basis.N_Gx = N_S_Gx;
	strain_basis.N_Gy = N_S_Gy;
	strain_basis.N_Gz = N_S_Gz;

	strain_basis.N_Gx_T = N_S_Gx;
	strain_basis.N_Gy_T = N_S_Gy;
	strain_basis.N_Gz_T = N_S_Gz;

	strain_basis.Gx = ( double* ) malloc( strain_basis.N_Gx_T*sizeof( double ) );
	strain_basis.Gy = ( double* ) malloc( strain_basis.N_Gy_T*sizeof( double ) );
	strain_basis.Gz = ( double* ) malloc( strain_basis.N_Gz_T*sizeof( double ) );

	// Compute basis set for strain supercell
	basis_set( strain_basis, supercell, 1 );
	
	// Write electronic structure supercell and strain supercell basis set to file
	write_basis_set( basis, strain_basis, supercell );


	////////////////////////////////////////////////////
	//// Interpolate material parameters for alloys ////
	////////////////////////////////////////////////////

	struct material_params *parameters = ( _material_params* ) malloc( N_shape*sizeof( _material_params ) ); // Holds material parameters of supercell and shapes

	// Computes material parameters of each shape in supercell
	material_parameters( N_shape, T,             supercell, parameters );

	// Writes supercell and shape parameters to file
	write_supercell    ( N_shape, T, Fx, Fy, Fz, supercell, parameters );

	printf( "Done.\n\n");


	//////////////////////////////////////////////////////////////////////
	//// Calculate shape characteristic functions in reciprocal space ////
	//////////////////////////////////////////////////////////////////////

	gettimeofday( &start, 0 );

	printf( "Calculating shape characteristic function Fourier transforms...\n");

	// Stores shape characteristic functions for electronic structure supercell basis set
	struct characteristic_functions *characteristic_function_G = ( _characteristic_functions* ) malloc( N_shape*sizeof( _characteristic_functions ) );

	for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
	{

		characteristic_function_G[idx_shape].chi = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );

	}

	// Compute shape charcteristic functions in reciprocal space for electronic structure supercell
	characteristic_functions_reciprocal_space( N_shape, basis, supercell, characteristic_function_G, 0 ); // Compute shape charcteristic functions in reciprocal space
	
	// Stores shape characteristic functions for the strain supercell basis set
	struct characteristic_functions *characteristic_function_strain_G = ( _characteristic_functions* ) malloc( N_shape*sizeof( _characteristic_functions ) );

	int shape_index_start; // Indicates whether to calculate Fourier tranform of electronic structure supercell characteristic function inside the strain supercell, if using efficient FFT-based method (commensurate supercell lengths) to downsample strain
	if ( strain_region_type == 0 ) // Free strain supercell lengths
	{

		shape_index_start = 1;

	}
	else // Commensurate strain supercell lengths with electronic structure supercell lengths
	{

		shape_index_start = 0;

	}

	for( idx_shape = shape_index_start; idx_shape < N_shape; idx_shape++ )
	{

		characteristic_function_strain_G[idx_shape].chi = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );

	}

	// Compute shape charcteristic functions in reciprocal space for strain supercell
	characteristic_functions_reciprocal_space( N_shape, strain_basis, supercell, characteristic_function_strain_G, strain_region_type + 1 ); // Compute shape charcteristic functions in reciprocal space

	gettimeofday( &finish, 0 );
	long sec = finish.tv_sec  - start.tv_sec;
	long usec = finish.tv_usec - start.tv_usec;
	double elapsed = sec + usec*1.0e-06;

	printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );


	//////////////////////////////////////////////////////////////////////////
	//// Calculate strain tensor in reciprocal space for strain supercell ////
	//////////////////////////////////////////////////////////////////////////

	printf( "Calculating strain in reciprocal space for strain supercell...\n");
	
	struct strain_tensor strain_region_G; // Stores strain tensor in reciprocal space

	strain_region_G.xx = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
	strain_region_G.yy = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
	strain_region_G.zz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
	strain_region_G.yz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
	strain_region_G.xz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
	strain_region_G.xy = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );

	gettimeofday( &start, 0 );

	// Compute strain tensor components in reciprocal space for strain supercell
	strain_region_reciprocal_space( N_shape, strain_basis, supercell, parameters, strain_region_G, characteristic_function_strain_G ); // Compute strain tensor components in reciprocal space

	gettimeofday( &finish, 0 );
	sec = finish.tv_sec  - start.tv_sec;
	usec = finish.tv_usec - start.tv_usec;
	elapsed = sec + usec*1.0e-06;

	printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );


	////////////////////////////////////////////////////////////////////////////////////
	//// Calculate piezoelectric potential in reciprocal space for strain supercell ////
	////////////////////////////////////////////////////////////////////////////////////

	complex double *phi_piezo_region_G; // Stores piezoelectric potential in recipriocal space

	if( piezo_order != 0 ) // If piezoelectric potential calculation is enabled
	{

		printf( "Calculating piezoelectric potential in reciprocal space for strain supercell...\n");

		phi_piezo_region_G = ( complex double* ) calloc( strain_basis.N_Gx*strain_basis.N_Gy*strain_basis.N_Gz, sizeof( complex double ) );

		gettimeofday( &start, 0 );

		// Compute piezoelectric potential in reciprocal space for strain supercell
		piezo_region_reciprocal_space( piezo_order, N_shape, strain_basis, supercell, parameters, characteristic_function_strain_G, strain_region_G, phi_piezo_region_G );

		gettimeofday( &finish, 0 );
		sec = finish.tv_sec  - start.tv_sec;
		usec = finish.tv_usec - start.tv_usec;
		elapsed = sec + usec*1.0e-06;

		printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );

	}


	////////////////////////////////////////////////////////////////////////////////////////
	//// Calculate strain tensor in reciprocal space for electronic structure supercell ////
	////////////////////////////////////////////////////////////////////////////////////////

	// Downsampling strain tensor components for electronic structure supercell

	printf( "Downsampling strain for electronic structure supercell...\n");

	struct strain_tensor strain_G; // Strain tensor in reciprocal space

	strain_G.xx = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
	strain_G.yy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
	strain_G.zz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
	strain_G.yz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
	strain_G.xz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
	strain_G.xy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );

	gettimeofday( &start, 0 );

	if ( strain_region_type == 0 ) // If using free strain supercell lengths
	{

		strain_reciprocal_space( basis, strain_basis, supercell, strain_G, strain_region_G ); // Compute strain tensor components in reciprocal space

	}
	else // If using commensurate strain and electronic structure supercell lengths
	{

		strain_reciprocal_space_fft( basis, strain_basis, supercell, strain_G, strain_region_G, characteristic_function_strain_G ); // Compute strain tensor components in reciprocal space

	}

	gettimeofday( &finish, 0 );
	sec = finish.tv_sec  - start.tv_sec;
	usec = finish.tv_usec - start.tv_usec;
	elapsed = sec + usec*1.0e-06;

	printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );


	free( strain_region_G.xx );
	free( strain_region_G.yy );
	free( strain_region_G.zz );
	free( strain_region_G.yz );
	free( strain_region_G.xz );
	free( strain_region_G.xy );


	//////////////////////////////////////////////////////////////////////////////////////////////////
	//// Calculate piezoelectric potential in reciprocal space for electronic structure supercell ////
	//////////////////////////////////////////////////////////////////////////////////////////////////

	// Downsampling piezoelectric potential for electronic structure supercell

	complex double *phi_piezo_G; // Piezoelectric potential in reciprocal space


	if( piezo_order != 0 ) // If piezoelectric potential is included in the calculation
	{

		printf( "Downsampling piezoelectric potential for electronic structure supercell...\n");

		phi_piezo_G = ( complex double* ) calloc( basis.N_Gx*basis.N_Gy*basis.N_Gz, sizeof( complex double ) );

		gettimeofday( &start, 0 );

		if ( strain_region_type == 0 ) // If using free strain supercell lengths
		{

			piezo_reciprocal_space( basis, strain_basis, supercell, phi_piezo_G, phi_piezo_region_G ); // Compute strain tensor components in reciprocal space

		}
		else // If using commensurate strain and electronic structure supercell lengths
		{

			piezo_reciprocal_space_fft( basis, strain_basis, supercell, phi_piezo_G, phi_piezo_region_G, characteristic_function_strain_G ); // Compute strain tensor components in reciprocal space

		}

		gettimeofday( &finish, 0 );
		sec = finish.tv_sec  - start.tv_sec;
		usec = finish.tv_usec - start.tv_usec;
		elapsed = sec + usec*1.0e-06;

		printf( "Done. Time to compute = %9e seconds.\n\n", elapsed );

		// Free piezoelectric potential for strain supercell as no longer required
		free( phi_piezo_region_G );

	}

	// Free characteristic functions for strain supercell as no longer required
	for( idx_shape = shape_index_start; idx_shape < N_shape; idx_shape++ )
	{

		free( characteristic_function_strain_G[idx_shape].chi );

	}
	free( characteristic_function_strain_G );


	////////////////////////////////////////
	//// Calculate real space linescans ////
	////////////////////////////////////////

	printf( "Calculating real space linescans...\n");

	// Direction indicator
	int direction;

	// Real space grid
	double *x = ( double* ) malloc( N_grid*sizeof( double ) );
	double *y = ( double* ) malloc( N_grid*sizeof( double ) );
	double *z = ( double* ) malloc( N_grid*sizeof( double ) );

	// Strain tensor in real space
	struct strain_tensor strain_r;

	strain_r.xx = ( complex double* ) calloc( N_grid, sizeof( complex double ) );
	strain_r.yy = ( complex double* ) calloc( N_grid, sizeof( complex double ) );
	strain_r.zz = ( complex double* ) calloc( N_grid, sizeof( complex double ) );
	strain_r.yz = ( complex double* ) calloc( N_grid, sizeof( complex double ) );
	strain_r.xz = ( complex double* ) calloc( N_grid, sizeof( complex double ) );
	strain_r.xy = ( complex double* ) calloc( N_grid, sizeof( complex double ) );

	// Piezoelectric potential in real space
	complex double *piezo_r = ( complex double* ) malloc( N_grid*sizeof( complex double ) );

	// Band-edge energies in real space
	struct band_edge_energies band_edges;

	band_edges.E_so = ( double* ) calloc( N_grid, sizeof( double ) );
	band_edges.E_lh = ( double* ) calloc( N_grid, sizeof( double ) );
	band_edges.E_hh = ( double* ) calloc( N_grid, sizeof( double ) );
	band_edges.E_cb = ( double* ) calloc( N_grid, sizeof( double ) );


	// Loop over directions
	for( direction = 0; direction < N_DIRECT; direction++ )
	{

		if( directions[direction] == 1 ) // Compute real space linescans along current direction
		{

			gettimeofday( &start, 0 );

			// Set real space grid for current linescan direction
			real_space_grid_linescan( direction, N_grid, supercell, x, y, z );


			// Strain linescan
			strain_real_space    (   N_shape, N_grid, x, y, z, basis, supercell, parameters, strain_G, strain_r );
			write_strain_linescan( direction, N_grid, x, y, z,                                         strain_r );


			// Piezoelectric potential linescan
			if( piezo_order != 0 )
			{

				piezo_real_space    (   N_shape, N_grid, x, y, z, basis, supercell, parameters, phi_piezo_G, piezo_r );
				write_piezo_linescan( direction, N_grid, x, y, z,                               piezo_order, piezo_r );

			}


			// Band edge energies linescan
			band_edge_energies_real_space    (    N_grid, N_shape, x, y, z, Fx, Fy, Fz, basis, supercell, parameters, characteristic_function_G, strain_r, piezo_r, band_edges );
			write_band_edge_energies_linescan( direction,  N_grid, x, y, z,                                                                                         band_edges );
			

			// Copy linescan files to output folder
			copy_output( 0, direction );

			gettimeofday( &finish, 0 );
		    	 sec = finish.tv_sec  - start.tv_sec;
		    	usec = finish.tv_usec - start.tv_usec;
		    	elapsed = sec + usec*1.0e-06;

		    printf( "\tTime to compute for direction %d = %9e seconds.\n", direction, elapsed );

		}


	}
	
	
	/////////////////////////////////////////////
	//// Calculate real space cross-sections ////
	/////////////////////////////////////////////
	
	
	if ( cross_sections )
	{

		printf( "\nCalculating real space cross-sections...\n");
	
		int direction; // Direction of normal 
		
		int N_grid_cross_sections_squared = N_grid_cross_sections * N_grid_cross_sections;
		
		// Real space grid
		double *x = ( double* ) malloc( N_grid_cross_sections_squared*sizeof( double ) );
		double *y = ( double* ) malloc( N_grid_cross_sections_squared*sizeof( double ) );
		double *z = ( double* ) malloc( N_grid_cross_sections_squared*sizeof( double ) );

		// Strain tensor in real space
		struct strain_tensor strain_r;

		strain_r.xx = ( complex double* ) calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
		strain_r.yy = ( complex double* ) calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
		strain_r.zz = ( complex double* ) calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
		strain_r.yz = ( complex double* ) calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
		strain_r.xz = ( complex double* ) calloc( N_grid_cross_sections_squared, sizeof( complex double ) );
		strain_r.xy = ( complex double* ) calloc( N_grid_cross_sections_squared, sizeof( complex double ) );

		// Piezoelectric potential in real space
		complex double *piezo_r = ( complex double* ) malloc( N_grid_cross_sections_squared*sizeof( complex double ) );

		// Band-edge energies in real space
		struct band_edge_energies band_edges;

		band_edges.E_so = ( double* ) calloc( N_grid_cross_sections_squared, sizeof( double ) );
		band_edges.E_lh = ( double* ) calloc( N_grid_cross_sections_squared, sizeof( double ) );
		band_edges.E_hh = ( double* ) calloc( N_grid_cross_sections_squared, sizeof( double ) );
		band_edges.E_cb = ( double* ) calloc( N_grid_cross_sections_squared, sizeof( double ) );


		// Loop over normal directions
		for( direction = 0; direction < N_DIRECT_CS; direction++ )
		{

			if( directions_cross_sections[direction] == 1 ) // Compute real space cross sections with normal along current direction
			{

				gettimeofday( &start, 0 );

				// Set real space grid for cross-section with current normal direction
				real_space_grid_cross_section( direction, N_grid_cross_sections, supercell, x, y, z );


				// Strain cross section
				strain_real_space(   N_shape, N_grid_cross_sections_squared, x, y, z, basis, supercell, parameters, strain_G, strain_r );
				write_strain_cross_section( direction, N_grid_cross_sections, x, y, z,                                         strain_r );


				// Piezoelectric potential linescan
				if( piezo_order != 0 )
				{

					piezo_real_space    (   N_shape, N_grid_cross_sections_squared, x, y, z, basis, supercell, parameters, phi_piezo_G, piezo_r );
					write_piezo_cross_section( direction, N_grid_cross_sections, x, y, z,                               piezo_order, piezo_r );

				}


				// Band edge energies linescan
				band_edge_energies_real_space    (    N_grid_cross_sections_squared, N_shape, x, y, z, Fx, Fy, Fz, basis, supercell, parameters, characteristic_function_G, strain_r, piezo_r, band_edges );
				write_band_edge_energies_cross_section( direction,  N_grid_cross_sections, x, y, z,                                                                                         band_edges );
			

				// Copy cross-section files to output folder
				copy_output( 1, direction );

				gettimeofday( &finish, 0 );
			    	 sec = finish.tv_sec  - start.tv_sec;
			    	usec = finish.tv_usec - start.tv_usec;
			    	elapsed = sec + usec*1.0e-06;

			    	printf( "\tTime to compute for normal along direction %d = %9e seconds.\n", direction, elapsed );

			}


		}
		
		free( x );
		free( y );
		free( z );

		free( strain_r.xx );
		free( strain_r.yy );
		free( strain_r.zz );
		free( strain_r.yz );
		free( strain_r.xz );
		free( strain_r.xy );

		free( piezo_r );

		free( band_edges.E_so );
		free( band_edges.E_lh );
		free( band_edges.E_hh );
		free( band_edges.E_cb );
	
	}

	
    char copy_supercell[80] = "mv supercell.dat "; // Copies supercell.dat to graphs folder 
	strcat( copy_supercell, graph_folder );
	system( copy_supercell );

	char copy_basis_set[80] = "mv basis_set.dat "; // Copies basis_set.dat to graphs folder 
	strcat( copy_basis_set, graph_folder );
	system( copy_basis_set );

	///////////////////////////////////////////
	//// Free dynamically allocated memory ////
	///////////////////////////////////////////

	free( directions );
	free( directions_cross_sections );

	free( basis.Gx   );
	free( basis.Gy   );
	free( basis.Gz   );
	free( strain_basis.Gx   );
	free( strain_basis.Gy   );
	free( strain_basis.Gz   );
	free( supercell  );
	free( parameters );

	for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
	{

		free( characteristic_function_G[idx_shape].chi );

	}

	free( characteristic_function_G );
	

	free( strain_G.xx );
	free( strain_G.yy );
	free( strain_G.zz );
	free( strain_G.yz );
	free( strain_G.xz );
	free( strain_G.xy );



	free( x );
	free( y );
	free( z );

	free( strain_r.xx );
	free( strain_r.yy );
	free( strain_r.zz );
	free( strain_r.yz );
	free( strain_r.xz );
	free( strain_r.xy );

	if( piezo_order != 0 )
	{

		free( phi_piezo_G );

	}

	free( piezo_r );

	free( band_edges.E_so );
	free( band_edges.E_lh );
	free( band_edges.E_hh );
	free( band_edges.E_cb );

	printf( "\nRun complete.\n");

	return 0;

}
