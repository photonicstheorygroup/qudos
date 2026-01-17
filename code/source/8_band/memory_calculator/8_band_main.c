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

char *endptr, graph_folder[256];

int N_shape, isotropic, softening, include_strain, include_piezo, include_shape_elastic_constants, shape_elastic_constants_order, optimise, include_overlaps, include_localisation, include_CME, include_MME, PACKED, Num_Possible_Polarisations, Polarisation_Indicator;

double sigma_softening;

MKL_INT N_H;

long long local_mem;
long long local_mem_peak;

struct characteristic_functions *characteristic_function_G;
struct characteristic_functions *characteristic_function_strain_G;
struct strain_tensor strain_G;
struct strain_tensor strain_region_G;
struct convolution_char_strain *convolutions_char_strain;
struct piezo_potential phi_piezo_G;
struct piezo_potential phi_piezo_region_G;
struct timing component_run_times;
struct memory component_memory_usage;

complex double *dH_dk, *y;



///////////////////////////////
//// MPI & BLACS variables ////
///////////////////////////////

// Descriptions given in 8_band.h

int myid, nprocs;
MKL_INT na_cols, na_rows, na_colsz, na_rowsz;
MKL_INT np_cols, np_rows;
MKL_INT my_prow, my_pcol;
MKL_INT mpi_comm;

MKL_INT my_blacs_ctxt, sc_desc[9], sc_descz[9], info;

MKL_INT block_size;


int main( int argc, char *argv[] )
{
	
	double start, finish, elapsed; // Variables to measure runtime of the different calculation steps

	local_mem = 0;
	local_mem_peak = 0;

	component_memory_usage.char_mem = 0;
	component_memory_usage.strain_mem = 0;
	component_memory_usage.piezo_mem = 0;
	component_memory_usage.strain_downsample_mem = 0;
	component_memory_usage.convolutions_mem = 0;
	component_memory_usage.elec_struc_mem = 0;
	component_memory_usage.eig_work_mem = 0;
	component_memory_usage.mme_mem = 0;
	component_memory_usage.real_space_mem = 0;
	component_memory_usage.CME_mem = 0;

	
	////////////////////////////////////
	//// MPI & BLACS initialisation ////
	////////////////////////////////////
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	// Optimally setting numbers of process rows and columns
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


	if ( myid == 0 )
	{

		printf( "Calculating peak memory usage for each code stage...\n");

	}


	///////////////////////////////////////////////////////////////////////////////////////
	//// Read parameters.in file, command line parameters, and supercell.in input file ////
	///////////////////////////////////////////////////////////////////////////////////////

	// Parameters to be read from command line
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

	
	/////////////////////////////////////
	//// Set up plane wave basis set ////
	/////////////////////////////////////
	
	// Electronic structure calculation basis set
	struct wave_vectors basis;

	basis.N_Gx = N_Gx;
	basis.N_Gy = N_Gy;
	basis.N_Gz = N_Gz;
	
	if ( include_strain == 1)
	{

		basis.N_Gx_T = 2*N_Gx - 1;
		basis.N_Gy_T = 2*N_Gy - 1;
		basis.N_Gz_T = 2*N_Gz - 1;
	
	}
	
	else
	{
	
		basis.N_Gx_T = 2*N_Gx - 1;
		basis.N_Gy_T = 2*N_Gy - 1;
		basis.N_Gz_T = 2*N_Gz - 1;
	
	}
	
	basis.Gx = ( double* ) malloc( basis.N_Gx_T*sizeof( double ) );
	basis.Gy = ( double* ) malloc( basis.N_Gy_T*sizeof( double ) );
	basis.Gz = ( double* ) malloc( basis.N_Gz_T*sizeof( double ) );

	local_mem += ( basis.N_Gx_T + basis.N_Gy_T + basis.N_Gz_T ) *sizeof( double ); // Adding memory to store basis set

	basis.Kx = Kx;
	basis.Ky = Ky;
	basis.Kz = Kz;

	// Compute basis set for electronic structure supercell
	basis_set( basis, supercell, 0, 0 );
	
	// Strain region basis set
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

		local_mem += ( strain_basis.N_Gx_T + strain_basis.N_Gy_T + strain_basis.N_Gz_T ) *sizeof( double );

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
	
	DESCINIT(sc_desc, &N_H, &N_H, &block_size, &block_size, &f, &f, &my_blacs_ctxt, &na_rows, &info); // Initialising array desriptor for distributed Hamiltonian matrix
	DESCINIT(sc_descz, &N_H, &N_H, &block_size, &block_size, &f, &f, &my_blacs_ctxt, &na_rows, &info); // Initialising array desriptor for distributed eigenvector matrix
	

	

	struct material_params *parameters = ( _material_params* ) malloc( N_shape*sizeof( _material_params ) ); // Holds material parameters of supercell and shapes

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

	
	//////////////////////////////////////////////////////////////////////
	//// Calculate shape characteristic functions in reciprocal space ////
	//////////////////////////////////////////////////////////////////////

	if ( include_strain == 1 )
	{
	
		// Characteristic functions for the strain region basis set
		characteristic_function_strain_G = ( _characteristic_functions* ) malloc( N_shape*sizeof( _characteristic_functions ) ); // Shape characteristic functions in reciprocal space for larger strain region

		int shape_index_start;
		if ( strain_region_type == 0 )
		{

			shape_index_start = 1;

		}
		else
		{

			shape_index_start = 0;

		}

		for( int idx_shape = shape_index_start; idx_shape < N_shape; idx_shape++ )
		{

			//characteristic_function_strain_G[idx_shape].chi = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );

		}

		local_mem += ( N_shape - shape_index_start )*( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T ) *sizeof( complex double ); // Adding memory to store shape characteristic functions for strain supercell

		characteristic_functions_reciprocal_space( N_shape, strain_basis, supercell, strain_region_type + 1 ); // Compute shape charcteristic functions in reciprocal space

		component_memory_usage.char_mem = local_mem_peak;
	
	}
	
	// Characteristic functions for the electronic structure calculation basis set
	characteristic_function_G = ( _characteristic_functions* ) malloc( N_shape*sizeof( _characteristic_functions ) ); // Shape characteristic functions in reciprocal space

	for( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
	{
		
		//characteristic_function_G[idx_shape].chi = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		
	}
	local_mem += ( N_shape - 1 )*( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T ) *sizeof( complex double ); // Adding memory to store shape characteristic functions for electronic structure supercell

	if ( include_strain == 1 && basis.N_Gx == strain_basis.N_Gx && supercell[0].L1 == supercell[0].SL1 && basis.N_Gy == strain_basis.N_Gy && supercell[0].L2 == supercell[0].SL2 && basis.N_Gz == strain_basis.N_Gz && supercell[0].L3 == supercell[0].SL3 )
	{
	
		for( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
		{

			for ( int idx_basis = 0; idx_basis < basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T; idx_basis++ )
			{
			
				//characteristic_function_G[idx_shape].chi[idx_basis] = characteristic_function_strain_G[idx_shape].chi[idx_basis];
			
			}

		}
	
	}
	else
	{

		characteristic_functions_reciprocal_space( N_shape, basis, supercell, 0 ); // Compute shape charcteristic functions in reciprocal space
		if ( include_strain == 0 )
		{

			component_memory_usage.char_mem = local_mem_peak;

		}

	}
	
	
	//////////////////////////////////////////////////////////////////////////
	//// Calculate strain tensor in reciprocal space for strain supercell ////
	//////////////////////////////////////////////////////////////////////////
	
	if ( include_strain == 1)
	{

		// Calculating strain tensor components for strain region

		//strain_region_G.xx = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
		//strain_region_G.yy = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
		//strain_region_G.zz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
		//strain_region_G.yz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
		//strain_region_G.xz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );
		//strain_region_G.xy = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );

		local_mem += ( 6 )*( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T ) *sizeof( complex double );  // Adding memory to store strain tensor components for strain supercell

		strain_region_reciprocal_space( N_shape, strain_basis, supercell, parameters ); // Compute strain tensor components in reciprocal space

		component_memory_usage.strain_mem = local_mem_peak;

	}


	////////////////////////////////////////////////////////////////////////////////////
	//// Calculate piezoelectric potential in reciprocal space for strain supercell ////
	////////////////////////////////////////////////////////////////////////////////////

	if ( include_strain == 1)
	{

		if( include_piezo != 0 )
		{

			//phi_piezo_region_G.pz = ( complex double* ) calloc( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T, sizeof( complex double ) );

			local_mem += ( 1 )*( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T ) *sizeof( complex double ); // Adding memory to store piezoelectric potential for strain supercell

			piezo_region_reciprocal_space( N_shape, strain_basis, supercell, parameters );	// Compute piezoelectric potential for strain supercell

			component_memory_usage.piezo_mem = local_mem_peak;

		}

		local_mem -= ( N_shape - 1 )*( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T ) *sizeof( complex double );

	}
		

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Calculate strain tensor and piezoelectric potential in reciprocal space for electronic structure supercell ////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	local_mem_peak = local_mem;

	if ( include_strain == 1)
	{
		
		// Calculating strain tensor components for electronic structure calculation
		
		//strain_G.xx = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		//strain_G.yy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		//strain_G.zz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		//strain_G.yz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		//strain_G.xz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );
		//strain_G.xy = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );

		local_mem += ( 6 )*( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T ) *sizeof( complex double ); // Adding memory to store global strain tensor components for electronic structure supercell

		if ( strain_region_type == 0 )
		{

			//strain_reciprocal_space( basis, strain_basis, supercell ); // Compute strain tensor components in reciprocal space

		}
		else
		{

			strain_reciprocal_space_fft( basis, strain_basis, supercell ); // Compute strain tensor components in reciprocal space

		}

		component_memory_usage.strain_downsample_mem = local_mem_peak;

		//free( strain_region_G.xx );
		//free( strain_region_G.yy );
		//free( strain_region_G.zz );
		//free( strain_region_G.yz );
		//free( strain_region_G.xz );
		//free( strain_region_G.xy );

		local_mem -= ( 6 )*( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T ) *sizeof( complex double );


		if( include_piezo != 0 )
		{

			//phi_piezo_G.pz = ( complex double* ) calloc( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T, sizeof( complex double ) );

			local_mem += ( 1 )*( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T ) *sizeof( complex double );	// Adding memory to store global piezoelectric potential for electronic structure supercell


			if ( strain_region_type == 0 )
			{

				//piezo_reciprocal_space( basis, strain_basis, supercell ); 

			}
			else
			{

				//piezo_reciprocal_space_fft( basis, strain_basis, supercell ); // Don't need to check, will be less than that used to downsample strain tensor components using same method

			}

			//free( phi_piezo_region_G.pz );

			local_mem -= ( 1 )*( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T ) *sizeof( complex double );

		}

		if ( strain_region_type == 1 )
		{

			//free( characteristic_function_strain_G[0].chi );
			local_mem -= ( 1 )*( strain_basis.N_Gx_T*strain_basis.N_Gy_T*strain_basis.N_Gz_T ) *sizeof( complex double );

		}
		//free( characteristic_function_strain_G );
		
	}
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	//// Calculate convolution of shape characteristic functions with strain tensor compoenents ////
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	if ( include_strain == 1)
	{
		
		convolutions_char_strain = ( _convolution_char_strain* ) malloc( N_shape*sizeof( _convolution_char_strain ) ); // Convolution between characteristic functions and strain tensor components

		for ( int idx_shape = 0; idx_shape < N_shape; idx_shape++ )
		{
		
			//convolutions_char_strain[idx_shape].char_xx = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			//convolutions_char_strain[idx_shape].char_yy = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			//convolutions_char_strain[idx_shape].char_zz = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			//convolutions_char_strain[idx_shape].char_yz = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			//convolutions_char_strain[idx_shape].char_xz = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			//convolutions_char_strain[idx_shape].char_xy = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			//convolutions_char_strain[idx_shape].char_hydrostatic = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			//convolutions_char_strain[idx_shape].char_biaxial = ( complex double* ) calloc( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1), sizeof( complex double ) );
			
		}

		local_mem += ( N_shape )*( 8 )*( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1) ) *sizeof( complex double );	// Adding memory to store global convolutions between strain tensor components and shape characteristic functions

		convolutions( N_shape, basis, supercell ); // Compute convolutions in reciprocal space
		

		component_memory_usage.convolutions_mem = local_mem_peak;
		
	}
	
	
	//////////////////////////////
	//// Electronic structure ////
	//////////////////////////////


	// Allocate memory to store eigenstates
	
	int N_states[2]; N_states[0] = N_vb; N_states[1] = N_cb;
	
	struct eigenstate *eigenstate_vb = malloc( 2*N_states[0]*sizeof( _eigenstate ) );
	struct eigenstate *eigenstate_cb = malloc( 2*N_states[1]*sizeof( _eigenstate ) );
	
	for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
	{

			//eigenstate_vb[idx_state].fourier       = calloc( na_rows     , sizeof( complex double ) );
			//eigenstate_vb[idx_state].bloch         = calloc( N_BAND  , sizeof(         double ) );

	}

	

	for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
	{

		//eigenstate_cb[idx_state].fourier       = calloc( na_rows     , sizeof( complex double ) );
		//eigenstate_cb[idx_state].bloch         = calloc( N_BAND  , sizeof(         double ) );
		
	}

	local_mem += ( 2*N_states[0] + 2*N_states[1] )*( na_rows ) *sizeof( complex double ); // Adding memory to store local portions of eigenvector matrix
	
	// Calculate supercell eigenstates

	eigenstates( N_shape, N_states, F, basis, supercell, parameters, eigenstate_vb, eigenstate_cb );


	component_memory_usage.elec_struc_mem = local_mem_peak;
	
	
	
	////////////////////////////////////////////
	//// Calculate momentum matrix elements ////
	////////////////////////////////////////////
	
	int Num_CB_States = N_states[1]; // Not including degeneracies
	int Num_VB_States = N_states[0]; // Not including degeneracies
	
	Num_Possible_Polarisations = 3;
	int polarisations[Num_Possible_Polarisations];
	polarisations[0] = 1;
	polarisations[1] = 1;
	polarisations[2] = 1;
	int Num_Polarisations = 0;
	for ( int i = 0; i < Num_Possible_Polarisations; i++ )
	{
	
		if ( polarisations[i] == 1 )
		{

			Num_Polarisations++;
			
		}
		
	}
	
	//complex double *momentum_matrix_elements = ( complex double* ) malloc( Num_Polarisations*2*Num_CB_States*2*Num_VB_States*sizeof(complex double) );

	local_mem += ( Num_Polarisations*2*Num_CB_States*2*Num_VB_States )*sizeof( complex double ); // Adding memory to store momentum matrix elements
	
	if (PACKED && nprocs == 1)
	{
	
		//dH_dk = ( complex double* ) calloc( N_H*(N_H + 1)/2, sizeof( complex double ) );
	
	}
	
	else
	{
	
		//dH_dk = ( complex double* ) malloc( na_rows*na_cols*sizeof( complex double ) );
		
		local_mem += ( na_rows*na_cols ) *sizeof( complex double );	// Adding memory to store local part of derivative of Hamiltonian matrix
		
	}
	
	//y = (complex double *) calloc( na_rowsz, sizeof( complex double ) );

	local_mem += ( na_rowsz ) *sizeof( complex double );	// Adding memory to store local part of y vector
	
	Polarisation_Indicator = -1;

	component_memory_usage.mme_mem = local_mem;
	
	
	//free( dH_dk );
	//free( y );

	local_mem -= ( na_rowsz ) *sizeof( complex double );
	local_mem -= ( na_rows*na_cols ) *sizeof( complex double );
	
	
	
	//////////////////////////////////////////////////////////////////////
	//// Calculate and write real space probability density linescans ////
	//////////////////////////////////////////////////////////////////////

	local_mem_peak = local_mem;

	int direction;

	double *x;// = ( double* ) malloc( N_grid*sizeof( double ) );
	double *y;// = ( double* ) malloc( N_grid*sizeof( double ) );
	double *z;// = ( double* ) malloc( N_grid*sizeof( double ) );
	local_mem += ( 3 )*( N_grid ) *sizeof( double );	// Adding memory to store real-space grid
	
	if ( myid == 0 )
	{
	
		/*for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
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

		}*/

		local_mem += ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid ) *sizeof( complex double ); // Adding memory to store real-space envelope functions
		local_mem += ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid ) *sizeof( double );			// Adding memory to store real-space probability densities
		local_mem_peak = local_mem;
	
	}
	
	if ( myid == 0 )
	{
	
		/*for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
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

		}*/

		local_mem -= ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid ) *sizeof( complex double );
		local_mem -= ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid ) *sizeof( double );
        
    }

	//free( x );
	//free( y );
	//free( z );

	local_mem -= ( 3 )*( N_grid ) *sizeof( double );
        
        
    ///////////////////////////////////////////////////////////////////////////
	//// Calculate and write real space probability density cross-sections ////
	///////////////////////////////////////////////////////////////////////////
	
	if (cross_sections)
	{
	
		int direction; // Direction of normal
		
		int N_grid_cross_sections_squared = N_grid_cross_sections * N_grid_cross_sections;

		double *x;// = ( double* ) malloc( N_grid_cross_sections_squared*sizeof( double ) );
		double *y;// = ( double* ) malloc( N_grid_cross_sections_squared*sizeof( double ) );
		double *z;// = ( double* ) malloc( N_grid_cross_sections_squared*sizeof( double ) );

		local_mem += ( 3 )*( N_grid_cross_sections_squared ) *sizeof( double );	// Adding memory to store real space grid
		
		if ( myid == 0 )
		{
		
			/*for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
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

			}*/

			local_mem += ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid_cross_sections_squared ) *sizeof( complex double );	// Adding memory to store real-space envelope functions
			local_mem += ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid_cross_sections_squared ) *sizeof( double );			// Adding memory to store real-space probability densities
			if ( local_mem > local_mem_peak )
			{
				
				local_mem_peak = local_mem;

			}
		
		}
		
		
		//free( x );
		//free( y );
		//free( z );

		local_mem -= ( 3 )*( N_grid_cross_sections_squared ) *sizeof( double );
		
		if( myid == 0 )
		{
		
			/*for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
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

			}*/

			local_mem -= ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid_cross_sections_squared ) *sizeof( complex double );
			local_mem -= ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid_cross_sections_squared ) *sizeof( double );
		
		}
        
    }


	//////////////////////////////////////////////////////////////////
	//// Calculate and write 3D real space probability densities  ////
	//////////////////////////////////////////////////////////////////
	
	if (enable_3D)
	{
		
		int N_grid_3D_cubed = N_grid_3D * N_grid_3D * N_grid_3D;

		double *x;// = ( double* ) malloc( N_grid_3D_cubed*sizeof( double ) );
		double *y;// = ( double* ) malloc( N_grid_3D_cubed*sizeof( double ) );
		double *z;// = ( double* ) malloc( N_grid_3D_cubed*sizeof( double ) );

		local_mem += ( 3 )*( N_grid_3D_cubed ) *sizeof( double ); // Adding memory to store real-space grid
		
		if ( myid == 0 )
		{
		
			/*for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
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

			}*/

			local_mem += ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid_3D_cubed ) *sizeof( complex double );	// Adding memory to store real-space envelope functions
			local_mem += ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid_3D_cubed ) *sizeof( double );			// Adding memory to store real-space probability densities
			if ( local_mem > local_mem_peak )
			{
				
				local_mem_peak = local_mem;

			}
		
		}
		
		//free( x );
		//free( y );
		//free( z );

		local_mem -= ( 3 )*( N_grid_3D_cubed ) *sizeof( double );
		
		if( myid == 0 )
		{
		
			/*for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
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

			}*/

			local_mem -= ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid_3D_cubed ) *sizeof( complex double );
			local_mem -= ( 16 )*( 2*N_states[0] + 2*N_states[1])*( N_grid_3D_cubed ) *sizeof( double );
		
		}
        
    }

	component_memory_usage.real_space_mem = local_mem_peak;


	///////////////////////////////////////////////////////////////////////////////////////////////
	//// Free dynamically allocated memory that is no longer required prior to CME calculation ////
	///////////////////////////////////////////////////////////////////////////////////////////////
	

	for( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
	{

		//free( characteristic_function_G[idx_shape].chi );

	}
	local_mem -= ( N_shape - 1 )*( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T ) *sizeof( complex double );

	free( characteristic_function_G );
	
	if (include_strain == 1)
	{
	
		free( strain_basis.Gx   );
		free( strain_basis.Gy   );
		free( strain_basis.Gz   );

		local_mem -= ( strain_basis.N_Gx_T + strain_basis.N_Gy_T + strain_basis.N_Gz_T ) *sizeof( double );
	
		//free( strain_G.xx );
		//free( strain_G.yy );
		//free( strain_G.zz );
		//free( strain_G.yz );
		//free( strain_G.xz );
		//free( strain_G.xy );

		local_mem -= ( 6 )*( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T ) *sizeof( complex double );


		for( int idx_shape = 0; idx_shape < N_shape; idx_shape++ )
		{
		
			//free( convolutions_char_strain[idx_shape].char_xx );
			//free( convolutions_char_strain[idx_shape].char_yy );
			//free( convolutions_char_strain[idx_shape].char_zz );
			//free( convolutions_char_strain[idx_shape].char_yz );
			//free( convolutions_char_strain[idx_shape].char_xz );
			//free( convolutions_char_strain[idx_shape].char_xy );
			//free( convolutions_char_strain[idx_shape].char_hydrostatic );
			//free( convolutions_char_strain[idx_shape].char_biaxial );
		
		}

		local_mem -= ( N_shape )*( 8 )*( (2*basis.N_Gx - 1)*(2*basis.N_Gy - 1)*(2*basis.N_Gz - 1) ) *sizeof( complex double );

		if ( include_piezo != 0 )
		{

			//free( phi_piezo_G.pz );

			local_mem -= ( 1 )*( basis.N_Gx_T*basis.N_Gy_T*basis.N_Gz_T ) *sizeof( complex double );

		}

		free( convolutions_char_strain );
	
	}


	////////////////////////////////////////////////
	//// Calculate Coulomb matrix elements /////////
	////////////////////////////////////////////////

	if ( include_CME >= 1 )
	{
	
		int N_CME_ee_unique = ( ( ( 2*N_states[1] )*( 2*N_states[1] + 1 ) ) / 2 ) * ( ( ( 2*N_states[1] )*( 2*N_states[1] + 1 ) ) / 2 ) - 2*N_states[1]*2*N_states[1]*( 2*N_states[1] - 1 ) / 2;
		int N_CME_hh_unique = ( ( ( 2*N_states[0] )*( 2*N_states[0] + 1 ) ) / 2 ) * ( ( ( 2*N_states[0] )*( 2*N_states[0] + 1 ) ) / 2 ) - 2*N_states[0]*2*N_states[0]*( 2*N_states[0] - 1 ) / 2;
		int N_CME_eh_unique = ( 2*N_states[0]*2*N_states[1] )*( 2*N_states[0]*2*N_states[1] + 1 ) / 2;

		int N_CME_ee = 2*N_states[1]*2*N_states[1]*2*N_states[1]*2*N_states[1];
		int N_CME_hh = 2*N_states[0]*2*N_states[0]*2*N_states[0]*2*N_states[0];
		int N_CME_eh = 2*N_states[0]*2*N_states[1]*2*N_states[0]*2*N_states[1];

		struct coulomb_matrix_element CME;

		if ( myid == 0 )
		{
		
			//CME.ee = calloc( N_CME_ee , sizeof( complex double ) );
			//CME.hh = calloc( N_CME_hh , sizeof( complex double ) );
			//CME.eh_direct = calloc( N_CME_eh , sizeof( complex double ) );
			//CME.eh_exchange = calloc( N_CME_eh , sizeof( complex double ) );

			local_mem += ( N_CME_ee + N_CME_hh + 2*N_CME_eh )*sizeof( complex double ); // Adding memory to store Coulomb matrix elements

		}
		
		local_mem_peak = local_mem;

		// Coulomb supercell basis set
		struct wave_vectors coulomb_basis;

		coulomb_basis.N_Gx = N_C_Gx;
		coulomb_basis.N_Gy = N_C_Gy;
		coulomb_basis.N_Gz = N_C_Gz;

		coulomb_basis.N_Gx_T = 2*N_C_Gx - 1;
		coulomb_basis.N_Gy_T = 2*N_C_Gy - 1;
		coulomb_basis.N_Gz_T = 2*N_C_Gz - 1;

		//coulomb_basis.Gx = ( double* ) malloc( coulomb_basis.N_Gx_T*sizeof( double ) );
		//coulomb_basis.Gy = ( double* ) malloc( coulomb_basis.N_Gy_T*sizeof( double ) );
		//coulomb_basis.Gz = ( double* ) malloc( coulomb_basis.N_Gz_T*sizeof( double ) );

		local_mem += (coulomb_basis.N_Gx_T + coulomb_basis.N_Gy_T + coulomb_basis.N_Gz_T )*sizeof( double ); // Adding memory to store Coulomb supercell basis set
			

		coulomb_matrix_elements( basis, coulomb_basis, supercell, parameters, N_states, eigenstate_vb, eigenstate_cb, CME, coulomb_region_type );

		
		component_memory_usage.CME_mem = local_mem_peak;


		if ( myid == 0 )
		{

			//write_coulomb_matrix_elements( 2*N_states[0], 2*N_states[1], CME );

			//free ( CME.ee );
			//free ( CME.hh );
			//free ( CME.eh_direct );
			//free ( CME.eh_exchange );

			local_mem -= ( N_CME_ee + N_CME_hh + 2*N_CME_eh )*sizeof( complex double );

		}

		//free( coulomb_basis.Gx );
		//free( coulomb_basis.Gy );
		//free( coulomb_basis.Gz );

		local_mem -= (coulomb_basis.N_Gx_T + coulomb_basis.N_Gy_T + coulomb_basis.N_Gz_T )*sizeof( double );

	}



	/////////////////////////////////
	//// Write code memory usage ////
	/////////////////////////////////

	write_memory_usage();
        
        
	// Moving output files to graph folder
	
	if ( myid == 0 )
	{

		char copy_basis_set[80] = "mv basis_set.dat "; // Copies basis_set.dat to graphs folder 
		strcat( copy_basis_set, graph_folder );
		system( copy_basis_set );

		char copy_memory_usage[80] = "mv memory_usage.dat "; // Copies basis_set.dat to graphs folder 
		strcat( copy_memory_usage, graph_folder );
		system( copy_memory_usage );
	
	}

	if ( myid == 0 )
	{

		printf( "Done.\n");

	}
	
	///////////////////////////////////////////
	//// Free dynamically allocated memory ////
	///////////////////////////////////////////

	free( directions );

	free( basis.Gx   );
	free( basis.Gy   );
	free( basis.Gz   );
	local_mem -= ( basis.N_Gx_T + basis.N_Gy_T + basis.N_Gz_T ) *sizeof( double );

	free( supercell  );
	free( parameters );
	
	for( int idx_state = 0; idx_state < 2*N_states[0]; idx_state++ )
	{

		//free( eigenstate_vb[idx_state].fourier );
		//free( eigenstate_vb[idx_state].bloch );

	}

	for( int idx_state = 0; idx_state < 2*N_states[1]; idx_state++ )
	{

		//free( eigenstate_cb[idx_state].fourier );
		//free( eigenstate_cb[idx_state].bloch );

	}

	local_mem -= ( 2*N_states[0] + 2*N_states[1] )*( na_rows ) *sizeof( complex double );
	
	free( eigenstate_cb );
	free( eigenstate_vb );
	
	//free( momentum_matrix_elements );

	local_mem -= ( Num_Polarisations*2*Num_CB_States*2*Num_VB_States )*sizeof( complex double );

	//printf( "Local mem end = %lld\n", local_mem ); // Check that memory equals zero at end
	
	BLACS_GRIDEXIT(&my_blacs_ctxt);
	MPI_Finalize();

	return 0;

}
