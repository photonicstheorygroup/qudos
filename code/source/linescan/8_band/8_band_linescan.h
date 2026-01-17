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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <yaml.h>
#include <malloc.h>
#include <math.h>
#include <complex.h>
#include <sys/time.h>
#include <time.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
typedef double _Complex DCOMPLEX;
#define MKL_Complex16 DCOMPLEX
#include <mkl.h>
#include <mkl_dfti.h>

///////////////////
//// Constants ////
///////////////////

#define N_BAND		8	    // # bands in k.p Hamiltonian
#define HBSQM 		7.61996  // ( ( Reduced Planck constant )^2 )/( Free electron mass  ) [eV (Angstrom)^2]
#define ESQEPS	  	180.9513 // ( ( Electronic charge )^2       )/( Vacuum permittivity ) [eV  Angstrom   ]

#define N_PARAM		36		// Number of input parameters

#define N_DIRECT	7        // Maximum number of directions along which to compute real space linescans
#define N_DIRECT_CS	8        // Maximum number of normal directions of which to compute real space cross sections

#define DELTA		2.0e-03  // Gaussian linewidth for density of states calculations
#define EPS 		1.0e-08  // "Small" number used for checking if variables are equal to some given number

#define EPS_ABS	1.0e-04
#define EPS_REL	1.0e-04
#define LIMIT		1000


//////////////////////////
//// Input parameters ////
//////////////////////////

extern char *endptr; // Used for reading in arguments to main( ) from the command line

extern int   N_shape;   // No. of shapes in supercell, including supercell shape
extern int   isotropic; // = 1 for isotropic approximation to strain tensor, = 0 for full anisotropic strain tensor
extern int   softening; // = 1 to soften interfaces via Gaussian convolution, = 0 for "hard" interfaces
extern int piezo_order; // = 0 for no piezoelectricity, = 1 for 1st order only, = 2 for 1st and 2nd order, = 3 for 2nd order only
extern int 	include_shape_elastic_constants;  // 0 for matrix material elastic constants, 1 for position-dependent (shape) elastic constants
extern int 	shape_elastic_constants_order;    // order of correction to strain tensor for shape elastic constants (suggest to = 2)

extern double sigma_softening; // Standard deviation of convoluted Gaussian

extern char graph_folder[80]; // Stores output file directory



/////////////////
//// Structs ////
/////////////////

struct wave_vectors { // Stores plane wave basis set information

	int   N_Gx, N_Gy, N_Gz;			// No. of plane waves for envelope functions along x, y, and, z, respectively
	int   N_Gx_T, N_Gy_T, N_Gz_T;	// No. of plane waves for Hamiltonian matrix element parameters along x, y, and, z, respectively, equal to 2*N_G(x/y/z) - 1
	double *Gx,  *Gy,  *Gz;			// Plane wave basis set of size N_G(x/y/z)_T

};

typedef struct characteristic_functions { complex double *chi; } _characteristic_functions;	// Stores shape characteristic function Fourier transforms
typedef struct characteristic_functions _characteristic_functions;

struct strain_tensor {	// Stores strain tensor Fourier transforms

	complex double *xx, *yy, *zz, *yz, *xz, *xy;

};

typedef struct convolution_char_strain {	// Stores convolutions between Fourier transforms of shape characteristic functions and strain tensor components

	complex double *char_xx, *char_yy, *char_zz, *char_yz, *char_xz, *char_xy, *char_hydrostatic, *char_biaxial, *xx_yz, *yy_yz, *zz_yz, *xy_xz, *yy_xz, *xx_xz, *zz_xz, *xy_yz, *zz_xy, *xx_xy, *yy_xy, *xz_yz, *char_xx_yz, *char_yy_yz, *char_zz_yz, *char_xy_xz, *char_yy_xz, *char_xx_xz, *char_zz_xz, *char_xy_yz, *char_zz_xy, *char_xx_xy, *char_yy_xy, *char_xz_yz;

} _convolution_char_strain;

struct piezo_potential {	// Stores Fourier transform of piezoelectric potential

	complex double *pz;

};

struct parameters {	// Stores parameters used for integrals calculated in some shape chactersitic function Fourier transforms

	double alpha, beta, R1, R2, f;

};


typedef struct supercell_geometry { int shape, embedded, material, material_type, RL1, RL2, RL3; double x, y, Ox, Oy, Oz, L1, L2, L3, L4, L5, SL1, SL2, SL3, alpha, beta, gamma; } _supercell_geometry; // Stores supercell parameters
typedef struct supercell_geometry _supercell_geometry;

typedef struct material_params { double  a, C11, C12, C44, e14, B114, B124, B156, Eg, Dso, Ep, P, sc, g1, g2, g3, ac, av, b, d, Ev_avg, E_so, E_lh, E_hh, E_cb, nr; } _material_params;	// Stores material parameters
typedef struct material_params _material_params;

typedef struct band_edge_energies { double *E_so, *E_lh, *E_hh, *E_cb; } _band_edge_energies;	// Stores band-edge energies
typedef struct band_edge_energies _band_edge_energies;

typedef struct interpolate_parameters {	double P_AB, P_AC, P_AD, P_BC, P_BD, P_CD, B_ABC, B_ABD, B_ACD, B_BCD, B_ABCD; } _interpolate_parameters; // Stores material interpolation parameters
typedef struct interpolate_parameters _interpolate_parameters;


/////////////////////////////
//// Function prototypes ////
/////////////////////////////

// Reads calculation parameters
#ifndef READ_PARAMETERS
#define READ_PARAMETERS

void read_parameters( int *N_shape, int *N_Gx, int *N_Gy, int *N_Gz, int *strain_region_type, int *N_S_Gx, int *N_S_Gy, int *N_S_Gz, double *T, int *piezo_order, double *Fx, double *Fy, double *Fz, int *N_grid, int *isotropic, int *include_shape_elastic_constants, int *shape_elastic_constants_order, int *softening, double *sigma_softening, int *directions, int* cross_sections, int *N_grid_cross_sections, int *directions_cross_sections, int *input_check );

#endif


// Reads command line inputs
#ifndef READ_COMMAND_LINE_INPUT
#define READ_COMMAND_LINE_INPUT

void read_command_line_input( int argc, char *argv[], int *N_shape, int *N_Gx, int *N_Gy, int *N_Gz, int *strain_region_type, int *N_S_Gx, int *N_S_Gy, int *N_S_Gz, double *T, int *piezo_order, double *Fx, double *Fy, double *Fz, int *N_grid, int *isotropic, int *include_shape_elastic_constants, int *shape_elastic_constants_order, int *softening, double *sigma_softening, int *directions, int* cross_sections, int *N_grid_cross_sections, int *directions_cross_sections, int *input_check );

#endif


// Reads supercell.in file
#ifndef READ_SUPERCELL
#define READ_SUPERCELL

void read_supercell( int N_shapes, int strain_region_type, struct supercell_geometry *geometry );

#endif


// Computes basis set
#ifndef BASIS_SET
#define BASIS_SET

void basis_set( struct wave_vectors basis, struct supercell_geometry *shape_info, int strain_region );

#endif


// Writes basis set to file
#ifndef WRITE_BASIS_SET
#define WRITE_BASIS_SET

void write_basis_set( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell );

#endif


// Calculates alloy material parameters
#ifndef MATERIAL_PARAMETERS
#define MATERIAL_PARAMETERS

void material_parameters( int N_shape, double T, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculates quaternary material parameter
#ifndef INTERPOLATE_QUATERNARY
#define INTERPOLATE_QUATERNARY

double interpolate_quaternary( int bowing, struct supercell_geometry supercell, struct interpolate_parameters alloy_parameters );

#endif


// Interface interpolation in the case of Gaussian softening
#ifndef INTERPOLATE_INTERFACE
#define INTERPOLATE_INTERFACE

void interpolate_interface( int N_shape, double *shape_frac, struct supercell_geometry *supercell, struct material_params *parameters, struct material_params *interface_params );

#endif


// Writes supercell parameters to file
#ifndef WRITE_SUPERCELL
#define WRITE_SUPERCELL

void write_supercell( int N_shape, double T, double Fx, double Fy, double Fz, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef KRONECKER_DELTA
#define KRONECKER_DELTA

double kronecker_delta( int idx1, int idx2 );

#endif


// Calculates Fourier transform of shape characteristic function
#ifndef CHARACTERISTIC_FUNCTION_FOURIER
#define CHARACTERISTIC_FUNCTION_FOURIER

complex double characteristic_function_fourier( double Gx, double Gy, double Gz, struct supercell_geometry *supercell, int idx_shape, int strain_region );

#endif


// Calculates Fourier transform of pyramid segment characteristic function
#ifndef CHARACTERISTIC_FUNCTION_FOURIER_PYRAMID
#define CHARACTERISTIC_FUNCTION_FOURIER_PYRAMID

complex double characteristic_function_fourier_pyramid_segment( double Gx0, double Gy0, double Gz0, int idx_side, int N_side, struct supercell_geometry *supercell, int idx_shape );

#endif


// Calculates Fourier transform of prism segment characteristic function
#ifndef CHARACTERISTIC_FUNCTION_FOURIER_PRISM
#define CHARACTERISTIC_FUNCTION_FOURIER_PRISM

complex double characteristic_function_fourier_prism_segment( double Gx0, double Gy0, double Gz0, int idx_side, int N_side, struct supercell_geometry *supercell, int idx_shape );

#endif


// Integration by trapezoidal rule
#ifndef TRAPEZOIDAL_RULE
#define TRAPEZOIDAL_RULE

double trapezoidal_rule( int N_x, double *x, double *f );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I_E0
#define INTEGRAL_I_E0

complex double I_e0( double a, double G );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_CHI_1
#define INTEGRAL_CHI_1

complex double chi_1( double Gx, double Gy, double Gz, double Lx, double Ly, double h );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I_J0
#define INTEGRAL_I_J0

double I_J0( double alpha, double beta );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I_J0_GSL
#define INTEGRAL_I_J0_GSL

double I_J0_GSL( double alpha, double beta );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I_J0_GSL_INTEGRAND
#define INTEGRAL_I_J0_GSL_INTEGRAND

double integrand_I_J0( double x, void *params );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2
#define INTEGRAL_I2

complex double I_2( double alpha, double beta );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_GSL
#define INTEGRAL_I2_GSL

complex double I_2_GSL(double alpha, double beta );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_GSL_INTEGRAND_REAL
#define INTEGRAL_I2_GSL_INTEGRAND_REAL

double integrand_I_2_real( double x, void *params );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_GSL_INTEGRAND_IMAG
#define INTEGRAL_I2_GSL_INTEGRAND_IMAG

double integrand_I_2_imag( double x, void *params );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_GZ
#define INTEGRAL_I2_GZ

complex double I_2_Gz( double alpha);

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_GZ_GSL
#define INTEGRAL_I2_GZ_GSL

complex double I_2_Gz_GSL( double alpha );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_GZ_GSL_INTEGRAND_REAL
#define INTEGRAL_I2_GZ_GSL_INTEGRAND_REAL

double integrand_I_2_Gz_real( double x, void *params );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_TC
#define INTEGRAL_I2_TC

complex double I_2_TC( double alpha, double R1, double R2, double beta );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_TC_GSL
#define INTEGRAL_I2_TC_GSL

complex double I_2_TC_GSL(double alpha, double R1, double R2, double beta );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_TC_GSL_INTEGRAND_REAL
#define INTEGRAL_I2_TC_GSL_INTEGRAND_REAL

double integrand_I_2_TC_real( double x, void *params );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_TC_GSL_INTEGRAND_IMAG
#define INTEGRAL_I2_TC_GSL_INTEGRAND_IMAG

double integrand_I_2_TC_imag( double x, void *params );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_TC_GZ
#define INTEGRAL_I2_TC_GZ

complex double I_2_TC_Gz( double alpha, double R1, double R2);

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_TC_GZ_GSL
#define INTEGRAL_I2_TC_GZ_GSL

complex double I_2_TC_Gz_GSL(double alpha, double R1, double R2);

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_TC_GZ_GSL_INTEGRAND_REAL
#define INTEGRAL_I2_TC_GZ_GSL_INTEGRAND_REAL

double integrand_I_2_TC_Gz_real( double x, void *params );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_ETC
#define INTEGRAL_I2_ETC

complex double I_2_ETC( double alpha, double beta, double f);

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_ETC_GSL
#define INTEGRAL_I2_ETC_GSL

complex double I_2_ETC_GSL(double alpha, double beta, double f );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_ETC_GSL_INTEGRAND_REAL
#define INTEGRAL_I2_ETC_GSL_INTEGRAND_REAL

double integrand_I_2_ETC_real( double x, void *params );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_ETC_GSL_INTEGRAND_IMAG
#define INTEGRAL_I2_ETC_GSL_INTEGRAND_IMAG

double integrand_I_2_ETC_imag( double x, void *params );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_ETC_GZ
#define INTEGRAL_I2_ETC_GZ

complex double I_2_ETC_Gz( double alpha, double f);

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_ETC_GZ_GSL
#define INTEGRAL_I2_ETC_GZ_GSL

complex double I_2_ETC_Gz_GSL( double alpha, double f );

#endif


// Shape characeteristic function Fourier transform integral
#ifndef INTEGRAL_I2_ETC_GZ_GSL_INTEGRAND_REAL
#define INTEGRAL_I2_ETC_GZ_GSL_INTEGRAND_REAL

double integrand_I_2_ETC_Gz_real( double x, void *params );

#endif


// Calculate shape characteristic function Fourier transforms
#ifndef CHARACTERISTIC_FUNCTIONS_RECIPROCAL_SPACE
#define CHARACTERISTIC_FUNCTIONS_RECIPROCAL_SPACE

void characteristic_functions_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct characteristic_functions *cf_G, int strain_region );

#endif


// Calculate strain tensor components in reciprocal space for strain supercell
#ifndef STRAIN_REGION_RECIPROCAL_SPACE
#define STRAIN_REGION_RECIPROCAL_SPACE

void strain_region_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct strain_tensor strain_G, struct characteristic_functions *characteristic_function_G );

#endif


// Downsample strain tensor components for electronic structure supercell (free supercell lengths)
#ifndef STRAIN_RECIPROCAL_SPACE
#define STRAIN_RECIPROCAL_SPACE

void strain_reciprocal_space( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell, struct strain_tensor strain_G, struct strain_tensor strain_G_expanded );

#endif


// Downsample strain tensor components for electronic structure supercell (commensurate supercell lengths)
#ifndef STRAIN_RECIPROCAL_SPACE_FFT
#define STRAIN_RECIPROCAL_SPACE_FFT

void strain_reciprocal_space_fft( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell, struct strain_tensor strain_G, struct strain_tensor strain_G_expanded, struct characteristic_functions *characteristic_function_strain_G );

#endif


// Calculate strain tensor component fourier coefficient
#ifndef STRAIN_FOURIER_COEFFICIENT
#define STRAIN_FOURIER_COEFFICIENT

complex double strain_fourier_coefficient( int component, int isotropic, int idx_shape, int idx_basis, int G_prime, double Gx, double Gy, double Gz, struct supercell_geometry *supercell, struct material_params *parameters, struct characteristic_functions characteristic_function_G, int correction_order, struct convolution_char_strain *convolutions_char_strain_correction );

#endif


// Green's function tensor for strain calculation
#ifndef GREENS_FUNCTION_TENSOR
#define GREENS_FUNCTION_TENSOR

double greens_function_tensor( struct material_params *parameters, int component, double Gx, double Gy, double Gz, double G, double sum );

#endif


// Calculate piezoelectric potential in reciprocal space for strain supercell
#ifndef PIEZO_REGION_RECIPROCAL_SPACE
#define PIEZO_REGION_RECIPROCAL_SPACE

void piezo_region_reciprocal_space( int piezo_order, int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct characteristic_functions *characteristic_function_G, struct strain_tensor strain_G, complex double *piezo_G );

#endif


// Downsample piezoelectric potential for electronic structure supercell (free supercell lengths)
#ifndef PIEZO_RECIPROCAL_SPACE
#define PIEZO_RECIPROCAL_SPACE

void piezo_reciprocal_space( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell, complex double *phi_piezo_G, complex double *phi_piezo_region_G );

#endif


// Downsample piezoelectric potential for electronic structure supercell (commensurate supercell lengths)
#ifndef PIEZO_RECIPROCAL_SPACE_FFT
#define PIEZO_RECIPROCAL_SPACE_FFT

void piezo_reciprocal_space_fft( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell, complex double *phi_piezo_G, complex double *phi_piezo_region_G, struct characteristic_functions *characteristic_function_strain_G );

#endif


// Calculate piezoelectric potential Fourier coefficient
#ifndef PIEZO_FOURIER_COEFFICIENT
#define PIEZO_FOURIER_COEFFICIENT

complex double piezo_fourier_coefficient( int piezo_order, int idx_shape, int idx_Gx, int idx_Gy, int idx_Gz, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct characteristic_functions *characteristic_function_G, struct strain_tensor strain_G, struct convolution_char_strain convolutions_char_strain_region );

#endif


// Calculate convolutions required for position-dependent elastic constants in strain calculation
#ifndef CONVOLUTIONS_STRAIN_CORRECTION
#define CONVOLUTIONS_STRAIN_CORRECTION

void convolutions_strain_correction( int idx_shape, struct wave_vectors basis, struct characteristic_functions *characteristic_function_strain_G, struct strain_tensor *strain_correction_region_G, struct convolution_char_strain *convolutions_char_strain_correction );

#endif


// Set up linescan grid
#ifndef REAL_SPACE_GRID_LINESCAN
#define REAL_SPACE_GRID_LINESCAN

void real_space_grid_linescan( int direction, int N_grid, struct supercell_geometry *supercell, double *x, double *y, double *z );

#endif


// Set up cross-section grid
#ifndef REAL_SPACE_GRID_CROSS_SECTION
#define REAL_SPACE_GRID_CROSS_SECTION

void real_space_grid_cross_section( int direction, int N_grid, struct supercell_geometry *supercell, double *x, double *y, double *z );

#endif


// Calculate strain in real space
#ifndef STRAIN_REAL_SPACE
#define STRAIN_REAL_SPACE

void strain_real_space( int N_shape, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct strain_tensor strain_G, struct strain_tensor strain_r );

#endif


// Calculate piezoelectric potential in real space
#ifndef PIEZO_REAL_SPACE
#define PIEZO_REAL_SPACE

void piezo_real_space( int N_shape, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, complex double *piezo_G, complex double *piezo_r );

#endif


// Write strain linescans to file
#ifndef WRITE_STRAIN_LINESCAN
#define WRITE_STRAIN_LINESCAN

void write_strain_linescan( int direction, int N_grid, double *x, double *y, double *z, struct strain_tensor strain_r );

#endif


// Write strain cross-sections to file
#ifndef WRITE_STRAIN_CROSS_SECTION
#define WRITE_STRAIN_CROSS_SECTION

void write_strain_cross_section( int direction, int N_grid, double *x, double *y, double *z, struct strain_tensor strain_r );

#endif


// Write piezo linescans to file
#ifndef WRITE_PIEZO_LINESCAN
#define WRITE_PIEZO_LINESCAN

void write_piezo_linescan( int direction, int N_grid, double *x, double *y, double *z, int piezo_order, complex double *piezo_r );

#endif


// Write piezo cross-sections to file
#ifndef WRITE_PIEZO_CROSS_SECTION
#define WRITE_PIEZO_CROSS_SECTION

void write_piezo_cross_section( int direction, int N_grid, double *x, double *y, double *z, int piezo_order, complex double *piezo_r );

#endif


// Calculate band-edge energies in real space
#ifndef BAND_EDGE_ENERGIES_REAL_SPACE
#define BAND_EDGE_ENERGIES_REAL_SPACE

void band_edge_energies_real_space( int N_grid, int N_shape, double *x, double *y, double *z, double Fx, double Fy, double Fz, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct characteristic_functions *characteristic_function_G, struct strain_tensor strain_r, complex double *piezo_r, struct band_edge_energies band_edges );

#endif


// Determine shape real space grid point resides within
#ifndef DETERMINE_IDX_SHAPE
#define DETERMINE_IDX_SHAPE

int determine_idx_shape( int N_shape, double x, double y, double z, struct supercell_geometry *supercell );

#endif


// Real space characteristic functions
#ifndef CHARACTERTIC_FUNCTION_REAL_SPACE
#define CHARACTERTIC_FUNCTION_REAL_SPACE

double characteristic_function_real_space( double x, double y, double z, struct wave_vectors basis, struct supercell_geometry *supercell, struct characteristic_functions characteristic_function_G );

#endif


// Calculate bulk Hamiltonian eigenstates in real space
#ifndef EIGENSTATES_REAL_SPACE
#define EIGENSTATES_REAL_SPACE

void eigenstates_real_space( int idx_grid, int idx_shape, double x, double y, double z, double Fx, double Fy, double Fz, struct material_params parameters, struct strain_tensor strain_r, complex double *piezo_r, struct band_edge_energies band_edges );

#endif


// Calculate bulk Hamiltonian matrix elements
#ifndef HAMILTONIAN_8_BAND_BULK
#define HAMILTONIAN_8_BAND_BULK

complex double hamiltonian_8_band_bulk( int idx_row, int idx_col, double x, double y, double z, struct material_params parameters, struct strain_tensor epsilon, complex double phi_piezo, complex double phi_F );

#endif


// Write band-edge energy linescans
#ifndef WRITE_BAND_EDGE_ENERGIES_LINESCAN
#define WRITE_BAND_EDGE_ENERGIES_LINESCAN

void write_band_edge_energies_linescan( int direction, int N_grid, double *x, double *y, double *z, struct band_edge_energies band_edges );

#endif


// Write band-edge energy cross-sections
#ifndef WRITE_BAND_EDGE_ENERGIES_CROSS_SECTION
#define WRITE_BAND_EDGE_ENERGIES_CROSS_SECTION

void write_band_edge_energies_cross_section( int direction, int N_grid, double *x, double *y, double *z, struct band_edge_energies band_edges );

#endif


#ifndef COPY_OUTPUT
#define COPY_OUTPUT

void copy_output( int type, int direction );

#endif


#ifndef COPY_OUTPUT_LINESCAN
#define COPY_OUTPUT_LINESCAN

void copy_output_linescan( int direction, char *input_name );

#endif


#ifndef COPY_OUTPUT_CROSS_SECTION
#define COPY_OUTPUT_CROSS_SECTION

void copy_output_cross_section( int direction, char *input_name );

#endif


// Check VSL function error status and abort if error detected
#ifndef CHECK_VSL
#define CHECK_VSL

void check_vsl( MKL_INT status, const char *func );

#endif
