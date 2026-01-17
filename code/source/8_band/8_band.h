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
#include <mpi.h>
#include <omp.h>
#include <mkl_blacs.h>
#include <mkl_scalapack.h>
#include <mkl_pblas.h>
#include <mkl_cdft.h>


///////////////////
//// Constants ////
///////////////////

#define N_BAND		8	    // # bands in k.p Hamiltonian
#define HBSQM 		7.61996  // ( ( Reduced Planck constant )^2 )/( Free electron mass  ) [eV (Angstrom)^2]
#define ESQEPS	  	180.9513 // ( ( Electronic charge )^2       )/( Vacuum permittivity ) [eV  Angstrom   ]

#define N_PARAM		53		// Number of input parameters

#define N_DIRECT	7        // Maximum number of directions along which to compute real space linescans
#define N_DIRECT_CS	8        // Maximum number of normal directions of which to compute real space cross sections

#define EPS 		1.0e-08  // "Small" number used for checking if variables are equal to some given number

#define EPS_ABS	1.0e-04 	// Absolute error for GSL QAG integration function
#define EPS_REL	1.0e-04 	// Relative error for GSL QAG integration function
#define LIMIT		1000	// Number of subintervaks for GSL QAG integration function


//////////////////////////
//// Input parameters ////
//////////////////////////

extern char *endptr; // Used for reading in arguments to main( ) from the command line


//////////////////////////
//// Global variables ////
//////////////////////////

extern int  N_shape;							// No. of shapes in supercell, including supercell shape
extern int  isotropic;							// Isotropic approximation for strain - 0 (1): anisotropic (isotropic)
extern int  softening;							// Include Gaussian softening of material parameter profiles at interfaces. Note: not perfect, does not accurately incorporate bowing terms.
extern int 	include_strain;						// Include strain - 0 (1): exclude (include)
extern int 	include_piezo;						// Include piezoelectric potential - 0 to exclude, 1 for 1st order, 2 for 1st and 2nd order, 3 for 2nd order only
extern int 	include_shape_elastic_constants;	// 0 for matrix material elastic constants, 1 for position-dependent (shape) elastic constants
extern int 	shape_elastic_constants_order;		// order of correction to strain tensor for shape elastic constants (suggest to = 2)
extern int 	optimise;							// Set to 1 to allow eigensolver to determine optimal solver workspace. WARNING: Results in large RAM usage by each process. Workspace memory usage in total is equal to roughly 4 times that used to store the Hamiltonian matrix (increases minimally with number of MPI processes). Set to 0 to assign the minimum workspace size required (usually negligible relative to memory used to store matrix elements). Optimal memory usage results in ~3 to 4 times faster eigensolver calculation.
extern int 	include_overlaps;					// 0 for no overlaps, 1 for overlaps
extern int 	include_localisation;				// 0 for no carrier shape localisation, 1 for carrier shape localisation
extern int 	include_MME;						// 0 for no momentum matrix elements, 1 for momentum matrix elements
extern int 	include_CME;						// 0 for no coulomb matrix elements, 1 for coulomb matrix elements with no corrections, 2 for coulomb matrix elements with monopole corrections, 3 for coulomb matrix elements with dipole corrections, 4 for coulomb matrix elements with quadrupole corrections

extern double sigma_softening;	// Standard deviation of Gaussian used for interface softening

extern MKL_INT N_H;			// Order of Hamiltonian matrix

extern char graph_folder[80];	// Stores output file directory

extern int PACKED;			// Enable packed matrix storage for eigensolver to decrease memory usage at cost of computation time - 0 (1): disable (enable). Note: only supported for single process mode, using OpenMP only.


///////////////////////////////
//// MPI & BLACS variables ////
///////////////////////////////

extern int myid, nprocs;  	// MPI process ID; No. MPI processes
extern MKL_INT na_cols, na_rows, na_colsz, na_rowsz;	// No. local Hamiltonian matrix columns; No. local Hamiltonian matrix rows; No. local eigenvector matrix columns; No. local eigenvector matrix rows
extern MKL_INT np_cols, np_rows;	// No. BLACS process columns; No. BLACS process rows 
extern MKL_INT my_prow, my_pcol;	// BLACS process row of my MPI process; BLACS process column of my MPI process

extern MKL_INT my_blacs_ctxt, sc_desc[9], sc_descz[9], info;	// BLACS system context; Array descriptor of distributed Hamiltonian matrix; Array descriptor of distributed eigenvector matrix

extern MKL_INT block_size;	// Blocking factor used to distribute matrix rows and columns among processes


/////////////////
//// Structs ////
/////////////////

struct wave_vectors { // Stores plane wave basis set information

	int   N_Gx, N_Gy, N_Gz; 		// No. of plane waves for envelope functions along x, y, and, z, respectively
	int   N_Gx_T, N_Gy_T, N_Gz_T; 	// No. of plane waves for Hamiltonian matrix element parameters along x, y, and, z, respectively, equal to 2*N_G(x/y/z) - 1
	double *Gx,  *Gy,  *Gz;			// Plane wave basis set of size N_G(x/y/z)_T
	double Kx, Ky, Kz;				// Supercell wave vector component for each direction (for miniband dispersion)

};

typedef struct characteristic_functions { complex double *chi; } _characteristic_functions; // Stores shape characteristic function Fourier transforms
typedef struct characteristic_functions _characteristic_functions;


typedef struct vector { double x, y, z; } _vector; // Vector

typedef struct vector_complex { complex double x, y, z; } _vector_complex; // Complex vector

typedef struct strain_tensor { // Stores strain tensor Fourier transforms

	complex double *xx, *yy, *zz, *yz, *xz, *xy;

} _strain_tensor;

typedef struct convolution_char_strain { // Stores convolutions between Fourier transforms of shape characteristic functions and strain tensor components

	complex double *char_xx, *char_yy, *char_zz, *char_yz, *char_xz, *char_xy, *char_hydrostatic, *char_biaxial, *xx_yz, *yy_yz, *zz_yz, *xy_xz, *yy_xz, *xx_xz, *zz_xz, *xy_yz, *zz_xy, *xx_xy, *yy_xy, *xz_yz, *char_xx_yz, *char_yy_yz, *char_zz_yz, *char_xy_xz, *char_yy_xz, *char_xx_xz, *char_zz_xz, *char_xy_yz, *char_zz_xy, *char_xx_xy, *char_yy_xy, *char_xz_yz;

} _convolution_char_strain;

typedef struct piezo_potential { // Stores Fourier transform of piezoelectric potential

	complex double *pz;

} _piezo_potential;

typedef struct convolution_char_piezo { // Stores convolution between shape characteristic function and piezoelectric potential

	complex double *char_piezo;

} _convolution_char_piezo;

typedef struct convolution_state { // Stores convolution between two state envelope functions

	complex double *conv_state;

} _convolution_state;

typedef struct convolution_char_state { // Stores convolution between state envelope functions and shape characteristic function

	complex double *conv_char_state;

} _convolution_char_state;

typedef struct coulomb_matrix_element { // Stores Coulomb matrix elements

	complex double *ee, *hh, *eh_direct, *eh_exchange;

} _coulomb_matrix_element;

typedef struct coulomb_corrections { // Stores parameters used for Coulomb matrix element multipole expansion corrections

	double a_mad;
	complex double *dx, *dy, *dz, *Q;
	complex double *Ix, *Iy, *Iz, *J;

} _coulomb_corrections;

struct parameters { // Stores parameters used for integrals calculated in some shape chactersitic function Fourier transforms

	double alpha, beta, R1, R2, f;

};

typedef struct timing { // Stores process calculation times of code stages

	double char_s_grid, char_es_grid, strain_s_grid, piezo_s_grid, strain_es_grid, piezo_es_grid, conv, ham_elements, ham_diag, eig, bloch, overlaps, car_loc, mme, linescan, cross_section, real_3D, cme_B, cme_I, cme_total, total;

} _timing;


typedef struct supercell_geometry { int shape, embedded, material, material_type, RL1, RL2, RL3, RCL1, RCL2, RCL3; double x, y, Ox, Oy, Oz, L1, L2, L3, L4, L5, SL1, SL2, SL3, CL1, CL2, CL3, alpha, beta, gamma; } _supercell_geometry; // Stores supercell parameters
typedef struct supercell_geometry _supercell_geometry;

typedef struct material_params { double  a, C11, C12, C44, e14, B114, B124, B156, Eg, Dso, Ep, P, sc, g1, g2, g3, ac, av, b, d, Ev_avg, E_so, E_lh, E_hh, E_cb, nr; } _material_params; // Stores material parameters
typedef struct material_params _material_params;

typedef struct band_edge_energies { double *E_so, *E_lh, *E_hh, *E_cb; } _band_edge_energies; // Stores band-edge energies
typedef struct band_edge_energies _band_edge_energies;

typedef struct interpolate_parameters {	double P_AB, P_AC, P_AD, P_BC, P_BD, P_CD, B_ABC, B_ABD, B_ACD, B_BCD, B_ABCD; } _interpolate_parameters; // Stores material interpolation parameters
typedef struct interpolate_parameters _interpolate_parameters;

typedef struct eigenstate { int idx_state; double E, *bloch, *probability_1, *probability_2, *probability_3, *probability_4, *probability_5, *probability_6, *probability_7, *probability_8;
                            complex double *fourier, *F_1, *F_2, *F_3, *F_4, *F_5, *F_6, *F_7, *F_8; } _eigenstate; // Stores eigenstate energies, Bloch characters, Fourier coefficients, and probability densities 
typedef struct eigenstate _eigenstate;


////////////////////////////////////////////////////////////////////////////////////
//// Global characteristic function, strain, and piezoelectric potential arrays ////
////////////////////////////////////////////////////////////////////////////////////

extern struct characteristic_functions *characteristic_function_G; 			// Stores shape characteristic function Fourier transforms for electronic structure supercell
extern struct characteristic_functions *characteristic_function_strain_G; 	// Stores shape characteristic function Fourier transforms for strain supercell
extern struct strain_tensor strain_G;			// Stores strain tensor Fourier transforms for electronic structure supercell
extern struct strain_tensor strain_region_G;	// Stores strain tensor Fourier transforms for strain supercell
extern struct convolution_char_strain *convolutions_char_strain; 	// Stores convolutions between strain tensor components and shape characteristic functions
extern struct piezo_potential phi_piezo_G;	// Stores piezoelectric potential for electronic structure supercell
extern struct piezo_potential phi_piezo_region_G;	// Stores piezoelectric potential for strain supercell
extern struct timing component_run_times;	// Stores code stage run times

////////////////////////////////////////////////////////////
//// Global parameters used for momentum matrix elements////
////////////////////////////////////////////////////////////

extern complex double *dH_dk; // Stores local part of Hamiltonian matrix derivative
extern complex double *yV;	  // Stores local part of intermediate matrix-vector product as part of vector-matrix-vector product calculation
extern int Num_Possible_Polarisations;	// Number of polarisations for calculation
extern int Polarisation_Indicator;		// Indicator for current polarisation


/////////////////////////////
//// Function prototypes ////
/////////////////////////////

// Reads calculation parameters
#ifndef READ_PARAMETERS
#define READ_PARAMETERS

void read_parameters( int *N_shape, int *N_Gx, int *N_Gy, int *N_Gz, int *strain_region_type, int *N_S_Gx, int *N_S_Gy, int *N_S_Gz, int *N_vb, int *N_cb, double *T, double *Fx, double *Fy, double *Fz, double *Kx, double *Ky, double *Kz, int *N_grid, int *include_strain, int *isotropic, int *include_shape_elastic_constants, int *shape_elastic_constants_order, int *optimise, int *include_overlaps, int *include_localisation, int *include_piezo, int *include_MME, int *include_CME, int *coulomb_region_type, int *N_C_Gx, int *N_C_Gy, int *N_C_Gz, int *softening, double *sigma_softening, int* directions, int* cross_sections, int *N_grid_cross_sections, int *directions_cross_sections, int *enable_3D, int *N_grid_3D, int *input_check );

#endif


// Reads command line inputs
#ifndef READ_COMMAND_LINE_INPUT
#define READ_COMMAND_LINE_INPUT

void read_command_line_input( int argc, char *argv[], int *N_shape, int *N_Gx, int *N_Gy, int *N_Gz, int *strain_region_type, int *N_S_Gx, int *N_S_Gy, int *N_S_Gz, int *N_vb, int *N_cb, double *T, double *Fx, double *Fy, double *Fz, double *Kx, double *Ky, double *Kz, int *N_grid, int *include_strain, int *isotropic, int *include_shape_elastic_constants, int *shape_elastic_constants_order, int *optimise, int *include_overlaps, int *include_localisation, int *include_piezo, int *include_MME, int *include_CME, int *coulomb_region_type, int *N_C_Gx, int *N_C_Gy, int *N_C_Gz, int *softening, double *sigma_softening, int* directions, int* cross_sections, int *N_grid_cross_sections, int *directions_cross_sections, int *enable_3D, int *N_grid_3D, int *input_check );

#endif


// Reads supercell.in file
#ifndef READ_SUPERCELL
#define READ_SUPERCELL

void read_supercell( int N_shapes, int strain_region_type, int coulomb_region_type, struct supercell_geometry *geometry );

#endif


// Computes basis set
#ifndef BASIS_SET
#define BASIS_SET

void basis_set( struct wave_vectors basis, struct supercell_geometry *shape_info, int strain_region, int coulomb_region );

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

complex double characteristic_function_fourier( double Gx, double Gy, double Gz, struct supercell_geometry *supercell, int idx_shape, int supercell_type );

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

void characteristic_functions_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, int strain_region );

#endif


// Calculate electronic structure supercell inside Coulomb supercell characteristic function Fourier transform
#ifndef CHARACTERISTIC_FUNCTION_COULOMB_SUPERCELL
#define CHARACTERISTIC_FUNCTION_COULOMB_SUPERCELL

void characteristic_function_coulomb_supercell( struct wave_vectors coulomb_basis, struct supercell_geometry *supercell, struct characteristic_functions characteristic_function_coulomb_G );

#endif


// Calculate strain tensor components in reciprocal space for strain supercell
#ifndef STRAIN_REGION_RECIPROCAL_SPACE
#define STRAIN_REGION_RECIPROCAL_SPACE

void strain_region_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Downsample strain tensor components for electronic structure supercell (free supercell lengths)
#ifndef STRAIN_RECIPROCAL_SPACE
#define STRAIN_RECIPROCAL_SPACE

void strain_reciprocal_space( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell );

#endif


// Downsample strain tensor components for electronic structure supercell (commensurate supercell lengths)
#ifndef STRAIN_RECIPROCAL_SPACE_FFT
#define STRAIN_RECIPROCAL_SPACE_FFT

void strain_reciprocal_space_fft( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell );

#endif


// Calculate strain tensor component fourier coefficient
#ifndef STRAIN_FOURIER_COEFFICIENT
#define STRAIN_FOURIER_COEFFICIENT

complex double strain_fourier_coefficient( int component, int idx_shape, int idx_basis, double Gx, double Gy, double Gz, struct supercell_geometry *supercell, struct material_params *parameters, int correction_order, struct convolution_char_strain *convolutions_char_strain_correction );

#endif


// Green's function tensor for strain calculation
#ifndef GREENS_FUNCTION_TENSOR
#define GREENS_FUNCTION_TENSOR

double greens_function_tensor( struct material_params *parameters, int component, double Gx, double Gy, double Gz, double G, double sum );

#endif


// Calculate piezoelectric potential in reciprocal space for strain supercell
#ifndef PIEZO_REGION_RECIPROCAL_SPACE
#define PIEZO_REGION_RECIPROCAL_SPACE

void piezo_region_reciprocal_space( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Downsample piezoelectric potential for electronic structure supercell (free supercell lengths)
#ifndef PIEZO_RECIPROCAL_SPACE
#define PIEZO_RECIPROCAL_SPACE

void piezo_reciprocal_space( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell );

#endif


// Downsample piezoelectric potential for electronic structure supercell (commensurate supercell lengths)
#ifndef PIEZO_RECIPROCAL_SPACE_FFT
#define PIEZO_RECIPROCAL_SPACE_FFT

void piezo_reciprocal_space_fft( struct wave_vectors basis, struct wave_vectors strain_basis, struct supercell_geometry *supercell );

#endif


// Calculate piezoelectric potential Fourier coefficient
#ifndef PIEZO_FOURIER_COEFFICIENT
#define PIEZO_FOURIER_COEFFICIENT

complex double piezo_fourier_coefficient( int piezo_order, int idx_shape, int idx_Gx_global, int idx_Gy, int idx_Gz, int idx_Gx_local, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters,  struct strain_tensor strain_region_G, struct convolution_char_strain convolutions_char_strain_region_reduced_local );

#endif


// Calculate electric field potential Fourier coefficient
#ifndef FIELD_FOURIER_COEFFICIENT
#define FIELD_FOURIER_COEFFICIENT

complex double field_fourier_coefficient( struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis );

#endif


// Calculate convolutions between reciprocal space strain tensor components and shape characteristic functions
#ifndef CONVOLUTIONS
#define CONVOLUTIONS

void convolutions( int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell );

#endif


// Calculate convolutions required for position-dependent elastic constants in strain calculation
#ifndef CONVOLUTIONS_STRAIN_CORRECTION
#define CONVOLUTIONS_STRAIN_CORRECTION

void convolutions_strain_correction( int idx_shape, struct wave_vectors basis, struct strain_tensor *strain_correction_region_G, struct convolution_char_strain *convolutions_char_strain_correction );

#endif


// Calculate eigenstates 
#ifndef EIGENSTATES
#define EIGENSTATES

void eigenstates( int N_shape, int *N_states, struct vector F, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb );

#endif


// Calculate 8-band Hamiltonian matrix elements
#ifndef HAMILTONIAN_8_BAND
#define HAMILTONIAN_8_BAND

void hamiltonian_8_band( int N_shape, struct vector F, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, complex double* H );

#endif


// Calculate material parameter contribution to Hamiltonian matrix element
#ifndef PARAMETER_MATRIX_ELEMENT
#define PARAMETER_MATRIX_ELEMENT

complex double parameter_matrix_element( int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, double *parameter );

#endif


// Calculate E_CB Hamiltonian matrix element
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_E_CB
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_E_CB

complex double E_CB( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate E_HH Hamiltonian matrix element
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_E_HH
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_E_HH

complex double E_HH( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate E_LH Hamiltonian matrix element
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_E_LH
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_E_LH

complex double E_LH( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate E_SO Hamiltonian matrix element
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_E_SO
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_E_SO

complex double E_SO( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate Q Hamiltonian matrix element parameter
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_Q
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_Q

complex double Q( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate R Hamiltonian matrix element parameter
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_R
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_R

complex double R( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate R_conj Hamiltonian matrix element parameter
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_R_CONJ
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_R_CONJ

complex double R_conj( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate S Hamiltonian matrix element parameter
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_S
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_S

complex double S( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate S_conj Hamiltonian matrix element parameter
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_S_CONJ
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_S_CONJ

complex double S_conj( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate T_plus Hamiltonian matrix element parameter
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_T_PLUS
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_T_PLUS

complex double T_plus( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate T_minus Hamiltonian matrix element parameter
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_T_MINUS
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_T_MINUS

complex double T_minus( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate U Hamiltonian matrix element parameter
#ifndef HAMILTONIAN_8_BAND_MATRIX_ELEMENT_U
#define HAMILTONIAN_8_BAND_MATRIX_ELEMENT_U

complex double U( int N_shape, struct vector F, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Diagonalise Hamiltonian
#ifndef DIAGONALISE_HAMILTONIAN
#define DIAGONALISE_HAMILTONIAN

void diagonalise_hamiltonian( int *idx_eig, complex double *H, double *eigval, complex double *eigvec );

#endif


// Collect eigenstates after diagonalisation
#ifndef COLLECT_EIGENSTATES
#define COLLECT_EIGENSTATES

void collect_eigenstates( int *N_states, struct wave_vectors basis, double *eigval, complex double *eigvec, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb );

#endif


// Calculate eigenstate Bloch character
#ifndef BLOCH_CHARACTER
#define BLOCH_CHARACTER

void bloch_character( int band, int *N_states, struct wave_vectors basis, struct eigenstate *eigenstate );

#endif


// Calculate Hamiltonian derivate 8-band matrix
#ifndef HAMILTONIAN_DERIVATIVE_8_BAND
#define HAMILTONIAN_DERIVATIVE_8_BAND

void hamiltonian_derivative_8_band( int polarisation, int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, complex double* H );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_E_CB
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_E_CB

complex double dE_CB_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_E_HH
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_E_HH

complex double dE_HH_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_E_LH
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_E_LH

complex double dE_LH_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_E_SO
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_E_SO

complex double dE_SO_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_Q
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_Q

complex double dQ_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_R
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_R

complex double dR_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_R_CONJ
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_R_CONJ

complex double dR_conj_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_S
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_S

complex double dS_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_S_CONJ
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_S_CONJ

complex double dS_conj_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_T_PLUS
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_T_PLUS

complex double dT_plus_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_T_MINUS
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_T_MINUS

complex double dT_minus_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


#ifndef HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_U
#define HAMILTONIAN_DERIVATIVE_8_BAND_MATRIX_ELEMENT_U

complex double dU_dk( int polarisation, int N_shape, int *idx_G_1, int *idx_G_2, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters );

#endif


// Calculate momentum matrix element between two states
#ifndef MOMENTUM_MATRIX_ELEMENT
#define MOMENTUM_MATRIX_ELEMENT

complex double momentum_matrix_element( int polarisation, int N_shape, struct wave_vectors basis, struct supercell_geometry *supercell, struct material_params *parameters, struct eigenstate eigenstate_1, MKL_INT eigenstate_1_idx, struct eigenstate eigenstate_2, MKL_INT eigenstate_2_idx, int *N_states );

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


// Set up 3D grid
#ifndef REAL_SPACE_GRID_3D
#define REAL_SPACE_GRID_3D

void real_space_grid_3D( int N_grid, struct supercell_geometry *supercell, double *x, double *y, double *z );

#endif


// Calculate eigenstate real-space envelope functions and probability densities
#ifndef ENVELOPE_FUNCTIONS
#define ENVELOPE_FUNCTIONS

void envelope_functions( struct supercell_geometry *supercell, int *N_states, int type, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct eigenstate *eigenstate, int direction );

#endif


// Calculate spatial overlap between states
#ifndef SPATIAL_OVERLAP
#define SPATIAL_OVERLAP

complex double spatial_overlap( struct wave_vectors basis, struct eigenstate eigenstate_1, MKL_INT eigenstate_1_idx, struct eigenstate eigenstate_2, MKL_INT eigenstate_2_idx );

#endif


// Calculate eigenstate carrier localisation
#ifndef CARRIER_LOCALISATION
#define CARRIER_LOCALISATION

complex double carrier_localisation( int idx_shape, struct wave_vectors basis, struct supercell_geometry *supercell, int *N_states, struct characteristic_functions *characteristic_function_G_local, int *shape_indicator, int state_type, struct eigenstate *eigenstate_stored, int idx_state );

#endif


// Gaussian lineshape
#ifndef LINESHAPE_GAUSSIAN
#define LINESHAPE_GAUSSIAN

double lineshape_gaussian( double x, double delta );

#endif


// Calculate Coulomb matrix elements
#ifndef COULOMB_MATRIX_ELEMENTS
#define COULOMB_MATRIX_ELEMENTS

void coulomb_matrix_elements( struct wave_vectors basis, struct wave_vectors coulomb_basis, struct supercell_geometry *supercell, struct material_params *parameters, int *N_states, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb, struct coulomb_matrix_element CME, int coulomb_region_type );

#endif


// Redistribute an eigenstate across MPI processes
#ifndef REDISTRIBUTE_EIGENSTATE
#define REDISTRIBUTE_EIGENSTATE

void redistribute_eigenstate( struct wave_vectors basis, int *N_states, int eigenstate_type, int idx_state, struct eigenstate *eigenstate_stored, struct eigenstate eigenstate_interim );

#endif


// Redistribute eigenstates across MPI processes
#ifndef STATE_REDISTRIBUTION
#define STATE_REDISTRIBUTION

void state_redistribution( struct wave_vectors basis, struct supercell_geometry *supercell, int *N_states, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb, struct eigenstate *eigenstate_vb_padded, struct eigenstate *eigenstate_cb_padded, struct eigenstate *eigenstate_vb_padded_conj, struct eigenstate *eigenstate_cb_padded_conj, MKL_INT process_limit, MKL_LONG *sizes, MKL_LONG local_memory_size, MKL_LONG local_start, MKL_LONG local_rows );

#endif


// Calculate convolutions between states for B components of Coulomb matrix element calculation
#ifndef B_CONVOLUTIONS
#define B_CONVOLUTIONS

void B_convolutions( struct wave_vectors basis, struct wave_vectors coulomb_basis, struct supercell_geometry *supercell, int *N_states, struct eigenstate *eigenstate_vb_padded, struct eigenstate *eigenstate_cb_padded, struct eigenstate *eigenstate_vb_padded_conj, struct eigenstate *eigenstate_cb_padded_conj, MKL_INT process_num, MKL_INT process_limit, DFTI_DESCRIPTOR_DM_HANDLE task, MKL_LONG local_memory_size, MKL_LONG local_size, struct convolution_state *convolution_states, struct convolution_state *convolution_states_reversed, int *sendcounts, int *displs_send, int coulomb_region_type );

#endif


// Upsample B components for Colulomb supercell (free supercell lengths)
#ifndef B_COULOMB_SUPERCELL
#define B_COULOMB_SUPERCELL

void B_coulomb_supercell( struct wave_vectors basis, struct wave_vectors coulomb_basis, struct supercell_geometry *supercell, int local_size_coulomb, struct convolution_state convolution_states_interim, struct convolution_state convolution_states_interim_local, int *recvcounts_coulomb, int *displs_coulomb, struct convolution_state convolution_states_global );

#endif


// Upsample B components for Colulomb supercell (commensurate supercell lengths)
#ifndef B_COULOMB_SUPERCELL_FFT
#define B_COULOMB_SUPERCELL_FFT

void B_coulomb_supercell_fft( struct wave_vectors basis, struct wave_vectors coulomb_basis, struct supercell_geometry *supercell, MKL_INT process_limit_B_supercell, MPI_Comm MPI_SUB_COMM_B_supercell, MKL_LONG *sizes_B_supercell, DFTI_DESCRIPTOR_DM_HANDLE task_B_supercell, MKL_LONG local_start_B_supercell, MKL_LONG local_rows_B_supercell, MKL_LONG local_size_B_supercell, MKL_LONG local_memory_size_B_supercell, struct convolution_state convolution_states_interim, struct convolution_state convolution_states_interim_local, struct characteristic_functions characteristic_function_coulomb_G_local, struct convolution_char_state convolution_char_state_interim_local, MKL_INT local_start_reduced_B_supercell, MKL_INT local_rows_reduced_B_supercell, MKL_INT local_size_reduced_B_supercell, int *recvcounts_convolutions_B_supercell, int *displs_convolutions_B_supercell, struct convolution_char_state convolution_char_state_interim_reduced_local, struct convolution_state convolution_states_global );

#endif


// Calculate Coulomb integrals between states
#ifndef COULOMB_INTEGRALS
#define COULOMB_INTEGRALS

void coulomb_integrals( struct supercell_geometry *supercell, struct material_params *parameters, int *N_states, struct convolution_state *convolution_states, struct convolution_state *convolution_states_reversed, int local_size_conv, struct coulomb_corrections CME_corrections, struct coulomb_matrix_element CME, int type );

#endif


#ifndef WRITE_CARRIER_LOCALISATION
#define WRITE_CARRIER_LOCALISATION

void write_carrier_localisation( int N_shape, int *N_state, complex double *localisation_vb, complex double *localisation_cb );

#endif


#ifndef WRITE_ENERGIES
#define WRITE_ENERGIES

void write_energies( int *N_states, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb );

#endif


#ifndef WRITE_ENVELOPE_FUNCTIONS
#define WRITE_ENVELOPE_FUNCTIONS

void write_envelope_functions( int *N_states, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb, int direction );

#endif


#ifndef WRITE_ENVELOPE_FUNCTIONS_CROSS_SECTION
#define WRITE_ENVELOPE_FUNCTIONS_CROSS_SECTION

void write_envelope_functions_cross_section( int *N_states, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb, int direction );

#endif


#ifndef WRITE_ENVELOPE_FUNCTIONS_3D
#define WRITE_ENVELOPE_FUNCTIONS_3D

void write_envelope_functions_3D( int *N_states, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb );

#endif


#ifndef WRITE_PROBABILITY_DENSITY
#define WRITE_PROBABILITY_DENSITY

void write_probability_density( int *N_states, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb, int direction );

#endif


#ifndef WRITE_PROBABILITY_DENSITY_CROSS_SECTION
#define WRITE_PROBABILITY_DENSITY_CROSS_SECTION

void write_probability_density_cross_section( int *N_states, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb, int direction );

#endif


#ifndef WRITE_PROBABILITY_DENSITY_3D
#define WRITE_PROBABILITY_DENSITY_3D

void write_probability_density_3D( int *N_states, int N_grid, double *x, double *y, double *z, struct wave_vectors basis, struct eigenstate *eigenstate_vb, struct eigenstate *eigenstate_cb );

#endif


#ifndef WRITE_OVERLAPS
#define WRITE_OVERLAPS

void write_overlaps( int *N_states, complex double *overlaps );

#endif


#ifndef WRITE_MOMENTUM_MATRIX_ELEMENTS
#define WRITE_MOMENTUM_MATRIX_ELEMENTS

void write_momentum_matrix_elements( int Num_CB_States, int Num_VB_States, int *polarisations, complex double *momentum_matrix_elements );

#endif


#ifndef WRITE_COULOMB_MATRIX_ELEMENTS
#define WRITE_COULOMB_MATRIX_ELEMENTS

void write_coulomb_matrix_elements( int Num_VB_States, int Num_CB_States, struct coulomb_matrix_element CME );

#endif


#ifndef WRITE_CALCULATION_TIMES
#define WRITE_CALCULATION_TIMES

void write_calculation_times();

#endif


#ifndef COPY_OUTPUT
#define COPY_OUTPUT

void copy_output( int type, int direction, int N_states_vb, int N_states_cb );

#endif


#ifndef COPY_OUTPUT_LINESCAN
#define COPY_OUTPUT_LINESCAN

void copy_output_linescan( int direction, char *input_name );

#endif


#ifndef COPY_OUTPUT_CROSS_SECTION
#define COPY_OUTPUT_CROSS_SECTION

void copy_output_cross_section( int direction, char *input_name );

#endif


// Mapping local matrix row/col to global row/col
#ifndef INDXL2G
#define INDXL2G

int indxl2g(int iloc, int nb, int iproc, int isrcproc, int nprocs);

#endif


// Mapping global matrix row/col to local row/col
#ifndef INDXG2L
#define INDXG2L

int indxg2l(int iglob, int nb, int iproc, int isrcproc, int nprocs);

#endif


// Check DFTI function error status and abort if error detected
#ifndef CHECK_DFTI
#define CHECK_DFTI

void check_dfti( MKL_LONG status, const char *func );

#endif
