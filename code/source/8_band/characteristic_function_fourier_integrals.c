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


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//// I_e0 integral for the pyramid and truncated pyramid characteristic function Fourier transforms ////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

complex double I_e0( double a, double G )
{

	complex double result;

	if( fabs( G ) < EPS )
	{

		result = a;

	}
	else
	{

		result = ( cexp( I*G*a ) - 1.0 )/( I*G );

	}

	return result;

}


//////////////////////////////////////////////////////////////////////////////////////////////
//// I_J0 integral for the hemisphere and lens characteristic function Fourier transforms ////
//////////////////////////////////////////////////////////////////////////////////////////////

double I_J0( double alpha, double beta ) // Using trapezoidal rule (for reference)
{

	int idx_x, N_x = 200;
	double *x, dx, *integrand, result = 0.0;


	// Set x grid

	x = malloc( N_x*sizeof( double ) );

	x[0      ] = 0.0; // Lower limit on integral
	x[N_x - 1] = 1.0; // Upper limit on integral

	dx = ( x[N_x - 1] - x[0] )/( ( double ) ( N_x - 1 ) );

	for( idx_x = 1; idx_x < N_x - 1; idx_x++ )
	{

		x[idx_x] = x[idx_x - 1] + dx;

	}


	// Set integrand

	integrand = malloc( N_x*sizeof( double ) );

	for( idx_x = 0; idx_x < N_x; idx_x++ )
	{

		integrand[idx_x] = x[idx_x]*cos( alpha*sqrt( 1.0 - x[idx_x]*x[idx_x] ) )*gsl_sf_bessel_J0( beta*x[idx_x] );

	}


	// Perform integral using trapezoidal rule

	result = trapezoidal_rule( N_x, x, integrand );


	// Free dynamically allocated memory

	free( x         );
	free( integrand );


	return result;

}


double I_J0_GSL( double alpha, double beta ) // Using GSL
{

	double result, error, lim_lo, lim_hi;

	lim_lo = 0.0; // Upper limit on integral
	lim_hi = 1.0; // Lower limit on integral

	int key = 6;


	// Use the GSL routine QAGS to evaluate the integral for the supplied alpha, beta and limits

	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( LIMIT );

	struct parameters integrand_params;
	integrand_params.alpha = alpha;
	integrand_params. beta =  beta;

	gsl_function integrand;
	integrand.function = &integrand_I_J0;
	integrand.params   = &integrand_params;
	//gsl_integration_qags( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, workspace, &result, &error );
	gsl_integration_qag( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, key, workspace, &result, &error );


	gsl_integration_workspace_free( workspace );

	return result;

}


double integrand_I_J0( double x, void *params )
{

	double result;
	struct parameters integrand_params = *( struct parameters * ) params;

	result = x*cos( integrand_params.alpha*sqrt( 1.0 - x*x ) )*gsl_sf_bessel_J0( integrand_params.beta*x );

	return result;

}



//////////////////////////////////////////////////////////////////////////////////////////////////
//// I_2 integral for the cone and elliptical cone characteristic function Fourier transforms ////
//////////////////////////////////////////////////////////////////////////////////////////////////

complex double I_2( double alpha, double beta ) // Using trapezoidal rule (for reference)
{

	int idx_x, N_x = 100;
	double *x, dx, *integrand_real, *integrand_imag, result_real, result_imag;
	complex double result;


	// Set x grid

	x = malloc( N_x*sizeof( double ) );

	x[0      ] = 0.0; // Lower limit on integral
	x[N_x - 1] = 1.0; // Upper limit on integral

	dx = ( x[N_x - 1] - x[0] )/( ( double ) ( N_x - 1 ) );

	for( idx_x = 1; idx_x < N_x - 1; idx_x++ )
	{

		x[idx_x] = x[idx_x - 1] + dx;

	}


	// Set integrands

	integrand_real = malloc( N_x*sizeof( double ) );
	integrand_imag = malloc( N_x*sizeof( double ) );

	for( idx_x = 0; idx_x < N_x; idx_x++ )
	{

		integrand_real[idx_x] = x[idx_x]*gsl_sf_bessel_J0( alpha*x[idx_x] )*cos( beta*x[idx_x] );
		integrand_imag[idx_x] = x[idx_x]*gsl_sf_bessel_J0( alpha*x[idx_x] )*sin( beta*x[idx_x] );

	}


	// Perform integral using trapezoidal rule

	result_real = trapezoidal_rule( N_x, x, integrand_real );
	result_imag = trapezoidal_rule( N_x, x, integrand_imag );

	result = result_real + I*result_imag;


	// Free dynamically allocated memory

	free( x              );
	free( integrand_real );
	free( integrand_imag );


	return result;

}


complex double I_2_GSL(double alpha, double beta ) // Using GSL
{

	double result_real, result_imag, error, lim_lo, lim_hi;
	complex double result;

	lim_lo = 0.0; // Upper limit on integral
	lim_hi = 1.0; // Lower limit on integral

	int key = 6;


	// Use the GSL routine QAGS to evaluate the integral for the supplied alpha, beta and limits

	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( LIMIT );

	struct parameters integrand_params;
	integrand_params.alpha = alpha;
	integrand_params. beta =  beta;

	gsl_function integrand;
	integrand.function = &integrand_I_2_real;
	integrand.params   = &integrand_params;
	//gsl_integration_qags( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, workspace, &result, &error );
	gsl_integration_qag( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, key, workspace, &result_real, &error );

	integrand.function = &integrand_I_2_imag;
	integrand.params   = &integrand_params;
	//gsl_integration_qags( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, workspace, &result, &error );
	gsl_integration_qag( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, key, workspace, &result_imag, &error );

	gsl_integration_workspace_free( workspace );

	result = result_real + I*result_imag;

	return result;

}


double integrand_I_2_real( double x, void *params )
{

	double result;
	struct parameters integrand_params = *( struct parameters * ) params;

	result = x*gsl_sf_bessel_J0( integrand_params.alpha*x )*cos( integrand_params.beta*x );

	return result;

}


double integrand_I_2_imag( double x, void *params )
{

	double result;
	struct parameters integrand_params = *( struct parameters * ) params;

	result = x*gsl_sf_bessel_J0( integrand_params.alpha*x )*sin( integrand_params.beta*x );

	return result;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// I_2_Gz integral for the cone and elliptical cone characteristic function Fourier transforms when Gz = 0 ////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

complex double I_2_Gz( double alpha) // Using trapezoidal rule (for reference)
{

	int idx_x, N_x = 100;
	double *x, dx, *integrand_real, result_real;
	complex double result;


	// Set x grid

	x = malloc( N_x*sizeof( double ) );

	x[0      ] = 0.0; // Lower limit on integral
	x[N_x - 1] = 1.0; // Upper limit on integral

	dx = ( x[N_x - 1] - x[0] )/( ( double ) ( N_x - 1 ) );

	for( idx_x = 1; idx_x < N_x - 1; idx_x++ )
	{

		x[idx_x] = x[idx_x - 1] + dx;

	}


	// Set integrands

	integrand_real = malloc( N_x*sizeof( double ) );

	for( idx_x = 0; idx_x < N_x; idx_x++ )
	{

		integrand_real[idx_x] = x[idx_x]*x[idx_x]*gsl_sf_bessel_J0( alpha*x[idx_x] );

	}


	// Perform integral using trapezoidal rule

	result_real = trapezoidal_rule( N_x, x, integrand_real );

	result = result_real;


	// Free dynamically allocated memory

	free( x              );
	free( integrand_real );


	return result;

}


complex double I_2_Gz_GSL( double alpha ) // Using GSL
{

	double result, result_real, error, lim_lo, lim_hi;

	lim_lo = 0.0; // Upper limit on integral
	lim_hi = 1.0; // Lower limit on integral

	int key = 6;


	// Use the GSL routine QAGS to evaluate the integral for the supplied alpha, beta and limits

	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( LIMIT );

	struct parameters integrand_params;
	integrand_params.alpha = alpha;

	gsl_function integrand;
	integrand.function = &integrand_I_2_Gz_real;
	integrand.params   = &integrand_params;
	//gsl_integration_qags( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, workspace, &result, &error );
	gsl_integration_qag( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, key, workspace, &result_real, &error );

	gsl_integration_workspace_free( workspace );

	result = result_real;

	return result;

}


double integrand_I_2_Gz_real( double x, void *params )
{

	double result;
	struct parameters integrand_params = *( struct parameters * ) params;

	result = x*x*gsl_sf_bessel_J0( integrand_params.alpha*x );

	return result;

}


//////////////////////////////////////////////////////////////////////////////////////////
//// I_2_TC integral for the truncated cone characteristic function Fourier transform ////
//////////////////////////////////////////////////////////////////////////////////////////

complex double I_2_TC(double alpha, double R1, double R2, double beta ) // Using trapezoidal rule (for reference)
{

	int idx_x, N_x = 100;
	double *x, dx, *integrand_real, *integrand_imag, result_real, result_imag;
	complex double result;


	// Set x grid

	x = malloc( N_x*sizeof( double ) );

	x[0      ] = R1 / (R2 - R1); // Lower limit on integral
	x[N_x - 1] = R2 / (R2 - R1); // Upper limit on integral

	dx = ( x[N_x - 1] - x[0] )/( ( double ) ( N_x - 1 ) );

	for( idx_x = 1; idx_x < N_x - 1; idx_x++ )
	{

		x[idx_x] = x[idx_x - 1] + dx;

	}


	// Set integrands

	integrand_real = malloc( N_x*sizeof( double ) );
	integrand_imag = malloc( N_x*sizeof( double ) );

	for( idx_x = 0; idx_x < N_x; idx_x++ )
	{

		integrand_real[idx_x] = x[idx_x]*gsl_sf_bessel_J0( (R2 - R1)*alpha*x[idx_x] )*cos( beta*x[idx_x] );
		integrand_imag[idx_x] = x[idx_x]*gsl_sf_bessel_J0( (R2 - R1)*alpha*x[idx_x] )*sin( beta*x[idx_x] );

	}


	// Perform integral using trapezoidal rule

	result_real = trapezoidal_rule( N_x, x, integrand_real );
	result_imag = trapezoidal_rule( N_x, x, integrand_imag );

	result = result_real + I*result_imag;


	// Free dynamically allocated memory

	free( x              );
	free( integrand_real );
	free( integrand_imag );


	return result;

}


complex double I_2_TC_GSL(double alpha, double R1, double R2, double beta ) // Using GSL
{

	double result_real, result_imag, error, lim_lo, lim_hi;
	complex double result;

	lim_lo = R1 / (R2 - R1); // Upper limit on integral
	lim_hi = R2 / (R2 - R1); // Lower limit on integral

	int key = 6;


	// Use the GSL routine QAGS to evaluate the integral for the supplied alpha, beta and limits

	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( LIMIT );

	struct parameters integrand_params;
	integrand_params.alpha = alpha;
	integrand_params. beta =  beta;
	integrand_params.R1 = R1;
	integrand_params.R2 = R2;

	gsl_function integrand;
	integrand.function = &integrand_I_2_TC_real;
	integrand.params   = &integrand_params;
	//gsl_integration_qags( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, workspace, &result, &error );
	gsl_integration_qag( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, key, workspace, &result_real, &error );

	integrand.function = &integrand_I_2_TC_imag;
	integrand.params   = &integrand_params;
	//gsl_integration_qags( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, workspace, &result, &error );
	gsl_integration_qag( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, key, workspace, &result_imag, &error );

	gsl_integration_workspace_free( workspace );

	result = result_real + I*result_imag;

	return result;

}


double integrand_I_2_TC_real( double x, void *params )
{

	double result;
	struct parameters integrand_params = *( struct parameters * ) params;

	result = x*gsl_sf_bessel_J0( (integrand_params.R2 - integrand_params.R1)*integrand_params.alpha*x )*cos( integrand_params.beta*x );

	return result;

}


double integrand_I_2_TC_imag( double x, void *params )
{

	double result;
	struct parameters integrand_params = *( struct parameters * ) params;

	result = x*gsl_sf_bessel_J0( (integrand_params.R2 - integrand_params.R1)*integrand_params.alpha*x )*sin( integrand_params.beta*x );

	return result;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//// I_2_TC_Gz integral for the truncated cone characteristic function Fourier transform when Gz = 0////
////////////////////////////////////////////////////////////////////////////////////////////////////////

complex double I_2_TC_Gz(double alpha, double R1, double R2) // Using trapezoidal rule (for reference)
{

	int idx_x, N_x = 100;
	double *x, dx, *integrand_real, result_real;
	complex double result;


	// Set x grid

	x = malloc( N_x*sizeof( double ) );

	x[0      ] = R1 / (R2 - R1); // Lower limit on integral
	x[N_x - 1] = R2 / (R2 - R1); // Upper limit on integral

	dx = ( x[N_x - 1] - x[0] )/( ( double ) ( N_x - 1 ) );

	for( idx_x = 1; idx_x < N_x - 1; idx_x++ )
	{

		x[idx_x] = x[idx_x - 1] + dx;

	}


	// Set integrands

	integrand_real = malloc( N_x*sizeof( double ) );

	for( idx_x = 0; idx_x < N_x; idx_x++ )
	{

		integrand_real[idx_x] = x[idx_x]*x[idx_x]*gsl_sf_bessel_J0( (R2 - R1)*alpha*x[idx_x] );

	}


	// Perform integral using trapezoidal rule

	result_real = trapezoidal_rule( N_x, x, integrand_real );

	result = result_real;

	// Free dynamically allocated memory

	free( x              );
	free( integrand_real );


	return result;

}


complex double I_2_TC_Gz_GSL(double alpha, double R1, double R2) // Using GSL
{

	double result, result_real, error, lim_lo, lim_hi;

	lim_lo = R1 / (R2 - R1); // Upper limit on integral
	lim_hi = R2 / (R2 - R1); // Lower limit on integral

	int key = 6;


	// Use the GSL routine QAGS to evaluate the integral for the supplied alpha, beta and limits

	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( LIMIT );

	struct parameters integrand_params;
	integrand_params.alpha = alpha;
	integrand_params.R1 = R1;
	integrand_params.R2 = R2;

	gsl_function integrand;
	integrand.function = &integrand_I_2_TC_Gz_real;
	integrand.params   = &integrand_params;
	//gsl_integration_qags( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, workspace, &result, &error );
	gsl_integration_qag( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, key, workspace, &result_real, &error );

	gsl_integration_workspace_free( workspace );

	result = result_real;

	return result;

}


double integrand_I_2_TC_Gz_real( double x, void *params )
{

	double result;
	struct parameters integrand_params = *( struct parameters * ) params;

	result = x*x*gsl_sf_bessel_J0( (integrand_params.R2 - integrand_params.R1)*integrand_params.alpha*x );

	return result;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//// I_2_ETC integral for the elliptical truncated cone characteristic function Fourier transform ////
//////////////////////////////////////////////////////////////////////////////////////////////////////

complex double I_2_ETC(double alpha, double beta, double f) // Using trapezoidal rule (for reference)
{

	int idx_x, N_x = 100;
	double *x, dx, *integrand_real, *integrand_imag, result_real, result_imag;
	complex double result;


	// Set x grid

	x = malloc( N_x*sizeof( double ) );

	x[0      ] = f; // Lower limit on integral
	x[N_x - 1] = 1.0; // Upper limit on integral

	dx = ( x[N_x - 1] - x[0] )/( ( double ) ( N_x - 1 ) );

	for( idx_x = 1; idx_x < N_x - 1; idx_x++ )
	{

		x[idx_x] = x[idx_x - 1] + dx;

	}


	// Set integrands

	integrand_real = malloc( N_x*sizeof( double ) );
	integrand_imag = malloc( N_x*sizeof( double ) );

	for( idx_x = 0; idx_x < N_x; idx_x++ )
	{

		integrand_real[idx_x] = x[idx_x]*gsl_sf_bessel_J0( alpha*x[idx_x] )*cos( beta*x[idx_x]/(1.0-f) );
		integrand_imag[idx_x] = x[idx_x]*gsl_sf_bessel_J0( alpha*x[idx_x] )*sin( beta*x[idx_x]/(1.0-f) );

	}


	// Perform integral using trapezoidal rule

	result_real = trapezoidal_rule( N_x, x, integrand_real );
	result_imag = trapezoidal_rule( N_x, x, integrand_imag );

	result = result_real + I*result_imag;


	// Free dynamically allocated memory

	free( x              );
	free( integrand_real );
	free( integrand_imag );


	return result;

}


complex double I_2_ETC_GSL(double alpha, double beta, double f ) // Using GSL
{

	double result_real, result_imag, error, lim_lo, lim_hi;
	complex double result;

	lim_lo = f; // Upper limit on integral
	lim_hi = 1.0; // Lower limit on integral

	int key = 6;


	// Use the GSL routine QAGS to evaluate the integral for the supplied alpha, beta and limits

	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( LIMIT );

	struct parameters integrand_params;
	integrand_params.alpha = alpha;
	integrand_params. beta =  beta;
	integrand_params.f =  f;

	gsl_function integrand;
	integrand.function = &integrand_I_2_ETC_real;
	integrand.params   = &integrand_params;
	//gsl_integration_qags( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, workspace, &result, &error );
	gsl_integration_qag( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, key, workspace, &result_real, &error );

	integrand.function = &integrand_I_2_ETC_imag;
	integrand.params   = &integrand_params;
	//gsl_integration_qags( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, workspace, &result, &error );
	gsl_integration_qag( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, key, workspace, &result_imag, &error );

	gsl_integration_workspace_free( workspace );

	result = result_real + I*result_imag;

	return result;

}


double integrand_I_2_ETC_real( double x, void *params )
{

	double result;
	struct parameters integrand_params = *( struct parameters * ) params;

	result = x*gsl_sf_bessel_J0( integrand_params.alpha*x )*cos( integrand_params.beta*x/(1.0-integrand_params.f) );

	return result;

}


double integrand_I_2_ETC_imag( double x, void *params )
{

	double result;
	struct parameters integrand_params = *( struct parameters * ) params;

	result = x*gsl_sf_bessel_J0( integrand_params.alpha*x )*sin( integrand_params.beta*x/(1.0-integrand_params.f) );

	return result;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//// I_2_ETC_Gz integral for the truncated cone characteristic function Fourier transform when Gz = 0 ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

complex double I_2_ETC_Gz(double alpha, double f) // Using trapezoidal rule (for reference)
{

	int idx_x, N_x = 100;
	double *x, dx, *integrand_real, result_real;
	complex double result;


	// Set x grid

	x = malloc( N_x*sizeof( double ) );

	x[0      ] = f; // Lower limit on integral
	x[N_x - 1] = 1.0; // Upper limit on integral

	dx = ( x[N_x - 1] - x[0] )/( ( double ) ( N_x - 1 ) );

	for( idx_x = 1; idx_x < N_x - 1; idx_x++ )
	{

		x[idx_x] = x[idx_x - 1] + dx;

	}


	// Set integrands

	integrand_real = malloc( N_x*sizeof( double ) );

	for( idx_x = 0; idx_x < N_x; idx_x++ )
	{

		integrand_real[idx_x] = x[idx_x]*x[idx_x]*gsl_sf_bessel_J0( alpha*x[idx_x] );

	}


	// Perform integral using trapezoidal rule

	result_real = trapezoidal_rule( N_x, x, integrand_real );

	result = result_real;

	// Free dynamically allocated memory

	free( x              );
	free( integrand_real );


	return result;

}


complex double I_2_ETC_Gz_GSL( double alpha, double f ) // Using GSL
{

	double result, result_real, error, lim_lo, lim_hi;

	lim_lo = f; // Upper limit on integral
	lim_hi = 1.0; // Lower limit on integral

	int key = 6;


	// Use the GSL routine QAGS to evaluate the integral for the supplied alpha, beta and limits

	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( LIMIT );

	struct parameters integrand_params;
	integrand_params.alpha = alpha;
	integrand_params.f = f;

	gsl_function integrand;
	integrand.function = &integrand_I_2_ETC_Gz_real;
	integrand.params   = &integrand_params;
	//gsl_integration_qags( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, workspace, &result, &error );
	gsl_integration_qag( &integrand, lim_lo, lim_hi, EPS_ABS, EPS_REL, LIMIT, key, workspace, &result_real, &error );

	gsl_integration_workspace_free( workspace );

	result = result_real;

	return result;

}


double integrand_I_2_ETC_Gz_real( double x, void *params )
{

	double result;
	struct parameters integrand_params = *( struct parameters * ) params;

	result = x*x*gsl_sf_bessel_J0( integrand_params.alpha*x );

	return result;

}
