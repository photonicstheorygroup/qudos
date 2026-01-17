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

void read_parameters( int *N_shape, int *N_Gx, int *N_Gy, int *N_Gz, int *strain_region_type, int *N_S_Gx, int *N_S_Gy, int *N_S_Gz, double *T, int *piezo_order, double *Fx, double *Fy, double *Fz, int *N_grid, int *isotropic, int *include_shape_elastic_constants, int *shape_elastic_constants_order, int *softening, double *sigma_softening, int *directions, int* cross_sections, int *N_grid_cross_sections, int *directions_cross_sections, int *input_check )
{

	// Reads parameters_linescan.in using YAML

	FILE *fp1;
	fp1 = fopen( "parameters_linescan.in", "r" );

	if ( fp1 == NULL )
	{

        printf("Error. Failed to open parameters_linescan.in file.\n");
        exit( 1 );

	}

	yaml_parser_t parser;
    yaml_event_t event;

    if ( !yaml_parser_initialize( &parser ) ) {
        printf("Error. Failed to initialize YAML parser.\n");
        exit( 1 );
    }
    
	yaml_parser_set_input_file( &parser, fp1 );

	int done = 0;			// Will indicate end of input file
	int expect_value = 0; 	// Will indicate whether parameter or value
	int len = 0;			// Will give length of 'value' pointer

	int size = 64; 			// Size of "key" holding parameter name
	char *key;
	key = malloc( sizeof(char)*size );

	char *value;
	value = malloc( sizeof(char)*size );

	while ( !done )
	{

		if ( !yaml_parser_parse( &parser, &event ) ) 
		{
            printf( "YAML parser error %d. Check parameters_linescan.in file.\n", parser.error );
            yaml_parser_delete( &parser );
            fclose( fp1 );
            exit( 1 );

        }

		switch ( event.type ) 
		{

			default:
				break;

            case YAML_STREAM_END_EVENT: // End of file
                done = 1;
                break;

            case YAML_SCALAR_EVENT:	// Detects parameter or value
				len = event.data.scalar.length;  // Length of 'value'

				if ( len >= size )
				{

					len = size - 1;

				}

				memcpy( value, event.data.scalar.value, len ); // Holds parameter or value
				value[len] = '\0';

				if ( !expect_value )  // If expecting parameter, "value" holds parameter
				{
					
                    strncpy( key, value, len + 1 ); // Copies parameter to "key"
                    expect_value = 1;

				}
				else // If expecting value, "value" holds value
				{

					if ( strcmp( key, "N_shape" ) == 0 ) // Number of shapes in supercell (including supercell)
					{

						*N_shape = ( int ) strtol( value, &endptr, 10 );
						input_check[0] = 1;
					
					}
					else if ( strcmp( key, "M_Gx" ) == 0 ) // Number of plane waves for electronic structure supercell
					{

						*N_Gx = 2*( ( int ) strtol( value, &endptr, 10 ) ) + 1; // N_Gx = 2*M_Gx + 1 = Total no. of plane waves along x
						input_check[1] = 1;
					
					}
					else if ( strcmp( key, "M_Gy" ) == 0 ) // Number of plane waves for electronic structure supercell
					{

						*N_Gy = 2*( ( int ) strtol( value, &endptr, 10 ) ) + 1; // N_Gy = 2*M_Gy + 1 = Total no. of plane waves along y
						input_check[2] = 1;
					
					}
					else if ( strcmp( key, "M_Gz" ) == 0 ) // Number of plane waves for electronic structure supercell
					{

						*N_Gz = 2*( ( int ) strtol( value, &endptr, 10 ) ) + 1; // N_Gz = 2*M_Gz + 1 = Total no. of plane waves along z
						input_check[3] = 1;
					
					}
					else if ( strcmp( key, "Strain_supercell_type" ) == 0 ) // Indicates whether to use free strain supercell lengths (0) or to set them to be integer multiples of the electronic structure supercell lengths (1) allowing much faster strain and piezoelectric potential conversion from the strain supercell to the electronic structure supercell.
					{

						*strain_region_type = ( int ) strtol( value, &endptr, 10 );
						input_check[4] = 1;
					
					}
					else if ( strcmp( key, "M_S_Gx" ) == 0 ) // Plane wave basis set size for the expanded strain supercell.
					{

						*N_S_Gx = 2*( ( int ) strtol( value, &endptr, 10 ) ) + 1; // N_S_Gx = 2*M_S_Gx + 1 = Total no. of plane waves along x
						input_check[5] = 1;
					
					}
					else if ( strcmp( key, "M_S_Gy" ) == 0 ) // Plane wave basis set size for the expanded strain supercell.
					{

						*N_S_Gy = 2*( ( int ) strtol( value, &endptr, 10 ) ) + 1; // N_S_Gy = 2*M_S_Gy + 1 = Total no. of plane waves along y
						input_check[6] = 1;
					
					}
					else if ( strcmp( key, "M_S_Gz" ) == 0 ) // Plane wave basis set size for the expanded strain supercell.
					{

						*N_S_Gz = 2*( ( int ) strtol( value, &endptr, 10 ) ) + 1; // N_S_Gz = 2*M_S_Gz + 1 = Total no. of plane waves along z
						input_check[7] = 1;
					
					}
					else if ( strcmp( key, "T" ) == 0 ) // System temperature [K]
					{

						*T = strtod( value, &endptr );
						input_check[8] = 1;
					
					}
					else if ( strcmp( key, "Fx" ) == 0 ) // Applied electric field x component - read in kV/cm and convert to V/Angstrom
					{

						*Fx = 1.0e-05*strtod( value, &endptr ); // Read in kV/cm and convert to V/Angstrom
						input_check[9] = 1;
					
					}
					else if ( strcmp( key, "Fy" ) == 0 ) // Applied electric field y component - read in kV/cm and convert to V/Angstrom
					{

						*Fy = 1.0e-05*strtod( value, &endptr ); // Read in kV/cm and convert to V/Angstrom
						input_check[10] = 1;
					
					}
					else if ( strcmp( key, "Fz" ) == 0 ) // Applied electric field z component - read in kV/cm and convert to V/Angstrom
					{

						*Fz = 1.0e-05*strtod( value, &endptr ); // Read in kV/cm and convert to V/Angstrom
						input_check[11] = 1;
					
					}
					else if ( strcmp( key, "N_line" ) == 0 ) // Number of real space grid points for linescans
					{

						*N_grid = ( int ) strtol( value, &endptr, 10 );
						input_check[12] = 1;
					
					}
					else if ( strcmp( key, "Isotropic" ) == 0 ) // Isotropic approximation for elastic tensor
					{

						*isotropic = ( int ) strtol( value, &endptr, 10 ); // = 0 for full anisotropic elastic tensor, = 1 for isotropic approximation to elastic tensor
						input_check[13] = 1;
					
					}
					else if ( strcmp( key, "Elastic_const" ) == 0 ) //  Material dependent elastic constants
					{

						*include_shape_elastic_constants = ( int ) strtol( value, &endptr, 10 ); // = 0 for matrix material elastic constants, = 1 for shape dependent elastic constants
						input_check[14] = 1;
					
					}
					else if ( strcmp( key, "Elastic_const_order" ) == 0 ) // Correction order for position-dependent elastic constants
					{

						*shape_elastic_constants_order = ( int ) strtol( value, &endptr, 10 ); // order of correction to strain tensor for shape elastic constants
						input_check[15] = 1;
					
					}
					else if ( strcmp( key, "Piezo" ) == 0 ) // Include piezoelectric potential in calculation
					{

						*piezo_order = ( int ) strtol( value, &endptr, 10 ); // = 0 for no piezoelectricity, = 1 for 1st order only, = 2 for 1st and 2nd order, = 3 for 2nd order only
						input_check[16] = 1;
					
					}
					else if ( strcmp( key, "Softening" ) == 0 ) // Interface softening via Gaussian convolution
					{

						*softening = ( int ) strtol( value, &endptr, 10 ); // = 1 for interface softening, = 0 for "hard" interface calculation
						input_check[17] = 1;
					
					}
					else if ( strcmp( key, "Sigma_softening" ) == 0 ) // Include state overlaps calculation
					{

						*sigma_softening = 10.0*strtod( value, &endptr     ); // Standard deviation of Gaussian used for interface softening, read in nm and convert to Angstrom (default: 1.0 Angstrom = 0.1 nm)
						input_check[18] = 1;
					
					}
					else if ( strcmp( key, "Linescan_100" ) == 0 ) // Set directions along which to compute real space linescans (linescans go through centre of supercell)
					{

						directions[0] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute linescan along [100]
						input_check[19] = 1;
					
					}
					else if ( strcmp( key, "Linescan_010" ) == 0 ) // Set directions along which to compute real space linescans (linescans go through centre of supercell)
					{

						directions[1] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute linescan along [010]
						input_check[20] = 1;
					
					}
					else if ( strcmp( key, "Linescan_001" ) == 0 ) // Set directions along which to compute real space linescans (linescans go through centre of supercell)
					{

						directions[2] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute linescan along [001]
						input_check[21] = 1;
					
					}
					else if ( strcmp( key, "Linescan_110" ) == 0 ) // Set directions along which to compute real space linescans (linescans go through centre of supercell)
					{

						directions[3] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute linescan along [110]
						input_check[22] = 1;
					
					}
					else if ( strcmp( key, "Linescan_101" ) == 0 ) // Set directions along which to compute real space linescans (linescans go through centre of supercell)
					{

						directions[4] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute linescan along [101]
						input_check[23] = 1;
					
					}
					else if ( strcmp( key, "Linescan_011" ) == 0 ) // Set directions along which to compute real space linescans (linescans go through centre of supercell)
					{

						directions[5] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute linescan along [011]
						input_check[24] = 1;
					
					}
					else if ( strcmp( key, "Linescan_111" ) == 0 ) // Set directions along which to compute real space linescans (linescans go through centre of supercell)
					{

						directions[6] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute linescan along [111]
						input_check[25] = 1;
					
					}
					else if ( strcmp( key, "Cross_sections" ) == 0 ) // Indicate whether to output cross-sections of envelope functions
					{

						*cross_sections = ( int ) strtol( value, &endptr, 10 ); // = 1 for cross sections, = 0 for no cross sections
						input_check[26] = 1;
					
					}
					else if ( strcmp( key, "N_cross" ) == 0 ) // Number of real-space grid points for each of the two dimensions for cross-sections
					{

						*N_grid_cross_sections = ( int ) strtol( value, &endptr, 10 );
						input_check[27] = 1;
					
					}
					else if ( strcmp( key, "Cross_section_100" ) == 0 ) // Set normals of which to compute real space cross-sections (cross-sections go through centre of supercell)
					{

						directions_cross_sections[0] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute cross-section with normal along [100]
						input_check[28] = 1;
					
					}
					else if ( strcmp( key, "Cross_section_010" ) == 0 ) // Set normals of which to compute real space cross-sections (cross-sections go through centre of supercell)
					{

						directions_cross_sections[1] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute cross-section with normal along [010]
						input_check[29] = 1;
					
					}
					else if ( strcmp( key, "Cross_section_001" ) == 0 ) // Set normals of which to compute real space cross-sections (cross-sections go through centre of supercell)
					{

						directions_cross_sections[2] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute cross-section with normal along [001]
						input_check[30] = 1;
					
					}
					else if ( strcmp( key, "Cross_section_m110" ) == 0 ) // Set normals of which to compute real space cross-sections (cross-sections go through centre of supercell)
					{

						directions_cross_sections[3] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute cross-section with normal along [m110]
						input_check[31] = 1;
					
					}
					else if ( strcmp( key, "Cross_section_110" ) == 0 ) // Set normals of which to compute real space cross-sections (cross-sections go through centre of supercell)
					{

						directions_cross_sections[4] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute cross-section with normal along [110]
						input_check[32] = 1;
					
					}
					else if ( strcmp( key, "Cross_section_m101" ) == 0 ) // Set normals of which to compute real space cross-sections (cross-sections go through centre of supercell)
					{

						directions_cross_sections[5] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute cross-section with normal along [m101]
						input_check[33] = 1;
					
					}
					else if ( strcmp( key, "Cross_section_0m11" ) == 0 ) // Set normals of which to compute real space cross-sections (cross-sections go through centre of supercell)
					{

						directions_cross_sections[6] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute cross-section with normal along [0m11]
						input_check[34] = 1;
					
					}
					else if ( strcmp( key, "Cross_section_111" ) == 0 ) // Set normals of which to compute real space cross-sections (cross-sections go through centre of supercell)
					{

						directions_cross_sections[7] = ( int ) strtol( value, &endptr, 10 ); // = 1 to compute cross-section with normal along [11]
						input_check[35] = 1;
					
					}
					
					expect_value = 0;

				}

				break;


		}

		yaml_event_delete(&event);


	}

	yaml_parser_delete(&parser);

	fclose( fp1 );

	free( key );
	free( value );

}
