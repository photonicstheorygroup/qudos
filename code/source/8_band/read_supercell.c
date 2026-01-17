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

void read_supercell( int N_shapes, int strain_region_type, int coulomb_region_type, struct supercell_geometry *geometry )
{

	// WARNING! Items in supercell.in must be ordered in the exact required format, or an error may occur.
	
	// Reads supercell.in using YAML

	FILE *fp1;
	fp1 = fopen( "supercell.in", "r" );

	if ( fp1 == NULL )
	{

        printf("Error. Failed to open supercell.in file.\n");
        MPI_Abort( MPI_COMM_WORLD, 1 );

	}

	yaml_parser_t parser;
    yaml_event_t event;

    if ( !yaml_parser_initialize( &parser ) ) {
        printf("Error. Failed to initialize YAML parser.\n");
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    
	yaml_parser_set_input_file( &parser, fp1 );

	int done = 0;			// Will indicate end of input file
	int idx1 = 0;			// Will iterate over supercell shapes
	int expect_value = 0; 	// Will indicate whether parameter or value
	int len = 0;			// Will give length of 'value' pointer

	int size = 64; 			// Size of "key" holding parameter name
	char *key;
	key = malloc( sizeof(char)*size );

	char *value;
	value = malloc( sizeof(char)*size );

	// Supercell matrix default elements that cannot be altered
	geometry[0].embedded 	= 0;
	geometry[0].Ox 			= 0;
	geometry[0].Oy 			= 0;
	geometry[0].Oz 			= 0;
	geometry[0].L4 			= 0;
	geometry[0].L5 			= 0;
	geometry[0].alpha 		= 0;
	geometry[0].beta 		= 0;
	geometry[0].gamma 		= 0;

	while ( !done )
	{

		if ( !yaml_parser_parse( &parser, &event ) ) 
		{
            printf( "YAML parser error %d. Check supercell.in file.\n", parser.error );
            yaml_parser_delete( &parser );
            fclose( fp1 );
            MPI_Abort( MPI_COMM_WORLD, 1 );

        }

		switch ( event.type ) 
		{

			default:
				break;

            case YAML_STREAM_END_EVENT:  // End of file
                done = 1;
                break;

            case YAML_SCALAR_EVENT: // Detects parameter or value
				if ( idx1 == N_shapes )
				{

					printf( "Error! The parameter 'N_shape' is less than the number of shapes given in supercell.in. Remember that the supercell shape is included in N_shape.\n" );
					yaml_parser_delete( &parser );
					fclose( fp1 );
					MPI_Abort( MPI_COMM_WORLD, 1 );

				}

				len = event.data.scalar.length;  // Length of 'value'

				if ( len >= size )
				{

					len = size - 1;

				}

				memcpy( value, event.data.scalar.value, len ); // Holds parameter or value
				value[len] = '\0';

				if ( !expect_value ) // If expecting parameter, "value" holds parameter
				{
					
                    strncpy( key, value, len + 1 ); // Copies parameter to "key"
                    expect_value = 1;

				}
				else // If expecting value, "value" holds value
				{

					if ( strcmp( key, "Embedded" ) == 0 ) // Shape index the current shape is embedded within
					{

                        geometry[idx1].embedded = ( int ) strtol( value, &endptr, 10 );
                    
					}
					else if ( strcmp( key, "QD shape" ) == 0 ) // Shape type
					{

                        geometry[idx1].shape = ( int ) strtol( value, &endptr, 10 );
                    
					}
					else if ( strcmp( key, "Material" ) == 0 ) // Shape material
					{

                        geometry[idx1].material = ( int ) strtol( value, &endptr, 10 );
                    
					}
					else if ( strcmp( key, "Compos x" ) == 0 ) // x composition fraction
					{

                        geometry[idx1].x = strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Compos y" ) == 0 ) // y composition fraction
					{

                        geometry[idx1].y = strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Centre x" ) == 0 ) // Location of shape x centre
					{

                        geometry[idx1].Ox = 10*strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Centre y" ) == 0 ) // Location of shape y centre
					{

                        geometry[idx1].Oy = 10*strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Centre z" ) == 0 ) // Location of shape z centre
					{

                        geometry[idx1].Oz = 10*strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Length 1" ) == 0 ) // Shape length 1 (see characteristic_function_fourier.c)
					{

                        geometry[idx1].L1 = 10*strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Length 2" ) == 0 ) // Shape length 2 (see characteristic_function_fourier.c)
					{

                        geometry[idx1].L2 = 10*strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Length 3" ) == 0 ) // Shape length 3 (see characteristic_function_fourier.c)
					{

                        geometry[idx1].L3 = 10*strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Length 4" ) == 0 ) // Shape length 4 (see characteristic_function_fourier.c)
					{

                        geometry[idx1].L4 = 10*strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Length 5" ) == 0 ) // Shape length 5 (see characteristic_function_fourier.c)
					{

                        geometry[idx1].L5 = 10*strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Strain length 1" ) == 0 ) // Strain supercell x length
					{

                        if ( strain_region_type == 0 ) // Free lengths
						{

							geometry[idx1].SL1 = 10*strtod( value, &endptr );

						}
						else // Relative lengths
						{

							geometry[idx1].RL1 = ( int ) strtol( value, &endptr, 10 );
							geometry[idx1].SL1 = geometry[idx1].L1*geometry[idx1].RL1;

						}
                    
					}
					else if ( strcmp( key, "Strain length 2" ) == 0 ) // Strain supercell y length
					{

                        if ( strain_region_type == 0 ) // Free lengths
						{

							geometry[idx1].SL2 = 10*strtod( value, &endptr );

						}
						else // Relative lengths
						{

							geometry[idx1].RL2 = ( int ) strtol( value, &endptr, 10 );
							geometry[idx1].SL2 = geometry[idx1].L2*geometry[idx1].RL2;

						}
                    
					}
					else if ( strcmp( key, "Strain length 3" ) == 0 ) // Strain supercell z length
					{

                        if ( strain_region_type == 0 ) // Free lengths
						{

							geometry[idx1].SL3 = 10*strtod( value, &endptr );

						}
						else // Relative lengths
						{

							geometry[idx1].RL3 = ( int ) strtol( value, &endptr, 10 );
							geometry[idx1].SL3 = geometry[idx1].L3*geometry[idx1].RL3;

						}
                    
					}
					else if ( strcmp( key, "Coulomb length 1" ) == 0 ) // Coulomb supercell x length
					{

                        if ( coulomb_region_type == 0 ) // Free lengths
						{

							geometry[idx1].CL1 = 10*strtod( value, &endptr );

						}
						else // Relative lengths
						{

							geometry[idx1].RCL1 = ( int ) strtol( value, &endptr, 10 );
							geometry[idx1].CL1 = geometry[idx1].L1*geometry[idx1].RCL1;

						}
                    
					}
					else if ( strcmp( key, "Coulomb length 2" ) == 0 ) // Coulomb supercell y length
					{

                        if ( coulomb_region_type == 0 ) // Free lengths
						{

							geometry[idx1].CL2 = 10*strtod( value, &endptr );

						}
						else // Relative lengths
						{

							geometry[idx1].RCL2 = ( int ) strtol( value, &endptr, 10 );
							geometry[idx1].CL2 = geometry[idx1].L2*geometry[idx1].RCL2;

						}
                    
					}
					else if ( strcmp( key, "Coulomb length 3" ) == 0 ) // Coulomb supercell z length
					{

                        if ( coulomb_region_type == 0 ) // Free lengths
						{

							geometry[idx1].CL3 = 10*strtod( value, &endptr );

						}
						else // Relative lengths
						{

							geometry[idx1].RCL3 = ( int ) strtol( value, &endptr, 10 );
							geometry[idx1].CL3 = geometry[idx1].L3*geometry[idx1].RCL3;

						}
						idx1++; // Increments shape for current shape index 0
                    
					}
					else if ( strcmp( key, "Euler alpha angle" ) == 0 ) // Shape alpha rotation angle
					{

                        geometry[idx1].alpha = strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Euler beta angle" ) == 0 ) // Shape beta rotation angle
					{

                        geometry[idx1].beta = strtod( value, &endptr );
                    
					}
					else if ( strcmp( key, "Euler gamma angle" ) == 0 ) // Shape gamma rotation angle
					{

                        geometry[idx1].gamma = strtod( value, &endptr );
						idx1++; // Increments shape for current shape index > 0
                    
					}
					
					expect_value = 0;

				}

				break;


		}

		yaml_event_delete(&event);


	}

	if ( idx1 < N_shapes )
	{

		printf( "Error! The parameter 'N_shape' is greater than the number of shapes given in supercell.in. Remember that the supercell shape is included in N_shape.\n" );
		yaml_parser_delete( &parser );
		fclose( fp1 );
		MPI_Abort( MPI_COMM_WORLD, 1 );

	}

	yaml_parser_delete(&parser);

	fclose( fp1 );

	free( key );
	free( value );

}
