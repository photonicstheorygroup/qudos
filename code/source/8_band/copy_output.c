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

// Move linescans, cross-sections, and 3D output to output files folder
void copy_output( int type, int direction, int N_states_vb, int N_states_cb )
{

	char create_graph_folder[80] = "mkdir -p ";
	strcat( create_graph_folder, graph_folder );
	system( create_graph_folder ); // Creates graphs folder if it does not exist
	
	if ( type == 0 ) // Linescans
	{
		
		int size = 160;
		
		char *input_name;
		input_name = malloc( sizeof(char)*size );


		// Valence band probability densities & envelope functions
		for (int i = 0; i < N_states_vb; i++)
		{

			char str1[size];
			sprintf( str1, "%d", i + 1 );

			strcpy( input_name, "probability_density_vb_" );
			strcat( input_name, str1 );
			strcat( input_name, "_" );
			copy_output_linescan( direction, input_name );

			strcpy( input_name, "envelope_function_vb_" );
			strcat( input_name, str1 );
			strcat( input_name, "_" );
			copy_output_linescan( direction, input_name );

		}

		// Conduction band probability densities & envelope functions
		for (int i = 0; i < N_states_cb; i++)
		{

			char str1[size];
			sprintf( str1, "%d", i + 1 );

			strcpy( input_name, "probability_density_cb_" );
			strcat( input_name, str1 );
			strcat( input_name, "_" );
			copy_output_linescan( direction, input_name );

			strcpy( input_name, "envelope_function_cb_" );
			strcat( input_name, str1 );
			strcat( input_name, "_" );
			copy_output_linescan( direction, input_name );

		}
		
		free( input_name );
		
	}
	else if ( type == 1 ) // Cross-sections
	{
		
		int size = 160;
		
		char *input_name;
		input_name = malloc( sizeof(char)*size );


		// Valence band probability densities & envelope functions
		for (int i = 0; i < N_states_vb; i++)
		{

			char str1[size];
			sprintf( str1, "%d", i + 1 );

			strcpy( input_name, "probability_density_cross_section_vb_" );
			strcat( input_name, str1 );
			strcat( input_name, "_" );
			copy_output_cross_section( direction, input_name );

			strcpy( input_name, "envelope_function_cross_section_vb_" );
			strcat( input_name, str1 );
			strcat( input_name, "_" );
			copy_output_cross_section( direction, input_name );

		}

		// Conduction band probability densities & envelope functions
		for (int i = 0; i < N_states_cb; i++)
		{

			char str1[size];
			sprintf( str1, "%d", i + 1 );

			strcpy( input_name, "probability_density_cross_section_cb_" );
			strcat( input_name, str1 );
			strcat( input_name, "_" );
			copy_output_cross_section( direction, input_name );

			strcpy( input_name, "envelope_function_cross_section_cb_" );
			strcat( input_name, str1 );
			strcat( input_name, "_" );
			copy_output_cross_section( direction, input_name );

		}
		
		free( input_name );
		
	}
	else // 3D
	{
		
		int size = 160;
		
		char *input_name;
		input_name = malloc( sizeof(char)*size );


		// Valence band probability densities & envelope functions
		for (int i = 0; i < N_states_vb; i++)
		{

			char *str1 = malloc( sizeof(char)*size );
			sprintf( str1, "%d", i + 1 );

			strcpy( input_name, "probability_density_3D_vb_" );
			strcat( input_name, str1 );
			strcat( input_name, ".dat " );

			char *copy_command = malloc( sizeof(char)*size );
			strcpy( copy_command, "mv ");
			strcat( copy_command, input_name );
			strcat( copy_command, graph_folder );
			system( copy_command );


			strcpy( input_name, "envelope_function_3D_vb_" );
			strcat( input_name, str1 );
			strcat( input_name, ".dat " );

			strcpy( copy_command, "mv ");
			strcat( copy_command, input_name );
			strcat( copy_command, graph_folder );
			system( copy_command );


			free( str1 );
			free( copy_command );
			
		}

		// Conduction band probability densities & envelope functions
		for (int i = 0; i < N_states_cb; i++)
		{

			char *str1 = malloc( sizeof(char)*size );
			sprintf( str1, "%d", i + 1 );

			strcpy( input_name, "probability_density_3D_cb_" );
			strcat( input_name, str1 );
			strcat( input_name, ".dat " );

			char *copy_command = malloc( sizeof(char)*size );
			strcpy( copy_command, "mv ");
			strcat( copy_command, input_name );
			strcat( copy_command, graph_folder );
			system( copy_command );


			strcpy( input_name, "envelope_function_3D_cb_" );
			strcat( input_name, str1 );
			strcat( input_name, ".dat " );

			strcpy( copy_command, "mv ");
			strcat( copy_command, input_name );
			strcat( copy_command, graph_folder );
			system( copy_command );


			free( str1 );
			free( copy_command );

		}
		
		free( input_name );
		
	}

}

void copy_output_linescan( int direction, char *input_name )
{

	int size = 160;
	char *str1, *str2;
	str1 = malloc( sizeof(char)*size );
	str2 = malloc( sizeof(char)*size );

	strcpy( str1, input_name );

	switch( direction )
	{

		case 0: // Linescan along [100]

			strcpy( str2, "100.dat " );
			
			break;

		case 1: // Linescan along [010]

			strcpy( str2, "010.dat " );

			break;

		case 2: // Linescan along [001]

			strcpy( str2, "001.dat " );

			break;

		case 3: // Linescan along [110]

			strcpy( str2, "110.dat " );

			break;

		case 4: // Linescan along [101]

			strcpy( str2, "101.dat " );

			break;

		case 5: // Linescan along [011]

			strcpy( str2, "011.dat " );

			break;
			
		case 6: // Linescan along [111]

			strcpy( str2, "111.dat " );

			break;

	}

	strcat( str1, str2 );
	
	
	char *copy_command = malloc( sizeof(char)*size );
	strcpy( copy_command, "mv ");
	strcat( copy_command, str1 );
	strcat( copy_command, graph_folder );
	system( copy_command );
	
	
	free( str1 );
	free( str2 );
	free( copy_command );

}


void copy_output_cross_section( int direction, char *input_name )
{

	int size = 160;
	char *str1, *str2;
	str1 = malloc( sizeof(char)*size );
	str2 = malloc( sizeof(char)*size );

	strcpy( str1, input_name );

	switch( direction )
	{

		case 0: // Normal along [100]

			strcpy( str2, "100.dat " );	

			break;

		case 1: // Normal along [010]

			strcpy( str2, "010.dat " );

			break;

		case 2: // Normal along [001]

			strcpy( str2, "001.dat " );

			break;

		case 3: // Normal along [-110]

			strcpy( str2, "-110.dat " );

			break;

		case 4: // Normal along [110]

			strcpy( str2, "110.dat " );

			break;

		case 5: // Normal along [-101]

			strcpy( str2, "-101.dat " );

			break;

		case 6: // Normal along [0-11]

			strcpy( str2, "0-11.dat " );

			break;

		case 7: // Normal along [111]

			strcpy( str2, "111.dat " );

			break;

	}

	strcat( str1, str2 );

	char *copy_command = malloc( sizeof(char)*size );
	strcpy( copy_command, "mv ");
	strcat( copy_command, str1 );
	strcat( copy_command, graph_folder );
	system( copy_command );
	
	
	free( str1 );
	free( str2 );
	free( copy_command );

}
