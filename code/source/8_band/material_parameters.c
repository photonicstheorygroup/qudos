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
#include "material_parameters.h"

void material_parameters( int N_shape, double T, struct supercell_geometry *supercell, struct material_params *parameters )
{

	//////////////////////////////////////////////////////////////////////////////////////
	//// Set up multinary interpolation for the material system comprising each shape ////
	//////////////////////////////////////////////////////////////////////////////////////

	int idx_shape, bowing, binary[N_shape][4], ternary[N_shape][4]; // Arrays for storing material indices for parameter interpolation in quaternary alloys

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		switch( supercell[idx_shape].material )
		{

			////////////////////////////////////////////////////////////////////////
			//// 1. Alloys of the form A(x)B(1-x)C(y)D(1-y) - material_type = 0 ////
			////////////////////////////////////////////////////////////////////////

			case 0 :

				// Material system 0 = Al(x)Ga(1-x)P(y)As(1-y) - A = Al, B = Ga, C = P, D = As

				supercell[idx_shape].material_type = 0;

				binary[idx_shape][0] =  1; // AC  = AlP
				binary[idx_shape][1] =  2; // AD  = AlAs
				binary[idx_shape][2] =  6; // BC  = GaP
				binary[idx_shape][3] =  7; // BD  = GaAs

				ternary[idx_shape][0] =  0; // ABC = AlGaP
				ternary[idx_shape][1] =  1; // ABD = AlGaAs
				ternary[idx_shape][2] =  9; // ACD = AlPAs
				ternary[idx_shape][3] = 12; // BCD = GaPAs

				break;

			case 1 :

				// Material system 1 = Al(x)Ga(1-x)P(y)Sb(1-y) - A = Al, B = Ga, C = P, D = Sb

				supercell[idx_shape].material_type = 0;

				binary[idx_shape][0] =  1; // AC  = AlP
				binary[idx_shape][1] =  3; // AD  = AlSb
				binary[idx_shape][2] =  6; // BC  = GaP
				binary[idx_shape][3] =  8; // BD  = GaSb

				ternary[idx_shape][0] =  0; // ABC = AlGaP
				ternary[idx_shape][1] =  2; // ABD = AlGaSb
				ternary[idx_shape][2] = 10; // ACD = AlPSb
				ternary[idx_shape][3] = 13; // BCD = GaPSb

				break;


			case 2 :

				// Material system 2 = Al(x)Ga(1-x)As(y)Sb(1-y) - A = Al, B = Ga, C = As, D = Sb

				supercell[idx_shape].material_type = 0;

				binary[idx_shape][0] =  2; // AC  = AlAs
				binary[idx_shape][1] =  3; // AD  = AlSb
				binary[idx_shape][2] =  7; // BC  = GaAs
				binary[idx_shape][3] =  8; // BD  = GaSb

				ternary[idx_shape][0] =  1; // ABC = AlGaAs
				ternary[idx_shape][1] =  2; // ABD = AlGaSb
				ternary[idx_shape][2] = 11; // ACD = AlAsSb
				ternary[idx_shape][3] = 14; // BCD = GaAsSb

				break;


			case 3 :

				// Material system 3 = Ga(x)In(1-x)P(y)As(1-y) - A = Ga, B = In, C = P, D = As

				supercell[idx_shape].material_type = 0;

				binary[idx_shape][0] =  6; // AC  = GaP
				binary[idx_shape][1] =  7; // AD  = GaAs
				binary[idx_shape][2] = 11; // BC  = InP
				binary[idx_shape][3] = 12; // BD  = InAs

				ternary[idx_shape][0] =  3; // ABC = GaInP
				ternary[idx_shape][1] =  4; // ABD = GaInAs
				ternary[idx_shape][2] = 12; // ACD = GaPAs
				ternary[idx_shape][3] = 15; // BCD = InPAs

				break;


			case 4 :

				// Material system 4 = Ga(x)In(1-x)P(y)Sb(1-y) - A = Ga, B = In, C = P, D = Sb

				supercell[idx_shape].material_type = 0;

				binary[idx_shape][0] =  6; // AC = GaP
				binary[idx_shape][1] =  8; // AD = GaSb
				binary[idx_shape][2] = 11; // BC = InP
				binary[idx_shape][3] = 13; // BD = InSb

				ternary[idx_shape][0] =  3; // ABC = GaInP
				ternary[idx_shape][1] =  5; // ABD = GaInSb
				ternary[idx_shape][2] = 13; // ACD = GaPSb
				ternary[idx_shape][3] = 16; // BCD = InPSb

				break;


			case 5 :

				// Material system 5 = Ga(x)In(1-x)As(y)Sb(1-y) - A = Ga, B = In, C = As, D = Sb

				supercell[idx_shape].material_type = 0;

				binary[idx_shape][0] =  7; // AC = GaAs
				binary[idx_shape][1] =  8; // AD = GaSb
				binary[idx_shape][2] = 12; // BC = InAs
				binary[idx_shape][3] = 13; // BD = InSb

				ternary[idx_shape][0] =  4; // ABC = GaInAs
				ternary[idx_shape][1] =  5; // ABD = GaInSb
				ternary[idx_shape][2] = 14; // ACD = GaAsSb
				ternary[idx_shape][3] = 17; // BCD = InAsSb

				break;


			case 6 :

				// Material system 6 = Al(x)In(1-x)P(y)As(1-y) - A = Al, B = In, C = P, D = As

				supercell[idx_shape].material_type = 0;

				binary[idx_shape][0] =  1; // AC  = AlP
				binary[idx_shape][1] =  2; // AD  = AlAs
				binary[idx_shape][2] = 11; // BC  = InP
				binary[idx_shape][3] = 12; // BD  = InAs

				ternary[idx_shape][0] =  6; // ABC = AlInP
				ternary[idx_shape][1] =  7; // ABD = AlInAs
				ternary[idx_shape][2] =  9; // ACD = AlPAs
				ternary[idx_shape][3] = 15; // BCD = InPAs

				break;


			case 7 :

				// Material system 7 = Al(x)In(1-x)P(y)Sb(1-y) - A = Al, B = In, C = P, D = Sb

				supercell[idx_shape].material_type = 0;

				binary[idx_shape][0] =  1; // AC  = AlP
				binary[idx_shape][1] =  3; // AD  = AlSb
				binary[idx_shape][2] = 11; // BC  = InP
				binary[idx_shape][3] = 13; // BD  = InSb

				ternary[idx_shape][0] =  6; // ABC = AlInP
				ternary[idx_shape][1] =  8; // ABD = AlInSb
				ternary[idx_shape][2] = 10; // ACD = AlPSb
				ternary[idx_shape][3] = 16; // BCD = InPSb

				break;


			case 8 :

				// Material system 8 = Al(x)In(1-x)As(y)Sb(1-y) - A = Al, B = In, C = As, D = Sb

				supercell[idx_shape].material_type = 0;

				binary[idx_shape][0] =  2; // AC  = AlAs
				binary[idx_shape][1] =  3; // AD  = AlSb
				binary[idx_shape][2] = 12; // BC  = InAs
				binary[idx_shape][3] = 13; // BD  = InSb

				ternary[idx_shape][0] =  7; // ABC = AlInAs
				ternary[idx_shape][1] =  8; // ABD = AlInSb
				ternary[idx_shape][2] = 11; // ACD = AlAsSb
				ternary[idx_shape][3] = 17; // BCD = InAsSb

				break;


			/////////////////////////////////////////////////////////////////////
			//// 2. Alloys of the form A(x)B(y)C(1-x-y)D - material_type = 1 ////
			/////////////////////////////////////////////////////////////////////

			case 9 :

				// Material system 9 = Al(x)Ga(y)In(1-x-y)P - A = Al, B = Ga, C = In, D = P

				supercell[idx_shape].material_type = 1;

				binary[idx_shape][0] =  1; // AD  = AlP
				binary[idx_shape][1] =  6; // BD  = GaP
				binary[idx_shape][2] = 11; // CD  = InP
				binary[idx_shape][3] = -1; // Not used

				ternary[idx_shape][0] =  0; // ABD = AlGaP
				ternary[idx_shape][1] =  6; // ACD = AlInP
				ternary[idx_shape][2] =  3; // BCD = GaInP
				ternary[idx_shape][3] = -1; // Not used

				break;


			case 10 :

				// Material system 10 = Al(x)Ga(y)In(1-x-y)As - A = Al, B = Ga, C = In, D = As

				supercell[idx_shape].material_type = 1;

				binary[idx_shape][0] =  2; // AD  = AlAs
				binary[idx_shape][1] =  7; // BD  = GaAs
				binary[idx_shape][2] = 12; // CD  = InAs
				binary[idx_shape][3] = -1; // Not used

				ternary[idx_shape][0] =  1; // ABD = AlGaAs
				ternary[idx_shape][1] =  2; // ACD = AlInAs
				ternary[idx_shape][2] =  4; // BCD = GaInAs
				ternary[idx_shape][3] = -1; // Not used

				break;


			case 11 :

				// Material system 11 = Al(x)Ga(y)In(1-x-y)Sb - A = Al, B = Ga, C = In, D = Sb

				supercell[idx_shape].material_type = 1;

				binary[idx_shape][0] =  3; // AD  = AlSb
				binary[idx_shape][1] =  8; // BD  = GaSb
				binary[idx_shape][2] = 13; // CD  = InSb
				binary[idx_shape][3] = -1; // Not used

				ternary[idx_shape][0] =  2; // ABD = AlGaSb
				ternary[idx_shape][1] =  8; // ACD = AlInSb
				ternary[idx_shape][2] =  5; // BCD = GaInSb
				ternary[idx_shape][3] = -1; // Not used

				break;


			/////////////////////////////////////////////////////////////////////
			//// 3. Alloys of the form AB(x)C(y)D(1-x-y) - material_type = 2 ////
			/////////////////////////////////////////////////////////////////////

			case 12 :

				// Material system 12 = AlP(x)As(y)Sb(1-x-y) - A = Al, B = P, C = As, D = Sb

				supercell[idx_shape].material_type = 2;

				binary[idx_shape][0] =  1; // AB  = AlP
				binary[idx_shape][1] =  2; // AC  = AlAs
				binary[idx_shape][2] =  3; // AD  = AlSb
				binary[idx_shape][3] = -1; // Not used

				ternary[idx_shape][0] =  9; // ABC = AlPAs
				ternary[idx_shape][1] = 10; // ABD = AlPSb
				ternary[idx_shape][2] = 11; // ACD = AlAsSb
				ternary[idx_shape][3] = -1; // Not used

				break;


			case 13 :

				// Material system 13 = GaP(x)As(y)Sb(1-x-y) - A = Ga, B = P, C = As, D = Sb

				supercell[idx_shape].material_type = 2;

				binary[idx_shape][0] =  6; // AB  = GaP
				binary[idx_shape][1] =  7; // AC  = GaAs
				binary[idx_shape][2] =  8; // AD  = GaSb
				binary[idx_shape][3] = -1; // Not used

				ternary[idx_shape][0] = 12; // ABC = GaPAs
				ternary[idx_shape][1] = 13; // ABD = GaPSb
				ternary[idx_shape][2] = 14; // ACD = GaAsSb
				ternary[idx_shape][3] = -1; // Not used

				break;


			case 14 :

				// Material system 14 = InP(x)As(y)Sb(1-x-y) - A = In, B = P, C = As, D = Sb

				supercell[idx_shape].material_type = 2;

				binary[idx_shape][0] = 11; // AB = InP
				binary[idx_shape][1] = 12; // AC = InAs
				binary[idx_shape][2] = 13; // AD = InSb
				binary[idx_shape][3] = -1; // Not used

				ternary[idx_shape][0] = 15; // ABC = InPAs
				ternary[idx_shape][1] = 16; // ABD = InPSb
				ternary[idx_shape][2] = 17; // ACD = InAsSb
				ternary[idx_shape][3] = -1; // Not used

				break;

		} // End of switch statement


		// If the requested material is not one of the cases handled in the above switch statement then the requested material system is not defined

		if( ( supercell[idx_shape].material < 0 ) || ( supercell[idx_shape].material > 14 ) )
		{

			printf( "\n\tRequested material system not recognised - the calculation has been aborted\n" );

			MPI_Abort( MPI_COMM_WORLD, 1 );

		}

	}

	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Calculate the modified inverse electron effective mass and Luttinger parameters for each binary compound using the original band gaps and spin-orbit splittings ////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	int idx_binary;

	int  binary_size = sizeof(     Eg      )/sizeof( double ); // Size of arrays in material_parameters.h containing binary  material parameters
	int ternary_size = sizeof( bow_Eg_tern )/sizeof( double ); // Size of arrays in material_parameters.h containing ternary material parameters
	double sc[binary_size];
	
	
	double g1[binary_size], g2[binary_size], g3[binary_size];

	for ( idx_binary = 0; idx_binary < binary_size; idx_binary++ )
	{

		sc[idx_binary] = 1.0/ mc[idx_binary] - ( Ep[idx_binary]/  3.0 )*( 2.0/Eg[idx_binary] + 1.0/( Eg[idx_binary] + Dso[idx_binary] ) );
		
		g1[idx_binary] =     gL1[idx_binary] -   Ep[idx_binary]/( 3.0  *      Eg[idx_binary] + Dso[idx_binary]);
		g2[idx_binary] =     gL2[idx_binary] -   Ep[idx_binary]/( 6.0  *      Eg[idx_binary] + 2.0 * Dso[idx_binary] );
		g3[idx_binary] =     gL3[idx_binary] -   Ep[idx_binary]/( 6.0  *      Eg[idx_binary] + 2.0 * Dso[idx_binary] );

	}
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Adjust the binary compound band gaps for the user-specified temperature using the Varshni relation ////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	for( idx_binary = 0; idx_binary < binary_size; idx_binary++ )
	{
		
		alpha_T[idx_binary] *= 0.001; // Convert the alpha Varshni parameter from meV/K to eV/K

		if( T > 0.001 ) // Adjust band gaps only for temperatures > 1 mK to avoid dividing by zero in various functions
		{

			if( idx_binary != 6 ) // Deal with GaP separately - see material_parameters.h
			{

				Eg[idx_binary] -= ( alpha_T[idx_binary]*T*T )/( T + beta_T[idx_binary] );

			}
			else
			{

				Eg[idx_binary] += 0.1081*( 1.0 - ( cosh( 164.0/T )/sinh( 164.0/T ) ) );

			}

		}

	}



	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Interpolate strain-related parameters for each shape (material) in the calculational supercell ////
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	struct interpolate_parameters alloy_parameters; // Struct to pass binary compound parameters and bowing parameters to interpolate_quaternary( )

	// bowing = 0 ==> Linear interpolation
	// bowing = 1 ==> Interpolation with ternary and quaternary bowing parameters

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		// Lattice constant, a

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = a[binary[idx_shape][0]];
				alloy_parameters.P_AD = a[binary[idx_shape][1]];
				alloy_parameters.P_BC = a[binary[idx_shape][2]];
				alloy_parameters.P_BD = a[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = a[binary[idx_shape][0]];
				alloy_parameters.P_BD = a[binary[idx_shape][1]];
				alloy_parameters.P_CD = a[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = a[binary[idx_shape][0]];
				alloy_parameters.P_AC = a[binary[idx_shape][1]];
				alloy_parameters.P_AD = a[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].a = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Elastic constant, C11

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = C11[binary[idx_shape][0]];
				alloy_parameters.P_AD = C11[binary[idx_shape][1]];
				alloy_parameters.P_BC = C11[binary[idx_shape][2]];
				alloy_parameters.P_BD = C11[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = C11[binary[idx_shape][0]];
				alloy_parameters.P_BD = C11[binary[idx_shape][1]];
				alloy_parameters.P_CD = C11[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = C11[binary[idx_shape][0]];
				alloy_parameters.P_AC = C11[binary[idx_shape][1]];
				alloy_parameters.P_AD = C11[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].C11 = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Elastic constant, C12

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = C12[binary[idx_shape][0]];
				alloy_parameters.P_AD = C12[binary[idx_shape][1]];
				alloy_parameters.P_BC = C12[binary[idx_shape][2]];
				alloy_parameters.P_BD = C12[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = C12[binary[idx_shape][0]];
				alloy_parameters.P_BD = C12[binary[idx_shape][1]];
				alloy_parameters.P_CD = C12[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = C12[binary[idx_shape][0]];
				alloy_parameters.P_AC = C12[binary[idx_shape][1]];
				alloy_parameters.P_AD = C12[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].C12 = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Elastic constant, C44

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = C44[binary[idx_shape][0]];
				alloy_parameters.P_AD = C44[binary[idx_shape][1]];
				alloy_parameters.P_BC = C44[binary[idx_shape][2]];
				alloy_parameters.P_BD = C44[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = C44[binary[idx_shape][0]];
				alloy_parameters.P_BD = C44[binary[idx_shape][1]];
				alloy_parameters.P_CD = C44[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = C44[binary[idx_shape][0]];
				alloy_parameters.P_AC = C44[binary[idx_shape][1]];
				alloy_parameters.P_AD = C44[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].C44 = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// First order piezoelectric coefficient, e14

		bowing = 0;

		if( include_piezo == 1 )
		{

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = e14[binary[idx_shape][0]];
				alloy_parameters.P_AD = e14[binary[idx_shape][1]];
				alloy_parameters.P_BC = e14[binary[idx_shape][2]];
				alloy_parameters.P_BD = e14[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = e14[binary[idx_shape][0]];
				alloy_parameters.P_BD = e14[binary[idx_shape][1]];
				alloy_parameters.P_CD = e14[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = e14[binary[idx_shape][0]];
				alloy_parameters.P_AC = e14[binary[idx_shape][1]];
				alloy_parameters.P_AD = e14[binary[idx_shape][2]];

				break;

		}

		}

		else
		{

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = e14_calc[binary[idx_shape][0]];
				alloy_parameters.P_AD = e14_calc[binary[idx_shape][1]];
				alloy_parameters.P_BC = e14_calc[binary[idx_shape][2]];
				alloy_parameters.P_BD = e14_calc[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = e14_calc[binary[idx_shape][0]];
				alloy_parameters.P_BD = e14_calc[binary[idx_shape][1]];
				alloy_parameters.P_CD = e14_calc[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = e14_calc[binary[idx_shape][0]];
				alloy_parameters.P_AC = e14_calc[binary[idx_shape][1]];
				alloy_parameters.P_AD = e14_calc[binary[idx_shape][2]];

				break;

		}

		}

	     

		parameters[idx_shape].e14 = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Second order piezoelectric coefficient, B114

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = B114[binary[idx_shape][0]];
				alloy_parameters.P_AD = B114[binary[idx_shape][1]];
				alloy_parameters.P_BC = B114[binary[idx_shape][2]];
				alloy_parameters.P_BD = B114[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = B114[binary[idx_shape][0]];
				alloy_parameters.P_BD = B114[binary[idx_shape][1]];
				alloy_parameters.P_CD = B114[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = B114[binary[idx_shape][0]];
				alloy_parameters.P_AC = B114[binary[idx_shape][1]];
				alloy_parameters.P_AD = B114[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].B114 = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Second order piezoelectric coefficient, B124

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = B124[binary[idx_shape][0]];
				alloy_parameters.P_AD = B124[binary[idx_shape][1]];
				alloy_parameters.P_BC = B124[binary[idx_shape][2]];
				alloy_parameters.P_BD = B124[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = B124[binary[idx_shape][0]];
				alloy_parameters.P_BD = B124[binary[idx_shape][1]];
				alloy_parameters.P_CD = B124[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = B124[binary[idx_shape][0]];
				alloy_parameters.P_AC = B124[binary[idx_shape][1]];
				alloy_parameters.P_AD = B124[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].B124 = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Second order piezoelectric coefficient, B156

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = B156[binary[idx_shape][0]];
				alloy_parameters.P_AD = B156[binary[idx_shape][1]];
				alloy_parameters.P_BC = B156[binary[idx_shape][2]];
				alloy_parameters.P_BD = B156[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = B156[binary[idx_shape][0]];
				alloy_parameters.P_BD = B156[binary[idx_shape][1]];
				alloy_parameters.P_CD = B156[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = B156[binary[idx_shape][0]];
				alloy_parameters.P_AC = B156[binary[idx_shape][1]];
				alloy_parameters.P_AD = B156[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].B156 = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Conduction band edge hydrostatic deformation potential, ac

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = ac[binary[idx_shape][0]];
				alloy_parameters.P_AD = ac[binary[idx_shape][1]];
				alloy_parameters.P_BC = ac[binary[idx_shape][2]];
				alloy_parameters.P_BD = ac[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = ac[binary[idx_shape][0]];
				alloy_parameters.P_BD = ac[binary[idx_shape][1]];
				alloy_parameters.P_CD = ac[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = ac[binary[idx_shape][0]];
				alloy_parameters.P_AC = ac[binary[idx_shape][1]];
				alloy_parameters.P_AD = ac[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].ac = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Valence band edge hydrostatic deformation potential, av

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = av[binary[idx_shape][0]];
				alloy_parameters.P_AD = av[binary[idx_shape][1]];
				alloy_parameters.P_BC = av[binary[idx_shape][2]];
				alloy_parameters.P_BD = av[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = av[binary[idx_shape][0]];
				alloy_parameters.P_BD = av[binary[idx_shape][1]];
				alloy_parameters.P_CD = av[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = av[binary[idx_shape][0]];
				alloy_parameters.P_AC = av[binary[idx_shape][1]];
				alloy_parameters.P_AD = av[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].av = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Valence band edge axial deformation potential, b

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = b[binary[idx_shape][0]];
				alloy_parameters.P_AD = b[binary[idx_shape][1]];
				alloy_parameters.P_BC = b[binary[idx_shape][2]];
				alloy_parameters.P_BD = b[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = b[binary[idx_shape][0]];
				alloy_parameters.P_BD = b[binary[idx_shape][1]];
				alloy_parameters.P_CD = b[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = b[binary[idx_shape][0]];
				alloy_parameters.P_AC = b[binary[idx_shape][1]];
				alloy_parameters.P_AD = b[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].b = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Valence band edge shear deformation potential, d

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = d[binary[idx_shape][0]];
				alloy_parameters.P_AD = d[binary[idx_shape][1]];
				alloy_parameters.P_BC = d[binary[idx_shape][2]];
				alloy_parameters.P_BD = d[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = d[binary[idx_shape][0]];
				alloy_parameters.P_BD = d[binary[idx_shape][1]];
				alloy_parameters.P_CD = d[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = d[binary[idx_shape][0]];
				alloy_parameters.P_AC = d[binary[idx_shape][1]];
				alloy_parameters.P_AD = d[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].d = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Refractive index, nr

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = nr[binary[idx_shape][0]];
				alloy_parameters.P_AD = nr[binary[idx_shape][1]];
				alloy_parameters.P_BC = nr[binary[idx_shape][2]];
				alloy_parameters.P_BD = nr[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = nr[binary[idx_shape][0]];
				alloy_parameters.P_BD = nr[binary[idx_shape][1]];
				alloy_parameters.P_CD = nr[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = nr[binary[idx_shape][0]];
				alloy_parameters.P_AC = nr[binary[idx_shape][1]];
				alloy_parameters.P_AD = nr[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].nr = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );

	}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Ternary bowing parameters for the average energy of the three highest energy valence bands ////
	////////////////////////////////////////////////////////////////////////////////////////////////////

	// For the calculation of these bowing parameters the reference lattice constant is chosen to be that of the barrier (matrix) material: parameters[0].a

	double bow_Ev_avg_tern[ternary_size];

	bow_Ev_avg_tern[0]  = ( 3.0/parameters[0].a )*( av[1]  - av[6]  )*( a[1]  - a[6]  ); // AlGaP
	bow_Ev_avg_tern[1]  = ( 3.0/parameters[0].a )*( av[2]  - av[7]  )*( a[2]  - a[7]  ); // AlGaAs
	bow_Ev_avg_tern[2]  = ( 3.0/parameters[0].a )*( av[3]  - av[8]  )*( a[3]  - a[8]  ); // AlGaSb

	bow_Ev_avg_tern[3]  = ( 3.0/parameters[0].a )*( av[6]  - av[11] )*( a[6]  - a[11] ); // GaInP
	bow_Ev_avg_tern[4]  = ( 3.0/parameters[0].a )*( av[7]  - av[12] )*( a[7]  - a[12] ); // GaInAs
	bow_Ev_avg_tern[5]  = ( 3.0/parameters[0].a )*( av[8]  - av[13] )*( a[8]  - a[13] ); // GaInSb

	bow_Ev_avg_tern[6]  = ( 3.0/parameters[0].a )*( av[1]  - av[11] )*( a[1]  - a[11] ); // AlInP
	bow_Ev_avg_tern[7]  = ( 3.0/parameters[0].a )*( av[2]  - av[12] )*( a[2]  - a[12] ); // AlInAs
	bow_Ev_avg_tern[8]  = ( 3.0/parameters[0].a )*( av[3]  - av[13] )*( a[3]  - a[13] ); // AlInSb

	bow_Ev_avg_tern[9]  = ( 3.0/parameters[0].a )*( av[1]  - av[2]  )*( a[1]  - a[2]  ); // AlPAs
	bow_Ev_avg_tern[10] = ( 3.0/parameters[0].a )*( av[1]  - av[3]  )*( a[1]  - a[3]  ); // AlPSb
	bow_Ev_avg_tern[11] = ( 3.0/parameters[0].a )*( av[2]  - av[3]  )*( a[2]  - a[3]  ); // AlAsSb

	bow_Ev_avg_tern[12] = ( 3.0/parameters[0].a )*( av[6]  - av[7]  )*( a[6]  - a[7]  ); // GaPAs
	bow_Ev_avg_tern[13] = ( 3.0/parameters[0].a )*( av[6]  - av[8]  )*( a[6]  - a[8]  ); // GaPSb
	bow_Ev_avg_tern[14] = ( 3.0/parameters[0].a )*( av[7]  - av[8]  )*( a[7]  - a[8]  ); // GaAsSb

	bow_Ev_avg_tern[15] = ( 3.0/parameters[0].a )*( av[11] - av[12] )*( a[11] - a[12] ); // InPAs
	bow_Ev_avg_tern[16] = ( 3.0/parameters[0].a )*( av[11] - av[13] )*( a[11] - a[13] ); // InPSb
	bow_Ev_avg_tern[17] = ( 3.0/parameters[0].a )*( av[12] - av[13] )*( a[12] - a[13] ); // InAsSb


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Interpolate band structure parameters for each shape (material) in the calculational supercell ////
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	// bowing = 0 ==> Linear interpolation
	// bowing = 1 ==> Interpolation with ternary and quaternary bowing parameters

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

		// Average energy of the three highest energy valence bands - i.e. the light-hole, heavy-hole and spin-split-off-hole bands

		bowing = 1;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC   = Ev_avg[binary[idx_shape][0]];
				alloy_parameters.P_AD   = Ev_avg[binary[idx_shape][1]];
				alloy_parameters.P_BC   = Ev_avg[binary[idx_shape][2]];
				alloy_parameters.P_BD   = Ev_avg[binary[idx_shape][3]];

				alloy_parameters.B_ABC  = bow_Ev_avg_tern[ternary[idx_shape][0]];
				alloy_parameters.B_ABD  = bow_Ev_avg_tern[ternary[idx_shape][1]];
				alloy_parameters.B_ACD  = bow_Ev_avg_tern[ternary[idx_shape][2]];
				alloy_parameters.B_BCD  = bow_Ev_avg_tern[ternary[idx_shape][3]];

				alloy_parameters.B_ABCD = 0.0;

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD   = Ev_avg[binary[idx_shape][0]];
				alloy_parameters.P_BD   = Ev_avg[binary[idx_shape][1]];
				alloy_parameters.P_CD   = Ev_avg[binary[idx_shape][2]];

				alloy_parameters.B_ABD  = bow_Ev_avg_tern[ternary[idx_shape][0]];
				alloy_parameters.B_ACD  = bow_Ev_avg_tern[ternary[idx_shape][1]];
				alloy_parameters.B_BCD  = bow_Ev_avg_tern[ternary[idx_shape][2]];

				alloy_parameters.B_ABCD = 0.0;

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB   = Ev_avg[binary[idx_shape][0]];
				alloy_parameters.P_AC   = Ev_avg[binary[idx_shape][1]];
				alloy_parameters.P_AD   = Ev_avg[binary[idx_shape][2]];

				alloy_parameters.B_ABC  = bow_Ev_avg_tern[ternary[idx_shape][0]];
				alloy_parameters.B_ABD  = bow_Ev_avg_tern[ternary[idx_shape][1]];
				alloy_parameters.B_ACD  = bow_Ev_avg_tern[ternary[idx_shape][2]];

				alloy_parameters.B_ABCD = 0.0;

				break;

		}

		parameters[idx_shape].Ev_avg = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Band gap

		bowing = 1;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC   = Eg[binary[idx_shape][0]];
				alloy_parameters.P_AD   = Eg[binary[idx_shape][1]];
				alloy_parameters.P_BC   = Eg[binary[idx_shape][2]];
				alloy_parameters.P_BD   = Eg[binary[idx_shape][3]];

				alloy_parameters.B_ABC  = bow_Eg_tern[ternary[idx_shape][0]];
				alloy_parameters.B_ABD  = bow_Eg_tern[ternary[idx_shape][1]];
				alloy_parameters.B_ACD  = bow_Eg_tern[ternary[idx_shape][2]];
				alloy_parameters.B_BCD  = bow_Eg_tern[ternary[idx_shape][3]];

				alloy_parameters.B_ABCD = bow_Eg_quat[supercell[idx_shape].material];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD   = Eg[binary[idx_shape][0]];
				alloy_parameters.P_BD   = Eg[binary[idx_shape][1]];
				alloy_parameters.P_CD   = Eg[binary[idx_shape][2]];

				alloy_parameters.B_ABD  = bow_Eg_tern[ternary[idx_shape][0]];
				alloy_parameters.B_ACD  = bow_Eg_tern[ternary[idx_shape][1]];
				alloy_parameters.B_BCD  = bow_Eg_tern[ternary[idx_shape][2]];

				alloy_parameters.B_ABCD = bow_Eg_quat[supercell[idx_shape].material];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB   = Eg[binary[idx_shape][0]];
				alloy_parameters.P_AC   = Eg[binary[idx_shape][1]];
				alloy_parameters.P_AD   = Eg[binary[idx_shape][2]];

				alloy_parameters.B_ABC  = bow_Eg_tern[ternary[idx_shape][0]];
				alloy_parameters.B_ABD  = bow_Eg_tern[ternary[idx_shape][1]];
				alloy_parameters.B_ACD  = bow_Eg_tern[ternary[idx_shape][2]];

				alloy_parameters.B_ABCD = bow_Eg_quat[supercell[idx_shape].material];

				break;

		}

		parameters[idx_shape].Eg = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Spin-orbit-splitting energy

		bowing = 1;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC   = Dso[binary[idx_shape][0]];
				alloy_parameters.P_AD   = Dso[binary[idx_shape][1]];
				alloy_parameters.P_BC   = Dso[binary[idx_shape][2]];
				alloy_parameters.P_BD   = Dso[binary[idx_shape][3]];

				alloy_parameters.B_ABC  = bow_Dso_tern[ternary[idx_shape][0]];
				alloy_parameters.B_ABD  = bow_Dso_tern[ternary[idx_shape][1]];
				alloy_parameters.B_ACD  = bow_Dso_tern[ternary[idx_shape][2]];
				alloy_parameters.B_BCD  = bow_Dso_tern[ternary[idx_shape][3]];

				alloy_parameters.B_ABCD = bow_Dso_quat[supercell[idx_shape].material];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD   = Dso[binary[idx_shape][0]];
				alloy_parameters.P_BD   = Dso[binary[idx_shape][1]];
				alloy_parameters.P_CD   = Dso[binary[idx_shape][2]];

				alloy_parameters.B_ABD  = bow_Dso_tern[ternary[idx_shape][0]];
				alloy_parameters.B_ACD  = bow_Dso_tern[ternary[idx_shape][1]];
				alloy_parameters.B_BCD  = bow_Dso_tern[ternary[idx_shape][2]];

				alloy_parameters.B_ABCD = bow_Dso_quat[supercell[idx_shape].material];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB   = Dso[binary[idx_shape][0]];
				alloy_parameters.P_AC   = Dso[binary[idx_shape][1]];
				alloy_parameters.P_AD   = Dso[binary[idx_shape][2]];

				alloy_parameters.B_ABC  = bow_Dso_tern[ternary[idx_shape][0]];
				alloy_parameters.B_ABD  = bow_Dso_tern[ternary[idx_shape][1]];
				alloy_parameters.B_ACD  = bow_Dso_tern[ternary[idx_shape][2]];

				alloy_parameters.B_ABCD = bow_Dso_quat[supercell[idx_shape].material];

				break;

		}

		parameters[idx_shape].Dso = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Modified inverse electron effective mass

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = sc[binary[idx_shape][0]];
				alloy_parameters.P_AD = sc[binary[idx_shape][1]];
				alloy_parameters.P_BC = sc[binary[idx_shape][2]];
				alloy_parameters.P_BD = sc[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = sc[binary[idx_shape][0]];
				alloy_parameters.P_BD = sc[binary[idx_shape][1]];
				alloy_parameters.P_CD = sc[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = sc[binary[idx_shape][0]];
				alloy_parameters.P_AC = sc[binary[idx_shape][1]];
				alloy_parameters.P_AD = sc[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].sc = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Kane parameter

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = Ep[binary[idx_shape][0]];
				alloy_parameters.P_AD = Ep[binary[idx_shape][1]];
				alloy_parameters.P_BC = Ep[binary[idx_shape][2]];
				alloy_parameters.P_BD = Ep[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = Ep[binary[idx_shape][0]];
				alloy_parameters.P_BD = Ep[binary[idx_shape][1]];
				alloy_parameters.P_CD = Ep[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = Ep[binary[idx_shape][0]];
				alloy_parameters.P_AC = Ep[binary[idx_shape][1]];
				alloy_parameters.P_AD = Ep[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].Ep = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Interband momentum matrix element (bulk) - compute directly using interpolated Kane parameter

		parameters[idx_shape].P = sqrt( 0.5*HBSQM*parameters[idx_shape].Ep );


		// First modified Luttinger parameter

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = g1[binary[idx_shape][0]];
				alloy_parameters.P_AD = g1[binary[idx_shape][1]];
				alloy_parameters.P_BC = g1[binary[idx_shape][2]];
				alloy_parameters.P_BD = g1[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = g1[binary[idx_shape][0]];
				alloy_parameters.P_BD = g1[binary[idx_shape][1]];
				alloy_parameters.P_CD = g1[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = g1[binary[idx_shape][0]];
				alloy_parameters.P_AC = g1[binary[idx_shape][1]];
				alloy_parameters.P_AD = g1[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].g1 = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Second modified Luttinger parameter

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = g2[binary[idx_shape][0]];
				alloy_parameters.P_AD = g2[binary[idx_shape][1]];
				alloy_parameters.P_BC = g2[binary[idx_shape][2]];
				alloy_parameters.P_BD = g2[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = g2[binary[idx_shape][0]];
				alloy_parameters.P_BD = g2[binary[idx_shape][1]];
				alloy_parameters.P_CD = g2[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = g2[binary[idx_shape][0]];
				alloy_parameters.P_AC = g2[binary[idx_shape][1]];
				alloy_parameters.P_AD = g2[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].g2 = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


		// Third modified Luttinger parameter

		bowing = 0;

		switch( supercell[idx_shape].material_type )
		{

			case 0 :

				// A(x)B(1-x)C(y)D(1-y)

				alloy_parameters.P_AC = g3[binary[idx_shape][0]];
				alloy_parameters.P_AD = g3[binary[idx_shape][1]];
				alloy_parameters.P_BC = g3[binary[idx_shape][2]];
				alloy_parameters.P_BD = g3[binary[idx_shape][3]];

				break;


			case 1 :

				// A(x)B(y)C(1-x-y)D

				alloy_parameters.P_AD = g3[binary[idx_shape][0]];
				alloy_parameters.P_BD = g3[binary[idx_shape][1]];
				alloy_parameters.P_CD = g3[binary[idx_shape][2]];

				break;


			case 2 :

				// AB(x)C(y)D(1-x-y)

				alloy_parameters.P_AB = g3[binary[idx_shape][0]];
				alloy_parameters.P_AC = g3[binary[idx_shape][1]];
				alloy_parameters.P_AD = g3[binary[idx_shape][2]];

				break;

		}

		parameters[idx_shape].g3 = interpolate_quaternary( bowing, supercell[idx_shape], alloy_parameters );


	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// Set the zero of energy and calculate the unstrained band edge energies in each shape (material) ////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	for( idx_shape = 0; idx_shape < N_shape; idx_shape++ )
	{

	    // Set the zero of energy at the unstrained valence band edge of the barrier material - i.e. the "matrix" material, idx_shape = 0

	    parameters[idx_shape].E_hh  = parameters[idx_shape].Ev_avg + parameters[idx_shape].Dso/3.0 - ( parameters[0].Ev_avg + parameters[0].Dso/3.0 );

		parameters[idx_shape].E_so = parameters[idx_shape].E_hh - parameters[idx_shape].Dso;
		parameters[idx_shape].E_lh = parameters[idx_shape].E_hh;
	    parameters[idx_shape].E_cb = parameters[idx_shape].E_hh + parameters[idx_shape].Eg;

	}

}
