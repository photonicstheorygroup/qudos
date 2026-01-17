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

void write_supercell( int N_shape, double T, double Fx, double Fy, double Fz, struct supercell_geometry *supercell, struct material_params *parameters )
{

	FILE *fp1;
	fp1 = fopen( "supercell.dat", "w" );

	fprintf( fp1, "Note: The zero of energy is taken at the unstrained valence band edge of the supercell matrix (barrier) material\n" );
	fprintf( fp1, "\n" );
	fprintf( fp1, "Temperature = %13e K\n", T );
	fprintf( fp1, "\n" );
	fprintf( fp1, "Applied electric field = ( %13e, %13e, %13e ) kV/cm \n", 1.0e+05*Fx, 1.0e+05*Fy, 1.0e+05*Fz );
	fprintf( fp1, "\n" );


	switch( include_strain )
	{

		case 0:

			fprintf( fp1, "Strain is not included in this calculation\n" );

			break;

		case 1:

			fprintf( fp1, "Strain is included in this calculation, using the " );

			if( isotropic == 0 )
			{

				fprintf( fp1, "full anisotropic elastic tensor\n" );

			}
			else
			{

				fprintf( fp1, "isotropic approximation to the elastic tensor" );

			}

			switch( include_shape_elastic_constants )
			{

				case 0:

					fprintf( fp1, "\nUsing elastic constants of matrix material\n" );

					break;

				case 1:

					fprintf( fp1, "\nUsing position-dependent elastic constants, with the %d order correction", shape_elastic_constants_order );

					break;

			}

			break;

	}

	fprintf( fp1, "\n" );


	switch( include_piezo )
	{

		case 0:

			fprintf( fp1, "The piezoelectric potential is not included in this calculation\n" );

			break;

		case 1:

			fprintf( fp1, "Including first order piezoelectric potential only\n" );

			break;

		case 2:

			fprintf( fp1, "Including first and second order piezoelectric potential\n" );

			break;

		case 3:

			fprintf( fp1, "Including second order piezoelectric potential only\n" );

			break;

	}

	fprintf( fp1, "\n" );


	///////////////////////////////////////////////////////////
	//// Write out information for matrix/barrier material ////
	///////////////////////////////////////////////////////////

	// 1. Material system and alloy composition

	fprintf( fp1, "Shape 0: Supercell matrix (barrier)\n" );
	fprintf( fp1, "\n" );

	switch( supercell[0].material )
	{

		// Alloys of the form A(x)B(1-x)C(y)D(1-y)

		case 0: fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)P(%.4f)As(%.4f) \n", supercell[0].x, 1.0 - supercell[0].x, supercell[0].y, 1.0 - supercell[0].y ); break;
		case 1: fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)P(%.4f)Sb(%.4f) \n", supercell[0].x, 1.0 - supercell[0].x, supercell[0].y, 1.0 - supercell[0].y ); break;
		case 2: fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)As(%.4f)Sb(%.4f)\n", supercell[0].x, 1.0 - supercell[0].x, supercell[0].y, 1.0 - supercell[0].y ); break;

		case 3: fprintf( fp1, "\tMaterial = Ga(%.4f)In(%.4f)P(%.4f)As(%.4f) \n", supercell[0].x, 1.0 - supercell[0].x, supercell[0].y, 1.0 - supercell[0].y ); break;
		case 4: fprintf( fp1, "\tMaterial = Ga(%.4f)In(%.4f)P(%.4f)Sb(%.4f) \n", supercell[0].x, 1.0 - supercell[0].x, supercell[0].y, 1.0 - supercell[0].y ); break;
		case 5: fprintf( fp1, "\tMaterial = Ga(%.4f)In(%.4f)As(%.4f)Sb(%.4f)\n", supercell[0].x, 1.0 - supercell[0].x, supercell[0].y, 1.0 - supercell[0].y ); break;

		case 6: fprintf( fp1, "\tMaterial = Al(%.4f)In(%.4f)P(%.4f)As(%.4f) \n", supercell[0].x, 1.0 - supercell[0].x, supercell[0].y, 1.0 - supercell[0].y ); break;
		case 7: fprintf( fp1, "\tMaterial = Al(%.4f)In(%.4f)P(%.4f)Sb(%.4f) \n", supercell[0].x, 1.0 - supercell[0].x, supercell[0].y, 1.0 - supercell[0].y ); break;
		case 8: fprintf( fp1, "\tMaterial = Al(%.4f)In(%.4f)As(%.4f)Sb(%.4f)\n", supercell[0].x, 1.0 - supercell[0].x, supercell[0].y, 1.0 - supercell[0].y ); break;


		// Alloys of the form A(x)B(y)C(1-x-y)D

		case 9 : fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)In(%.4f)P \n", supercell[0].x, supercell[0].y, 1.0 - supercell[0].x - supercell[0].y ); break;
		case 10: fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)In(%.4f)As\n", supercell[0].x, supercell[0].y, 1.0 - supercell[0].x - supercell[0].y ); break;
		case 11: fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)In(%.4f)Sb\n", supercell[0].x, supercell[0].y, 1.0 - supercell[0].x - supercell[0].y ); break;


		// Alloys of the form AB(x)C(y)D(1-x-y)

		case 12: fprintf( fp1, "\tMaterial = AlP(%.4f)As(%.4f)Sb(%.4f)\n", supercell[0].x, supercell[0].y, 1.0 - supercell[0].x - supercell[0].y ); break;
		case 14: fprintf( fp1, "\tMaterial = InP(%.4f)As(%.4f)Sb(%.4f)\n", supercell[0].x, supercell[0].y, 1.0 - supercell[0].x - supercell[0].y ); break;
		case 13: fprintf( fp1, "\tMaterial = GaP(%.4f)As(%.4f)Sb(%.4f)\n", supercell[0].x, supercell[0].y, 1.0 - supercell[0].x - supercell[0].y ); break;

	}

	fprintf( fp1, "\n" );


	// 2. Geometrical information - location and size

	if( supercell[0].shape == 0 )
	{

		fprintf( fp1, "\tShape = Cuboid\n" );

		fprintf( fp1, "\n" );
		fprintf( fp1, "\tCentre (x, y, z)         [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[0].Ox, 0.1*supercell[0].Oy, 0.1*supercell[0].Oz );
		fprintf( fp1, "\tSupercell length along x [nm] = %.4f\n", 0.1*supercell[0].L1 );
		fprintf( fp1, "\tSupercell length along y [nm] = %.4f\n", 0.1*supercell[0].L2 );
		fprintf( fp1, "\tSupercell length along z [nm] = %.4f\n", 0.1*supercell[0].L3 );
		fprintf( fp1, "\n" );
		
		fprintf( fp1, "\tStrain calculation supercell length along x [nm] = %.4f\n", 0.1*supercell[0].SL1 );
		fprintf( fp1, "\tStrain calculation supercell length along y [nm] = %.4f\n", 0.1*supercell[0].SL2 );
		fprintf( fp1, "\tStrain calculation supercell length along z [nm] = %.4f\n", 0.1*supercell[0].SL3 );
		fprintf( fp1, "\n" );

		fprintf( fp1, "\tCoulomb matrix element calculation supercell length along x [nm] = %.4f\n", 0.1*supercell[0].CL1 );
		fprintf( fp1, "\tCoulomb matrix element calculation supercell length along y [nm] = %.4f\n", 0.1*supercell[0].CL2 );
		fprintf( fp1, "\tCoulomb matrix element calculation supercell length along z [nm] = %.4f\n", 0.1*supercell[0].CL3 );
		fprintf( fp1, "\n" );

	}
	else
	{

		printf( "\n" );
		printf( "Only cuboid-shaped supercells are currently supported by this software - the calculation has been aborted\n" );
		printf( "\n" );
		printf( "Please edit supercell.in so that QD shape = 0 for the 1st shape and try again\n" );
		printf( "\n" );

		MPI_Abort( MPI_COMM_WORLD, 1 );

	}


	// 3. Lattice constant

	fprintf( fp1, "\tLattice constant, a [nm] = %13e\n", 0.1*parameters[0].a );
	fprintf( fp1, "\n" );


	// 4. Band structure parameters

	fprintf( fp1, "\tBand gap                                   [eV] = %13e\n", parameters[0].Eg  );
	fprintf( fp1, "\tSpin-orbit-splitting energy                [eV] = %13e\n", parameters[0].Dso );
	fprintf( fp1, "\tModified inverse electron effective mass   [--] = %13e\n", parameters[0].sc  );
	fprintf( fp1, "\tKane parameter                             [eV] = %13e\n", parameters[0].Ep  );
	fprintf( fp1, "\tFirst  modified Luttinger parameter        [--] = %13e\n", parameters[0].g1  );
	fprintf( fp1, "\tSecond modified Luttinger parameter        [--] = %13e\n", parameters[0].g2  );
	fprintf( fp1, "\tThird  modified Luttinger parameter        [--] = %13e\n", parameters[0].g3  );
	fprintf( fp1, "\n" );
	fprintf( fp1, "\tUnstrained spin-split-off band edge energy [eV] = %13e\n", parameters[0].E_so );
	fprintf( fp1, "\tUnstrained     light-hole band edge energy [eV] = %13e\n", parameters[0].E_lh );
	fprintf( fp1, "\tUnstrained     heavy-hole band edge energy [eV] = %13e\n", parameters[0].E_hh );
	fprintf( fp1, "\tUnstrained     conduction band edge energy [eV] = %13e\n", parameters[0].E_cb );
	fprintf( fp1, "\n" );
	fprintf( fp1, "\n" );


	//////////////////////////////////////////////////////////////////////////////////////
	//// Write out information for all additional shapes (materials) in the supercell ////
	//////////////////////////////////////////////////////////////////////////////////////

	int idx_shape;

	for( idx_shape = 1; idx_shape < N_shape; idx_shape++ )
	{

		// 1. Material system and alloy composition

		fprintf( fp1, "Shape %d: Embedded in shape %d\n", idx_shape, supercell[idx_shape].embedded );
		fprintf( fp1, "\n" );

		switch( supercell[idx_shape].material )
		{

			// Alloys of the form A(x)B(1-x)C(y)D(1-y)

			case 0: fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)P(%.4f)As(%.4f) \n", supercell[idx_shape].x, 1.0 - supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].y ); break;
			case 1: fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)P(%.4f)Sb(%.4f) \n", supercell[idx_shape].x, 1.0 - supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].y ); break;
			case 2: fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)As(%.4f)Sb(%.4f)\n", supercell[idx_shape].x, 1.0 - supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].y ); break;
			case 3: fprintf( fp1, "\tMaterial = Ga(%.4f)In(%.4f)P(%.4f)As(%.4f) \n", supercell[idx_shape].x, 1.0 - supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].y ); break;
			case 4: fprintf( fp1, "\tMaterial = Ga(%.4f)In(%.4f)P(%.4f)Sb(%.4f) \n", supercell[idx_shape].x, 1.0 - supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].y ); break;
			case 5: fprintf( fp1, "\tMaterial = Ga(%.4f)In(%.4f)As(%.4f)Sb(%.4f)\n", supercell[idx_shape].x, 1.0 - supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].y ); break;
			case 6: fprintf( fp1, "\tMaterial = Al(%.4f)In(%.4f)P(%.4f)As(%.4f) \n", supercell[idx_shape].x, 1.0 - supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].y ); break;
			case 7: fprintf( fp1, "\tMaterial = Al(%.4f)In(%.4f)P(%.4f)Sb(%.4f) \n", supercell[idx_shape].x, 1.0 - supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].y ); break;
			case 8: fprintf( fp1, "\tMaterial = Al(%.4f)In(%.4f)As(%.4f)Sb(%.4f)\n", supercell[idx_shape].x, 1.0 - supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].y ); break;


			// Alloys of the form A(x)B(y)C(1-x-y)D

			case 9 : fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)In(%.4f)P \n", supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].x - supercell[idx_shape].y ); break;
			case 10: fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)In(%.4f)As\n", supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].x - supercell[idx_shape].y ); break;
			case 11: fprintf( fp1, "\tMaterial = Al(%.4f)Ga(%.4f)In(%.4f)Sb\n", supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].x - supercell[idx_shape].y ); break;


			// Alloys of the form AB(x)C(y)D(1-x-y)

			case 12: fprintf( fp1, "\tMaterial = AlP(%.4f)As(%.4f)Sb(%.4f)\n", supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].x - supercell[idx_shape].y ); break;
			case 14: fprintf( fp1, "\tMaterial = InP(%.4f)As(%.4f)Sb(%.4f)\n", supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].x - supercell[idx_shape].y ); break;
			case 13: fprintf( fp1, "\tMaterial = GaP(%.4f)As(%.4f)Sb(%.4f)\n", supercell[idx_shape].x, supercell[idx_shape].y, 1.0 - supercell[idx_shape].x - supercell[idx_shape].y ); break;

		}

		fprintf( fp1, "\n" );


		// 2. Geometrical information - Shape, location, size, and rotation

		switch( supercell[idx_shape].shape )
		{

			case 0:

				fprintf( fp1, "\tShape = Cuboid\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z)    [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tSide length along x [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tSide length along y [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\tSide length along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L3 );

				break;


			case 1:

				fprintf( fp1, "\tShape = Sphere\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tRadius           [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );

				break;


			case 2:

				fprintf( fp1, "\tShape = Ellipsoid\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tSemi-axis length along x [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tSemi-axis length along y [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\tSemi-axis length along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L3 );

				break;


			case 3:

				fprintf( fp1, "\tShape = Hemisphere\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tRadius           [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );

				break;


			case 4:

				fprintf( fp1, "\tShape = Hemiellipsoid (lens)\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tSemi-axis length along x [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tSemi-axis length along y [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\tSemi-axis length along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L3 );

				break;


			case 5:

				fprintf( fp1, "\tShape = Rectangular-base pyramid\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z)    [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tBase half-length along x [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tBase half-length along y [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\tHeight           along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L3 );

				break;


			case 6:

				fprintf( fp1, "\tShape = Rectangular-base truncated pyramid\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z)    [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tBase half-length along x [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tBase half-length along y [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\tTop  half-length along x [nm] = %.4f\n", 0.1*supercell[idx_shape].L3 );
				fprintf( fp1, "\tHeight           along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L4 );

				break;


			case 7:

				fprintf( fp1, "\tShape = N-sided pyramid\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tDistance from centre of pyramid base to side = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\t                              Height along z = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\t                                No. of sides = %.4f\n", round( 0.1*supercell[idx_shape].L3 ) );

				break;


			case 8:

				fprintf( fp1, "\tShape = N-sided truncated pyramid\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tDistance from centre of truncated pyramid base to side = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\t Distance from centre of truncated pyramid top to side = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\t                                        Height along z = %.4f\n", 0.1*supercell[idx_shape].L3 );
				fprintf( fp1, "\t                                          No. of sides = %.4f\n", round( 0.1*supercell[idx_shape].L4 ) );

				break;


			case 9:

				fprintf( fp1, "\tShape = Cylinder\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tRadius           [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tHeight along z   [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );

				break;


			case 10:

				fprintf( fp1, "\tShape = Cylindrical ring\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tOuter radius in x-y plane [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tInner radius in x-y plane [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\t           Height along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L3 );

				break;


			case 11:

				fprintf( fp1, "\tShape = Cone\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tBase radius in x-y plane [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\t          Height along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );

				break;


			case 12:

				fprintf( fp1, "\tShape = Truncated cone\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tBase radius in x-y plane [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tTop  radius in x-y plane [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\t          Height along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L3 );

				break;


			case 13:

				fprintf( fp1, "\tShape = Truncated conic ring\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tOuter base radius in x-y plane [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tOuter top  radius in x-y plane [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\tInner      radius in x-y plane [nm] = %.4f\n", 0.1*supercell[idx_shape].L3 );
				fprintf( fp1, "\t                Height along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L4 );

				break;


			case 14:

				fprintf( fp1, "\tShape = Elliptical-base cone\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tSemi-axis length along x [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tSemi-axis length along y [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\t          Height along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L3 );

				break;


			case 15:

				fprintf( fp1, "\tShape = Truncated elliptical base cone\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tSemi-axis base length along x [nm] = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\tSemi-axis base length along y [nm] = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\tSemi-axis top  length along x [nm] = %.4f\n", 0.1*supercell[idx_shape].L3 );
				fprintf( fp1, "\t               Height along z [nm] = %.4f\n", 0.1*supercell[idx_shape].L4 );

				break;


			case 16:

				fprintf( fp1, "\tShape = N-sided prism\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				fprintf( fp1, "\tDistance from centre of prism base to side = %.4f\n", 0.1*supercell[idx_shape].L1 );
				fprintf( fp1, "\t                            Height along z = %.4f\n", 0.1*supercell[idx_shape].L2 );
				fprintf( fp1, "\t                              No. of sides = %.4f\n", round( 0.1*supercell[idx_shape].L3 ) );

				break;


			case 17:

				fprintf( fp1, "\tShape = Quantum well\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				if ( round( 0.1*supercell[idx_shape].L1 ) == 0 )
				{
					
					fprintf( fp1, "\tQuantum confinement direction = [100]\n" );

				}
				else if ( round( 0.1*supercell[idx_shape].L1 ) == 1 )
				{
					
					fprintf( fp1, "\tQuantum confinement direction = [010]\n" );

				}
				else if ( round( 0.1*supercell[idx_shape].L1 ) == 2 )
				{
					
					fprintf( fp1, "\tQuantum confinement direction = [001]\n" );

				}
				else
				{
					
					fprintf( fp1, "\tQuantum confinement direction = undefined (check supercell.in)\n" );

				}
				fprintf( fp1, "\tLength = %.4f\n", 0.1*supercell[idx_shape].L2 );

				break;


			case 18:

				fprintf( fp1, "\tShape = Quantum wire\n" );

				fprintf( fp1, "\n" );
				fprintf( fp1, "\tCentre (x, y, z) [nm] = (%.4f, %.4f, %.4f)\n", 0.1*supercell[idx_shape].Ox, 0.1*supercell[idx_shape].Oy, 0.1*supercell[idx_shape].Oz );
				if ( round( 0.1*supercell[idx_shape].L1 ) == 0 )
				{
					
					fprintf( fp1, "\tQuantum confinement plane = (100)\n" );

				}
				else if ( round( 0.1*supercell[idx_shape].L1 ) == 1 )
				{
					
					fprintf( fp1, "\tQuantum confinement plane = (010)\n" );

				}
				else if ( round( 0.1*supercell[idx_shape].L1 ) == 2 )
				{
					
					fprintf( fp1, "\tQuantum confinement plane = (001)\n" );

				}
				else
				{
					
					fprintf( fp1, "\tQuantum confinement plane = undefined (check supercell.in)\n" );

				}
				
				if ( round( 0.1*supercell[idx_shape].L1 ) == 0 )
				{
					
					fprintf( fp1, "\tLength along [010] = %.4f\n", 0.1*supercell[idx_shape].L2 );
					fprintf( fp1, "\tLength along [001] = %.4f\n", 0.1*supercell[idx_shape].L3 );

				}
				else if ( round( 0.1*supercell[idx_shape].L1 ) == 1 )
				{
					
					fprintf( fp1, "\tLength along [100] = %.4f\n", 0.1*supercell[idx_shape].L2 );
					fprintf( fp1, "\tLength along [001] = %.4f\n", 0.1*supercell[idx_shape].L3 );

				}
				else if ( round( 0.1*supercell[idx_shape].L1 ) == 2 )
				{
					
					fprintf( fp1, "\tLength along [100] = %.4f\n", 0.1*supercell[idx_shape].L2 );
					fprintf( fp1, "\tLength along [010] = %.4f\n", 0.1*supercell[idx_shape].L3 );

				}

				break;


		}

		fprintf( fp1, "\tShape alpha rotation angle (degrees) = %.4f\n", supercell[idx_shape].alpha 	);
		fprintf( fp1, "\tShape beta rotation angle (degrees)  = %.4f\n", supercell[idx_shape].beta 		);
		fprintf( fp1, "\tShape gamma rotation angle (degrees) = %.4f\n", supercell[idx_shape].gamma 	);
		fprintf( fp1, "\n" );


		// 3. Lattice constant

		fprintf( fp1, "\tLattice constant, a [nm] = %13e\n", 0.1*parameters[idx_shape].a );
		fprintf( fp1, "\n" );


		// 4. Band structure parameters

		fprintf( fp1, "\tBand gap                                   [eV] = %13e\n", parameters[idx_shape].Eg  );
		fprintf( fp1, "\tSpin-orbit-splitting energy                [eV] = %13e\n", parameters[idx_shape].Dso );
		fprintf( fp1, "\tModified inverse electron effective mass   [--] = %13e\n", parameters[idx_shape].sc  );
		fprintf( fp1, "\tKane parameter                             [eV] = %13e\n", parameters[idx_shape].Ep  );
		fprintf( fp1, "\tFirst  modified Luttinger parameter        [--] = %13e\n", parameters[idx_shape].g1  );
		fprintf( fp1, "\tSecond modified Luttinger parameter        [--] = %13e\n", parameters[idx_shape].g2  );
		fprintf( fp1, "\tThird  modified Luttinger parameter        [--] = %13e\n", parameters[idx_shape].g3  );
		fprintf( fp1, "\n" );
		fprintf( fp1, "\tUnstrained spin-split-off band edge energy [eV] = %13e\n", parameters[idx_shape].E_so );
		fprintf( fp1, "\tUnstrained     light-hole band edge energy [eV] = %13e\n", parameters[idx_shape].E_lh );
		fprintf( fp1, "\tUnstrained     heavy-hole band edge energy [eV] = %13e\n", parameters[idx_shape].E_hh );
		fprintf( fp1, "\tUnstrained     conduction band edge energy [eV] = %13e\n", parameters[idx_shape].E_cb );
		fprintf( fp1, "\n" );
		fprintf( fp1, "\n" );


	}

	fclose( fp1 );

} // End of function
