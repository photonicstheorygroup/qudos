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
//
//
//
// Last modified: 27/05/2015
//
//
// III-V material parameters for plane wave quantum dot calculations
//
// The material parameters contained in this header file are used in the material_parameters.c subroutines to obtain
// (by interpolation) those for the user-specified quaternary alloys contained in the calculational supercell
//
//
// 	Unless otherwise specified, these material parameters have been taken from:
//
// 	[1] I. Vurgaftman, J. R. Meyer and L. R. Ram-Mohan, J. Appl. Phys. 89, 5815 (2001)
//
//	with the exception of the piezoelectric coefficients, which are taken from:
//
//	[2] M. A. Caro, S. Schulz and E. P. O'Reilly, Phys. Rev. B 91, 075203 (2015)
//
//   while parameters for III-Bi compounds are taken from
//
//	[3] M. Ferhat and A. Zaoui, Phys. Rev. B 73, 115107 (2006)
//	[4] A. Janotti, S.-H. Wei and S. B. Zhang, Phys. Rev. B 65, 115203 (2002)
//	[5] S. Q. Wang and H. Q. Ye, Phys. Rev. B 66, 235111 (2002)
//	[6] P. Carrier and S.-H. Wei, Phys. Rev. B 70, 035212 (2004)
//	[7] S. Q. Wang and H. Q. Ye, Phys. Stat. Sol. B 240, 45 (2003)
//
//
// 	Indices allocated for the storage of material parameters in the arrays contained this header file:
//
// 		Binaries are stored as:
//
// 			0	AlN
// 			1	AlP
// 			2	AlAs
// 			3	AlSb
// 			4	AlBi
//
// 			5	GaN
// 			6	GaP
// 			7	GaAs
// 			8	GaSb
// 			9	GaBi
//
// 			10	InN
// 			11	InP
// 			12	InAs
// 			13	InSb
// 			14	InBi
//
// 		Note: N/Bi containing compounds are typically with in the codes using band anti-crossing models. For this reason,
// 			 we do not need the band gaps, etc. of these binaries since they are produced at the end of calculations by
//			 the modified k.p Hamiltonian in use. However, we do need the values of the lattice and elastic constants,
//			 as well as the deformation potentials of these binaries, so that the effects of strain in the system due to
//			 the incorporation N and/or Bi can be determined.
//
//
// 		Ternaries are stored as:
//
// 			0	AlGaP
// 			1	AlGaAs
// 			2	AlGaSb
//
// 			3	GaInP
// 			4	GaInAs
// 			5	GaInSb
//
// 			6	AlInP
// 			7	AlInAs
//			8 	AlInSb
//
//			9	AlPAs
//			10	AlPSb
//			11	AlAsSb
//
//			12	GaPAs
//			13	GaPSb
//			14	GaAsSb
//
//			15	InPAs
//			16	InPSb
//			17	InAsSb
//
//
// 		Quaternaries are stored as:
//
// 			0	AlGaPAs
// 			1	AlGaPSb
// 			2	AlGaAsSb
//
// 			3	GaInPAs
// 			4	GaInPSb
// 			5	GaInAsSb
//
// 			6	AlInPAs
// 			7	AlInPSb
//			8 	AlInAsSb
//
//			9	AlGaInP
//			10	AlGaInAs
//			11	AlGaInSb
//
//			12	AlPAsSb
//			13	GaPAsSb
//			14	InPAsSb
//
//
// Units:
//		In order to remain consistent with the conventions in the codes, all energies should be in eV and all lengths
//		in Angstrom. Effective masses should be in units of the free electron mass.
//
//		CHECK THIS FOR QUANTUM DOT CASE AND FIX UNITS -----> It should also be noted that the strain dependent parts of all k.p Hamiltonians rely only on ratios of elastic constants,
//		and for this reason it is immaterial which units are used when specifying them (just be careful to ensure that they
//		are all in the same units).


//////////////////////////////////////////////////////////////////////
//// 1. Elastic and piezoelectric parameters for binary compounds ////
//////////////////////////////////////////////////////////////////////

// Note: Due to a lack of information regarding the piezoelectric coefficients for for III-Bi compounds,
//       we have set them equal to the corresponding III-Sb values for each III-Bi compound

double    a[15] = { 4.3800, 5.4672, 5.6611, 6.1355, 6.3417,
				4.5000, 5.4505, 5.6503, 6.0959, 6.3240,
				4.9800, 5.8697, 6.0553, 6.4794, 6.6860 }; // Lattice constant [Angstrom]

double  C11[15] = { 304.00, 133.00, 125.00,  87.69,  72.20,
				293.00, 140.50, 118.8,  88.42,  73.00,
				187.00, 101.10,  83.3,  68.47,  58.50 }; // Elastic constant [GPa]

double  C12[15] = { 160.00,  63.00,  53.40,  43.41,  36.20,
				159.00,  62.03,  53.8,  40.26,  32.70,
				125.00,  56.10,  45.3,  37.35,  30.40 }; // Elastic constant [GPa]

double  C44[15] = { 193.00,  61.50,  54.20,  40.76,  33.40,
				155.00,  70.33,  59.4,  43.22,  36.30,
				 86.00,  45.60,  39.6,  31.11,  27.50 }; // Elastic constant [GPa]

double  e14[15] = {  0.548,  0.014, -0.055, -0.094, -0.094,
				 0.366, -0.121, -0.160, -0.216, -0.216,   // 
				 0.593,  0.016, -0.045, -0.161, -0.161 }; // First-order piezoelectric coefficient [C/m^2]

double  e14_calc[15] = {  0.548,  0.014, -0.055, -0.094, -0.094,
				 0.366, -0.121, -0.230, -0.216, -0.216,   // 
				 0.593,  0.016, -0.115, -0.161, -0.161 }; // First-order piezoelectric coefficient [C/m^2]

double B114[15] = { -6.810, -2.020, -1.610, -0.760, -0.760,
				-5.380, -1.230, -0.439, -0.310, -0.310,   // 
				-5.960, -1.540, -0.531, -0.620, -0.620 }; // Second-order piezoelectric coefficient [C/m^2]

double B124[15] = { -5.040, -2.760, -2.590, -1.990, -1.990,
				-6.730, -3.270, -3.765, -2.770, -2.770,   // 
				-6.320, -3.620, -4.076, -4.040, -4.040 }; // Second-order piezoelectric coefficient [C/m^2]

double B156[15] = { -4.150, -1.430, -1.320, -0.820, -0.820,
				-3.180, -1.380, -0.492, -0.700, -0.700,   // 
				-2.000, -1.020, -0.120, -0.160, -0.160 }; // Second-order piezoelectric coefficient [C/m^2]

/*

// TEST PARAMETERS FROM SCHLIWA Phys. Rev. B - REMOVE AND UNCOMMENT ORIGINALS ABOVE WHEN FINISHED BENCHMARKING

double    a[15] = { 4.3800, 5.4672, 5.6611, 6.1355, 6.3417,
				4.5000, 5.4505, 5.6503, 6.0959, 6.3240,
				4.9800, 5.8697, 6.0553, 6.4794, 6.6860 }; // Lattice constant [Angstrom]

double  C11[15] = { 304.00, 133.00, 125.00,  87.69,  72.20,
				293.00, 140.50, 118.80,  88.42,  73.00,
				187.00, 101.10,  83.30,  68.47,  58.50 }; // Elastic constant [GPa]

double  C12[15] = { 160.00,  63.00,  53.40,  43.41,  36.20,
				159.00,  62.03,  53.80,  40.26,  32.70,
				125.00,  56.10,  45.30,  37.35,  30.40 }; // Elastic constant [GPa]

double  C44[15] = { 193.00,  61.50,  54.20,  40.76,  33.40,
				155.00,  70.33,  59.40,  43.22,  36.30,
				 86.00,  45.60,  39.60,  31.11,  27.50 }; // Elastic constant [GPa]

double  e14[15] = {  0.548,  0.014, -0.055, -0.094, -0.094,
				 0.366, -0.121, -0.230, -0.216, -0.216,
				 0.593,  0.016, -0.115, -0.161, -0.161 }; // First-order piezoelectric coefficient [C/m^2]

double B114[15] = { -6.810, -2.020, -1.610, -0.760, -0.760,
				-5.380, -1.230, -0.439, -0.310, -0.310,
				-5.960, -1.540, -0.531, -0.620, -0.620 }; // Second-order piezoelectric coefficient [C/m^2]

double B124[15] = { -5.040, -2.760, -2.590, -1.990, -1.990,
				-6.730, -3.270, -3.765, -2.770, -2.770,
				-6.320, -3.620, -4.076, -4.040, -4.040 }; // Second-order piezoelectric coefficient [C/m^2]

double B156[15] = { -4.150, -1.430, -1.320, -0.820, -0.820,
				-3.180, -1.380, -0.492, -0.700, -0.700,
				-2.000, -1.020, -0.120, -0.160, -0.160 }; // Second-order piezoelectric coefficient [C/m^2]

// END OF TEST PARAMETERS

*/


////////////////////////////////////////////////////////
//// 3. Deformation potentials for binary compounds ////
////////////////////////////////////////////////////////

// Note: Deformation potentials follow the sign conventions of Krijn - Semicond. Sci. Technol. 6, 27 (1991) - so that all
//       deformation potentials are negative, except for the valence band hydrostatic deformation potential
//
//       Due to a lack of information regarding the deformation potentials for for III-Bi compounds,
//       we have set them equal to the corresponding III-Sb values for each III-Bi compound

double   ac[15] = { -6.000, -5.700, -5.640, -4.500, -4.500,
				-2.200, -8.200, -8.013, -7.500, -7.500,
				-1.850, -6.000, -5.080, -6.940, -6.940 }; // Conduction band deformation potential [eV]

double   av[15] = {  3.400,  3.000,  2.470,  1.400,  1.400,
				 5.200,  1.700,  0.220,  0.800,  0.800,
				 1.500,  0.600,  1.000,  0.360,  0.360 }; // Valence band deformation potential [eV]

double    b[15] = { -1.900, -1.500, -2.300, -1.350, -1.350,
				-2.200, -1.600, -1.824, -2.000, -2.000,
				-1.200, -2.000, -1.800, -2.000, -2.000 }; // Axial deformation potential [eV]

double    d[15] = { -10.00, -4.600, -3.400, -4.300, -4.300,
				-3.400, -4.600, -5.062, -4.700, -4.700,
				-9.300, -5.000, -3.600, -4.700, -4.700 }; // Shear deformation potential [eV]


///////////////////////////////////////////////////////////
//// 3. Band structure parameters for binary compounds ////
///////////////////////////////////////////////////////////

// Note: For N- and Bi-containing alloys we need only calculate the temperature dependence of the host matrix band gap, so the Varshni parameters are not used for III-N/Bi compounds
//
//       All that matters for the average energy of the top three valence bands are the differences in this value between different materials (i.e. the valence band offset). For that
//       reason the average valence band energies listed here are simple the valence band offsets provided in [1], adjusted so that the GaAs VB edge is taken as the zero of energy.
//       The valence band offsets for III-Bi compounds were obtained by extrapolating linear fits of valence band offset as a function of lattice constant for N-, P-, As- and Sb-
//       containing compounds (having the same species of cation).

double Ev_avg[15] = { -2.640, -0.940, -0.530, 0.3900, 1.1000,
				  -1.840, -0.470,  0.000, 0.7700, 1.1000,
				  -1.580, -0.140,  0.173, 0.8000, 1.1000 }; // Average valence band energy [eV]

double     Eg[15] = { 4.9000, 3.6300, 3.0990, 2.3860,  0.042,           
				  3.2990, 2.8860, 1.5180, 0.8120, -1.450,
				  1.9400, 1.4240, 0.4050, 0.2350, -1.630 }; // Band gap energy at T = 0 K [eV]

double    Dso[15] = { 0.0190, 0.0700, 0.2800, 0.6760, 2.1240,
				  0.0170, 0.0800, 0.3400, 0.7600, 2.1500,
				  0.0060, 0.1080, 0.3800, 0.8100, 2.1500 }; // Spin-orbit-splitting energy [eV]

double     mc[15] = { 0.2500, 0.2200, 0.1500, 0.1400, 99.999,
				  0.1500, 0.1300, 0.0667, 0.0390, 99.999,
				  0.1200, 0.0795, 0.02226, 0.0135, 99.999 }; // Electron effective mass [free electron mass]

double     Ep[15] = { 27.100, 17.700, 21.100, 18.700, 99.999,
				  25.000, 31.400, 25.7, 27.000, 99.999,
				  25.000, 20.700, 21.846, 23.300, 99.999 }; // Kane parameter [eV]

double    gL1[15] = { 1.9200, 3.3500, 3.7600, 5.1800, 99.999,
				  2.6700, 4.0500, 7.1, 13.400, 99.999,
				  3.7200, 5.0800, 19.67, 34.800, 99.999 }; // First Luttinger parameter [dimensionless]

double    gL2[15] = { 0.4700, 0.7100, 0.8200, 1.1900, 99.999,
				  0.7500, 0.4900, 2.0200, 4.7000, 99.999,
				  1.2600, 1.6000, 8.4000, 15.000, 99.999 }; // Second Luttinger parameter [dimensionless]

double    gL3[15] = { 0.8500, 1.2300, 1.4200, 1.9700, 99.999,
				  1.1000, 2.9300, 2.9100, 6.0000, 99.999,
				  1.6300, 2.1000, 9.3000, 16.500, 99.999 }; // Third Luttinger parameter [dimensionless]

double     nr[15] = { 2.1600, 2.7500, 3.0000, 3.1900, 3.1900,
				  2.4000, 2.9000, 3.6304, 3.7900, 3.7900, // 3.3000
				  2.5900, 3.1000, 3.8210, 3.9500, 3.9500 }; // Refractive index [dimensionless]


//////////////////////////////////////////////////////////////////////////////////////
//// 4. Band gap temperature dependence for binary compounds - Varshni parameters ////
//////////////////////////////////////////////////////////////////////////////////////

// Note: For GaP we use the functional form of the band gap temperature dependence given in [1] ==> Varshni parameters not used for GaP, so have set them = 0
//
//	    Due to a lack of information regarding the deformation potentials and piezoelectric coefficients for for III-Bi compounds, we have set them equal to the
//       corresponding III-Sb values for each III-Bi compound
//
//       For N- and Bi-containing alloys we need only calculate the temperature dependence of the host matrix band gap, so the Varshni parameters are not used
//	    in practice for III-N/Bi compounds

double alpha_T[15] = { 0.5930, 0.5771, 0.8850, 0.4200, 0.4200,
				   0.5930, 0.0000, 0.5405, 0.4170, 0.4170,
				   0.2450, 0.3630, 0.2760, 0.4200, 0.4200 }; // First Varshni parameter [meV/K]

double  beta_T[15] = { 600.00, 372.00, 530.00, 140.00, 140.00,
				   600.00, 0.0000, 204.00, 140.00, 140.00,
				   624.00, 162.00,  93.00, 170.00, 170.00 }; // Second Varshni parameter [K]


/////////////////////////////////////////////////
//// 5. Bowing parameters for ternary alloys ////
/////////////////////////////////////////////////

// Note: For AlGaAs and AlGaSb, [1] gives the composition dependent bowing parameters -0.127 + 1.31*x and -0.044 + 1.22*x, respectively - we have taken the average values of
//       these bowing parameters by setting x = 0.5
//
//       [1] does not provide bowing parameters for the spin-orbit-splitting energies in AlGaP, AlPAs, AlPSb, GaPAs and GaPSb, so have set them = 0

double  bow_Eg_tern[18] = { 0.000, 0.528, 0.566,
					   0.650, 0.477, 0.415,
					   -0.48, 0.700, 0.430,
					   0.220, 2.700, 0.800,
					   0.190, 2.700, 1.430,
					   0.100, 1.900, 0.670 }; // Bowing coefficient for band gap energy [eV]

double bow_Dso_tern[18] = { 0.000, 0.000, 0.300,
					   0.000, 0.150, 0.100,
					   -0.19, 0.150, 0.250,
					   0.000, 0.000, 0.150,
					   0.000, 0.000, 0.600,
					   0.160, 0.750, 1.200 }; // Bowing coefficient for spin-orbit-splitting energy [eV]


////////////////////////////////////////////////////
//// 6. Bowing parameters for quaternary alloys ////
////////////////////////////////////////////////////

double  bow_Eg_quat[15] = { 0.0000, 0.0000, 0.0000,
					   0.0000, 0.0000, 0.0000,
					   0.0000, 0.0000, 0.0000,
					   0.0000, 0.0000, 0.0000,
					   0.0000, 0.0000, 0.0000 };

double bow_Dso_quat[15] = { 0.0000, 0.0000, 0.0000,
					   0.0000, 0.0000, 0.0000,
					   0.0000, 0.0000, 0.0000,
					   0.0000, 0.0000, 0.0000,
					   0.0000, 0.0000, 0.0000 };
