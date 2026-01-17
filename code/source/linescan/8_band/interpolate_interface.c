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

void interpolate_interface( int N_shape, double *shape_frac, struct supercell_geometry *supercell, struct material_params *parameters, struct material_params *interface_params )
{

	interface_params->a 	= parameters[0].a;
	interface_params->C11 	= parameters[0].C11;
	interface_params->C12 	= parameters[0].C12;
	interface_params->C44 	= parameters[0].C44;
	interface_params->e14 	= parameters[0].e14;
	interface_params->B114	= parameters[0].B114;
	interface_params->B124	= parameters[0].B124;
	interface_params->B156 	= parameters[0].B156;
	interface_params->Eg 	= parameters[0].Eg;
	interface_params->Dso 	= parameters[0].Dso;
	interface_params->Ep 	= parameters[0].Ep;
	interface_params->P 	= parameters[0].P;
	interface_params->sc 	= parameters[0].sc;
	interface_params->g1 	= parameters[0].g1;
	interface_params->g2 	= parameters[0].g2;
	interface_params->g3 	= parameters[0].g3;
	interface_params->ac 	= parameters[0].ac;
	interface_params->av 	= parameters[0].av;
	interface_params->b 	= parameters[0].b;
	interface_params->d 	= parameters[0].d;
	interface_params->Ev_avg = parameters[0].Ev_avg;
	interface_params->E_so 	= parameters[0].E_so;
	interface_params->E_lh 	= parameters[0].E_lh;
	interface_params->E_hh 	= parameters[0].E_hh;
	interface_params->E_cb 	= parameters[0].E_cb;
	interface_params->nr 	= parameters[0].nr;

	for( int idx_shape = 1; idx_shape < N_shape; idx_shape++ )
	{

		interface_params->a 		+= ( shape_frac[idx_shape]*parameters[idx_shape].a - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].a );
		interface_params->C11		+= ( shape_frac[idx_shape]*parameters[idx_shape].C11 - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].C11 );
		interface_params->C12 		+= ( shape_frac[idx_shape]*parameters[idx_shape].C12 - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].C12 );
		interface_params->C44 		+= ( shape_frac[idx_shape]*parameters[idx_shape].C44 - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].C44 );
		interface_params->e14 		+= ( shape_frac[idx_shape]*parameters[idx_shape].e14 - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].e14 );
		interface_params->B114 		+= ( shape_frac[idx_shape]*parameters[idx_shape].B114 - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].B114 );
		interface_params->B124 		+= ( shape_frac[idx_shape]*parameters[idx_shape].B124 - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].B124 );
		interface_params->B156 		+= ( shape_frac[idx_shape]*parameters[idx_shape].B156 - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].B156 );
		interface_params->Eg 		+= ( shape_frac[idx_shape]*parameters[idx_shape].Eg - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].Eg );
		interface_params->Dso 		+= ( shape_frac[idx_shape]*parameters[idx_shape].Dso - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].Dso );
		interface_params->Ep 		+= ( shape_frac[idx_shape]*parameters[idx_shape].Ep - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].Ep );
		interface_params->P 		+= ( shape_frac[idx_shape]*parameters[idx_shape].P - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].P );
		interface_params->sc 		+= ( shape_frac[idx_shape]*parameters[idx_shape].sc - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].sc );
		interface_params->g1 		+= ( shape_frac[idx_shape]*parameters[idx_shape].g1 - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].g1 );
		interface_params->g2 		+= ( shape_frac[idx_shape]*parameters[idx_shape].g2 - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].g2 );
		interface_params->g3 		+= ( shape_frac[idx_shape]*parameters[idx_shape].g3 - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].g3 );
		interface_params->ac 		+= ( shape_frac[idx_shape]*parameters[idx_shape].ac - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].ac );
		interface_params->av 		+= ( shape_frac[idx_shape]*parameters[idx_shape].av - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].av );
		interface_params->b 		+= ( shape_frac[idx_shape]*parameters[idx_shape].b - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].b );
		interface_params->d 		+= ( shape_frac[idx_shape]*parameters[idx_shape].d - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].d );
		interface_params->Ev_avg	+= ( shape_frac[idx_shape]*parameters[idx_shape].Ev_avg - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].Ev_avg );
		interface_params->E_so 		+= ( shape_frac[idx_shape]*parameters[idx_shape].E_so - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].E_so );
		interface_params->E_lh 		+= ( shape_frac[idx_shape]*parameters[idx_shape].E_lh - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].E_lh );
		interface_params->E_hh 		+= ( shape_frac[idx_shape]*parameters[idx_shape].E_hh - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].E_hh );
		interface_params->E_cb 		+= ( shape_frac[idx_shape]*parameters[idx_shape].E_cb - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].E_cb );
		interface_params->nr 		+= ( shape_frac[idx_shape]*parameters[idx_shape].nr - shape_frac[idx_shape] *parameters[supercell[idx_shape].embedded].nr );

	}

}
