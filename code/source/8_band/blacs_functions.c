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

int indxl2g(int iloc, int nb, int iproc, int isrcproc, int nprocs)
{

	return (nprocs * nb * ((iloc) / nb) + ((iloc) % nb) + ((nprocs + iproc - isrcproc) % nprocs) * nb);

}

int indxg2l(int iglob, int nb, int proc, int isrcproc, int nprocs)
{

	return (nb * ((iglob) / (nb * nprocs)) + ((iglob) % nb));

}
