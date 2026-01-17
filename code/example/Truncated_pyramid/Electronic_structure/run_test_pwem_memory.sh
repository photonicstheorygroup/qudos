#!/bin/bash

# QUantum Dot Open-source Simulator (QUDOS)
#
# Copyright (c) 2025
# Tommy Murphy and Christopher A. Broderick
#
# This file is part of QUDOS.
#
# QUDOS is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your discretion) any later version.
#
# QUDOS is distributed in the hope that it will be useful, but without any warranty and
# without an implied warranty of merchantability or fitness for a particular purpose.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with QUDOS.
# If not, see <https://www.gnu.org/licenses/>.
#
# Alternative commercial licensing terms may be available from the copyright holders upon request.
#
#
#
#
# Test run script for local installation for the electronic structure memory calculator code
#


export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

mpirun -np 4 ../../../build/8_band_pwem_memory.ex