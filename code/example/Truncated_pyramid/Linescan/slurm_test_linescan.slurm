#!/bin/bash -l

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
# Example SLURM job run script for linescan code that requests 1 node (tested on Meluxina)
#


#SBATCH --job-name=linescan
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --qos=default
#SBATCH -t 00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --ntasks-per-core=1


# Loading modules

module load env/release/2023.1
module load intel-compilers/2023.1.0
module load imkl/2023.1.0
module load impi/2021.9.0-intel-compilers-2023.1.0
module load GSL/2.7-GCC-12.3.0


# Setting 128 OpenMP threads per process

export OMP_NUM_THREADS=128
export MKL_NUM_THREADS=$OMP_NUM_THREADS


# SLURM boilerplate

cd $SLURM_SUBMIT_DIR # job submit directory


# Running calculation
	
srun "code folder root"/code/build/8_band_line.ex
