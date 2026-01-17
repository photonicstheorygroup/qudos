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
# Example SLURM job run script to compile QUDOS code (tested on UCC cluster)
#


#SBATCH --job-name=8_band
#SBATCH -p cpuq-long
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --ntasks-per-core=1


# Loading modules

module load slurm oneapi oneapi-hpc spack
module load gsl libyaml


# SLURM boilerplate

cd $SLURM_SUBMIT_DIR # job submit directory


# Compiling QUDOS

make
