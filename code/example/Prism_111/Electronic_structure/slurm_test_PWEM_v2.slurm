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
# Example SLURM job run script that requests 1 node with 32 tasks per node and 2 cpus per task (tested on UCC cluster)
#


#SBATCH --job-name=8_band
#SBATCH -p cpuq-long
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --threads-per-core=1
#SBATCH --ntasks-per-core=1


# Loading modules

module load slurm oneapi oneapi-hpc spack
module load gsl libyaml


# Explicitly setting path variables for Intel MPI (may not be necessary)

export FI_PROVIDER_PATH=/cm/shared/apps/oneapi-hpc/2025.2.1.44/mpi/2021.16/opt/mpi/libfabric/lib/prov
export LD_LIBRARY_PATH=/cm/shared/apps/oneapi-hpc/2025.2.1.44/mpi/2021.16/opt/mpi/libfabric/lib:$LD_LIBRARY_PATH


# Setting 2 OpenMP threads per process

export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=$OMP_NUM_THREADS


# SLURM boilerplate

cd $SLURM_SUBMIT_DIR # job submit directory


# Run calculation

mpirun -n 32 "QUDOS root"/code/build/8_band_pwem.ex


# Modifying number of plane waves

#mpirun -n 32 "QUDOS root"/code/build/8_band_pwem.ex -M_Gx=7 -M_Gy=7 -M_Gz=10 -M_S_Gx=21 -M_S_Gy=21 -M_S_Gz=40 -M_C_Gx=7 -M_C_Gy=7 -M_C_Gz=20




