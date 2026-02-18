.. QUDOS documentation master file, created by
   sphinx-quickstart on Mon Feb  9 15:40:01 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
============

Installation of QUDOS requires compilation of the C source code using an MPI-supported compiler. 
It is recommended to use the MPI-wrapped Intel C compiler (mpiicx) that is provided as part of the Intel oneAPI HPC Toolkit, as the code employs routines from the Intel oneAPI Math Kernel Library (oneMKL).
The code has been tested to compile and run correctly on the following systems:

.. list-table:: Systems tested with QUDOS
   :header-rows: 1
   :stub-columns: 1
   :widths: 18 20 20 20 20

   * -  
     - System 1
     - System 2
     - System 3
     - System 4
   * - OS
     - Ubuntu 20.04
     - CentOS 7
     - Rocky 8.10
     - RHEL 9.6
   * - CPU hardware
     - Intel Core i7-10700K
     - Intel Xeon Platinum 8260
     - AMD EPYC 7H12
     - Intel Xeon Gold 6430
   * - Compiler version
     - Intel mpiicx 2021.2.0
     - Intel mpiicc 2019.5
     - Intel mpiicx 2023.1
     - Intel mpiicx 2025.2.1
   * - MPI implementation
     - Intel MPI 2021.2.0
     - Intel MPI 2019.5
     - Intel MPI 2021.9
     - Intel MPI 2021.16
   * - Scheduler/launcher
     - Local; mpirun (Hydra)
     - PBS; mpirun
     - Slurm; srun
     - Slurm; mpirun




Prerequisites
~~~~~~~~~~~~~

1. **LibYAML**

  LibYAML is used to read input files. To install, follow the instructions in `the Github repo <https://github.com/yaml/libyaml>`__.

2. **GNU Scientific Library (GSL)**

  `GSL <https://www.gnu.org/software/gsl/>`__ is used for numerical integration. To install, follow the instructions in the included README and INSTALL files in the GSL release. 

3. **Intel oneAPI Math Kernal Library (oneMKL)**

  Functions from oneMKL are used to compute fast-Fourier transforms, and oneMKL also provides optimised Basic Linear Algebra Subprograms (BLAS) and ScaLAPACK routines. 
  oneMKL is available as part of the `Intel oneAPI Base Toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html>`__. 

  Note that it is recommended to use version 2025.1 or later as a bug which can cause the employed ScaLAPACK eigensolver to hang was fixed in the 2025.1 release.
  If such a version is unavailable, versions 2023.1 and earlier are also likely to be used without issue.

4. **Intel MPI Library (recommended)**

  It is recommended to compile the code using Intel's mpiicx C compiler wrapper, which is provided as part of the Intel MPI Library. The Intel MPI library is available as part of the `Intel oneAPI HPC Toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html>`__.
  QUDOS provides example compilation scripts that employ the mpiicx compiler.


.. _Compilation:

Compilation instructions
~~~~~~~~~~~~~~~~~~~~~~~~

1.  Navigate to qudos/code/build and create a symbolic link to the makefile to be employed via::

      ln -s makefile_icx makefile

    where makefile_icx is an example makefile provided in the build folder and compiles the code using the mpiicx compiler wrapper.
    An example makefile is also included which uses the mpiicc compiler wrapper to support older Intel MKL versions (makefile_icc).
    Furthermore, an example makefile is also included for cluster systems that make use of the Spack package manager and that may require specific linking to LibYAML and GSL libraries (makefile_icx_spack).

2.  Prepare a script for your cluster's job scheduler to build QUDOS using the makefile via::

      make

    It is recommended to compile QUDOS on the back end of the cluster to enable QUDOS to be as optimised as possible.
    Example job scripts are provided for both PBS and Slurm.
    In your job script, ensure the environment is set up correctly by loading modules relating to LibYAML, GSL, and oneAPI (and oneAPI HPC if required).
    If compiling on a local system, ensure the Intel oneAPI environment is set up correctly by sourcing the setvars.sh script as described `here <https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/2025-1/use-the-setvars-and-oneapi-vars-scripts-with-linux.html>`__ before running the 'make' command.


Test calculations
~~~~~~~~~~~~~~~~~

Once QUDOS has compiled successfully, perform the following two tests:

1.  **Linescan code** 

    Navigate to qudos/code/example/Truncated_pyramid/Linescan and prepare the job script. Example PBS, Slurm, and local (.sh) job scripts have been provided.
    The environment should be set up in the same manner as when compiling the code. To execute the linescan code, include the following in the job script::

      "QUDOS root"/code/build/8_band_line.ex

    where "QUDOS root" is replaced with the path to the QUDOS root directory.
    As the linescan code is parallelised via OpenMP, ensure the number of OpenMP threads is set to the number of cores requested via the job scheduler by setting the OMP_NUM_THREADS and MKL_NUM_THREADS environment variables.

    After code execution, the output files will be in the folder "Output_files...". Verify that the output files are the same as those in the folder "Linescan_test_ref" to ensure that the linescan code is working correctly.

2.  **Electronic structure code**

    Navigate to qudos/code/example/Truncated_pyramid/Electronic_structure and prepare the job script. Example PBS, Slurm, and local (.sh) job scripts have been provided.

    The environment should be set up in the same manner as when compiling the code. To execute the linescan code, include the following in the job script::

      mpirun -n N "QUDOS root"/code/build/8_band_pwem.ex

    where N is the number of MPI processes and "QUDOS root" is replaced with the path to the QUDOS root directory. As the linescan code is parallelised via MPI and OpenMP, ensure the number of MPI processes times the number of OpenMP threads is equal to the number of cores requested via the job scheduler, and there is at least one MPI process per node. The number of OpenMP threads is set via the OMP_NUM_THREADS and MKL_NUM_THREADS environment variables.

    After code execution, the output files will be in the folder "Output_files...". Verify that the output files are the same as those in the folder "PWEM_test_ref" to ensure that the electronic structure code is working correctly.
    Note that the individual complex values of the output momentum and Coulomb matrix elements will differ from run-to-run, and thus from those in the "PWEM_test_ref" folder, due to how the eigenvectors for degenerate energy states are output from the eigensolver.
    However, in the case of the momentum matrix elements, the sum of the squared magnitudes of the four individual values describing optical transitions between pairs of spin-degenerate states will always be constant. Similarly, configuration interaction method calculations performed using the Coulomb matrix elements obtained via distinct runs of the eigensolver will always produce identical excitonic energies.