
QUDOS
=====

.. start-docs

QUDOS is the QUantum Dot Open-source Simulator. QUDOS is an free, open-source software package that underpins predictive analysis and design of structures incorporating quantum dots for materials with zinc-blende crystal structures.

QUDOS implements a semi-analytical reciprocal space plane wave expansion method (PWEM) to compute the electronic structure of QD heterostructures, including strain and
piezoelectric effects, and the effect of an arbitrarily oriented electric field. The PWEM allows for the analytical
computation of strain fields, piezoelectric potential, electric field potential, and the Hamiltonian
matrix, for a user-defined calculational supercell. Semiconductors containing
zinc blende crystal structures are treated using an 8-band k·p Hamiltonian.
Once constructed for the supplied input geometry, the supercell Hamiltonian is numerically
diagonalised to obtain the single-particle QD eigenstates. These single-particle eigenstates are
then employed to compute optical (momentum) and Coulomb matrix elements, expressions for which are analytical in the Fourier coefficients of the eigenstate envelope
functions, allowing for efficient numerical evaluation without recourse to numerical differentiation
or integration. By enabling rigorous and efficient calculation of the single-particle QD
eigenstates, and of matrix elements between those eigenstates, QUDOS provides all of the inputs
required for a user to implement semi-classical or many-body calculations of classical and
quantum optical properties.

QUDOS provides a high degree of flexibility, allowing to simulate a variety of heterostructures
including numerous QD shapes both with and without wetting layers, QD molecules, and QD in--a--well structures,
including the ability to treat composition-graded or arbitrarily oriented
QDs. The core physics library is written in C, employs hybrid MPI-OpenMP parallelisation,
and exploits fast-Fourier transform and eigensolver routines from Intel's oneMKL
to accelerate performance.



QUDOS is released under the terms of the GNU General Public License (GPL) v3.0. Alternative commercial licensing terms may be available from the copyright holders upon request.

.. end-docs
