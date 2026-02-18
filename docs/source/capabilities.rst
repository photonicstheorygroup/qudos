
Capabilities
============

QUDOS contains two main executables: **8_band_line.ex** and **8_band_pwem.ex**, which we respectively refer to as the linescan code and electronic structure code. The linescan
code computes linescans and cross-sections of the strain fields, piezoelectric potential, and bulk
band edge energies (i.e. QD band offsets), allowing to visualise these quantities. The purpose
of the linescan code is therefore initial visualisation and sanity checking, and is designed for
use prior to the execution of the electronic structure code. In the electronic structure code,
the strain fields and piezoelectric potential are directly included in the Hamiltonian, which is
exactly diagonalised to obtain the electronic structure. The optical and Coulomb matrix elements are then analytically determined using the calculated single-particle eigenstates.

Theoretical background
~~~~~~~~~~~~~~~~~~~~~~

QUDOS implements a semi-analytical reciprocal space plane wave expansion method (PWEM) :footcite:`Reilly_Springer_2014` to compute the electronic structure of QD heterostructures. According to the PWEM, we define
the envelope functions :math:`F_{nb}(\mathbf{r})` as sums of plane waves, i.e. as a Fourier series:

.. math::

    F_{nb}(\mathbf{r})=\frac{1}{\sqrt{\Omega}}\sum_{\mathbf{G}_{\vec{m}}}\tilde{F}_{nb}(\mathbf{G}_{\vec{m}})e^{i\mathbf{G}_{\vec{m}}\cdot\mathbf{r}}=
    \frac{1}{\sqrt{\Omega}}\sum_{\vec{m}}a_{nb\vec{m}}e^{i\mathbf{G}_{\vec{m}}\cdot\mathbf{r}},

where :math:`\vec{m}=\{m^{(x)}, m^{(y)}, m^{(z)}\}`,  :math:`\mathbf{G}_{\vec{m}}=G_{m^{(x)}}\hat{\mathbf{x}} + G_{m^{(y)}}\hat{\mathbf{y}} + G_{m^{(z)}}\hat{\mathbf{z}}`, 
:math:`G_{m^{(\alpha)}}=\frac{2\pi m^{(\alpha)}}{L^{(\alpha)}}`, :math:`\Omega=L^{(x)}L^{(y)}L^{(z)}` is the supercell volume, and the eigenstates of the 8-band :math:`\mathbf{k}\cdot\mathbf{p}`
Hamiltonian :footcite:`Bahder_PRB_1990` :math:`\psi_{n}(\mathbf{r})` are expressed in terms of the envelope functions :math:`F_{nb}(\mathbf{r})` as

.. math::

    |\psi_{n}(\mathbf{r})\rangle=\sum_{b=1}^{8}F_{nb}(\mathbf{r})|u_{b}\rangle,

where :math:`n` is the subband, and :math:`|u_{b}\rangle` is the zone-centre Bloch state. Substituting the plane wave expression of :math:`F_{nb}(\mathbf{r})` into the Schrödinger equation
then allows the electronic structure calculation to be performed in reciprocal space:

.. math::

    \sum_{b}\sum_{\vec{m}}\left(\hat{H}_{b'b}(\mathbf{r})\right)_{\vec{m}'\vec{m}}a_{nb\vec{m}}=E_{n}a_{nb'\vec{m}'},

where :math:`\left(\hat{H}_{b'b}(\mathbf{r})\right)_{\vec{m}'\vec{m}}` is the matrix representation of the Hamiltonian:

.. math::

    \left(\hat{H}_{b'b}(\mathbf{r})\right)_{\vec{m}'\vec{m}} = \frac{1}{\Omega}\int_{\Omega}e^{-i\mathbf{G}_{\vec{m}'}\cdot\mathbf{r}}\hat{H}_{b'b}(\mathbf{r})e^{i\mathbf{G}_{\vec{m}}\cdot\mathbf{r}}
    \, \text{d}^{3}\mathbf{r}

and is of order :math:`8(2M^{(x)}+1)(2M^{(y)}+1)(2M^{(z)}+1)`, where :math:`2M^{(\alpha)}+1` is the number of plane waves for dimension :math:`\alpha=x,y,z`. The eigenvalues of the matrix produce
the state energies :math:`E_{n}` and the eigenvectors produce the Fourier coefficients :math:`a_{nb\vec{m}}` enabling the calculation of the envelope functions :math:`F_{nb}(\mathbf{r})`.

The operators :math:`\hat{H}_{b'b}(\mathbf{r})` are obtained from the matrix elements :math:`\hat{H}_{b'b}` of the bulk :math:`\mathbf{k}\cdot\mathbf{p}` Hamiltonian by symmetrising and quantising
the components :math:`k_{\alpha}` with respect to the position-dependent material parameters :math:`\mathcal{P}(\mathbf{r})`. For a supercell comprised of a total of :math:`S + 1` distinct shapes, 
:math:`\mathcal{P}(\mathbf{r})` is given by

.. math::
    \mathcal{P}(\mathbf{r}) = \mathcal{P}^{(0)}+\sum_{s=1}^{S}\left(\mathcal{P}^{(s)}-\mathcal{P}_{k}^{(s)}\right)\chi_{s}(\mathbf{r}) \, ,


where :math:`\mathcal{P}^{(0)}` is the value of :math:`\mathcal{P}(\mathbf{r})` in the matrix (i.e. underlying supercell) material, :math:`\mathcal{P}^{(s)}` is the value of :math:`\mathcal{P}(\mathbf{r})` in shape :math:`s`,
:math:`\mathcal{P}_{k}^{(s)}` is the value of :math:`\mathcal{P}(\mathbf{r})` in the shape within which shape :math:`s` is embedded, and :math:`\chi_{s}(\mathbf{r})` is the
characteristic function for shape :math:`s` which is given by

.. math::
    \chi_{s}(\mathbf{r})=
    \begin{cases}
    1, & \text{inside shape } s.\\
    0, & \text{outside shape } s.
    \end{cases}

Calculation of the :math:`\left(\hat{H}_{b'b}(\mathbf{r})\right)_{\vec{m}'\vec{m}}` matrix elements therefore requires calculation of the Fourier transforms of the shape characteristic functions
:math:`\chi_{s}(\mathbf{r})`:

.. math::
    \tilde{\chi}_{s}(\mathbf{G})=\frac{1}{\Omega}\int_{\Omega}\chi_{s}(\mathbf{r})e^{-i\mathbf{G}\cdot\mathbf{r}} \, \text{d}^{3}\mathbf{r} \, ,

expressions of which are analytical in the cases of many shapes. Shapes supported by QUDOS are given in the :ref:`Geometry` page.

As QUDOS utilises a reciprocal space method, the strain, piezoelectric potential, and applied external electric field potential are all calculated using reciprocal space techniques, enabling them
to be directly input when constructing the plane wave Hamiltonian matrix :math:`\left(\hat{H}_{b'b}(\mathbf{r})\right)_{\vec{m}'\vec{m}}`.


Linescan code
~~~~~~~~~~~~~

The linescan code provides the following capabilities:

1.  Calculation of the real-space strain tensor components :math:`\epsilon_{ij}(\mathbf{r})` through the initial use of
    a reciprocal space technique :footcite:`Andreev_JAP_1999` to analytically evaluate :math:`\tilde{\epsilon}_{ij}(\mathbf{G})`. The strain
    calculation can be performed using the common approximation of position-independent
    elastic moduli :math:`C_{ijkl}`, or optionally including the full position-dependent elastic tensor
    :math:`C_{ijkl}(\mathbf{r})` via a perturbative technique :footcite:`Andreev_PRB_2000`. In the case of position-independent
    elastic moduli, the elastic moduli of the matrix material are used. The linescan
    code also employs a larger supercell of volume :math:`\Omega_{S}` for strain calculations, to enable more
    accurate single-QD calculations :footcite:`Vukmirovc_JAP_2008`. In this case, the resulting large-supercell strain tensor component Fourier transforms
    :math:`\tilde{\epsilon}_{ij}^{S}(\mathbf{G})` are then downsampled to obtain :math:`\tilde{\epsilon}_{ij}(\mathbf{G})` on 
    the smaller supercell of volume :math:`\Omega` employed in the electronic structure calculations. The calculated strain
    fields :math:`\epsilon_{ij}(\mathbf{r})` are output as one-dimensional linescans and/or two-dimensional cross-sections.

2.  Calculation of the real space 1\ :sup:`st`\  and 2\ :sup:`nd`\  order piezoelectric potentials :math:`\phi_{\text{pz}}(\mathbf{r})` using
    position-dependent piezoelectric coefficients. The piezoelectric potential is initially calculated in reciprocal space
    :math:`\tilde{\phi}_{\text{pz}}^{S}(\mathbf{G})` using the larger strain supercell (:math:`\Omega_{S} > \Omega`), and then downsampled to 
    obtain :math:`\tilde{\phi}_{\text{pz}}(\mathbf{G})` in the electronic structure supercell. The calculated real-space piezoelectric
    potential :math:`\phi_{\text{pz}}(\mathbf{r})` is output in the electronic structure supercell using one-dimensional linescans and/or
    two-dimensional cross-sections.

3.  Calculation of the 8-band bulk band-edge energies, i.e. the QD band offsets. Using the  computed strain fields :math:`\epsilon_{ij}(\mathbf{r})`
    and piezoelectric potential :math:`\phi_{pz}(\mathbf{r})`, in conjunction with the position-dependent band structure parameters
    and external electrostatic potential :math:`\phi_{\text{ext}}(\mathbf{r})`, the bulk 8-band :math:`\mathbf{k}\cdot\mathbf{p}` Hamiltonian at :math:`\mathbf{k} = 0` is set up
    and diagonalised at real-space grid points :math:`\mathbf{r}` in the electronic structure supercell. The calculated real space
    band-edge energy profiles are output using one-dimensional linescans and/or two-dimensional cross-sections.

4.  Defining a heterostructure produces piecewise continuous material parameters in realspace which, when represented in reciprocal space via Fourier
    transforms of the shape characteristic functions, can produce real-space oscillations (Gibbs phenomenon) in piecewise continuous functions such as
    the strain fields :math:`\epsilon_{ij}(\mathbf{r})` and hence in the band-edge energies. The linescan code supports Gaussian broadening to soften abrupt material parameter profiles,
    allowing to suppress spurious oscillations in strain fields :math:`\epsilon_{ij}(\mathbf{r})` and band edge energies for the purpose of visualisation.


Electronic structure code
~~~~~~~~~~~~~~~~~~~~~~~~~

1.  Inclusion of strain fields :math:`\tilde{\epsilon}_{ij}(\mathbf{G})` and piezoelectric potential :math:`\tilde{\phi}_{\text{pz}}(\mathbf{G})` in the reciprocal space matrix
    representation of the quantised 8-band :math:`\mathbf{k}\cdot\mathbf{p}` Hamiltonian. :math:`\tilde{\epsilon}_{ij}(\mathbf{G})` and :math:`\tilde{\phi}_{\text{pz}}(\mathbf{G})` are
    calculated using the same method as in the linescan code, including the optional use of a larger supercell of volume :math:`\Omega_{S} > \Omega` to better represent isolated single-QD
    structures. The potential of the applied electric field is also included in the supercell Hamiltonian by computing the Fourier transform of the electric field potential
    :math:`\tilde{\phi}_{\text{ext}}(\mathbf{G})`. Furthermore, the analysis of miniband dispersion is also supported for QD superlattice structures via setting of the supercell wave vector :footcite:`Tomic_PRB_2010`.

2.  Calculation of the single-particle eigenstates of the supercell - the energies :math:`E_{n}` and associated envelope functions :math:`\tilde{\psi}_{n}(\mathbf{G})` - via exact
    diagonalisation of the matrix representation of the supercell Hamiltonian.

3.  Calculation of real-space envelope functions :math:`\psi_{n}(\mathbf{r})`, and probability densities :math:`|\psi_{n}(\mathbf{r})|^{2}`,
    as linescans and/or cross-sections along specific directions, plus full three-dimensional outputs on a user-defined real-space grid.

4.  Calculation of the interband momentum matrix elements :math:`\mathbf{p}_{mn}` :footcite:`Szmulowicz_PRB_1995` between states :math:`\psi_{m}(\mathbf{r})` and :math:`\psi_{n}(\mathbf{r})`, computed analytically
    in reciprocal space.

5.  Calculation of the electron-electron :math:`V^{\text{e}}_{ijkl}`, hole-hole :math:`V^{\text{h}}_{ijkl}`, direct electron-hole :math:`V^{\text{dir}}_{ijkl}`, and exchange electron-hole
    :math:`V^{\text{ex}}_{ijkl}` Coulomb matrix elements :footcite:`Saroj_Thesis` between states :math:`\psi_{i}(\mathbf{r})`, :math:`\psi_{j}(\mathbf{r})`, :math:`\psi_{k}(\mathbf{r})`, and :math:`\psi_{l}(\mathbf{r})`,
    computed analytically in reciprocal space. The calculation supports the inclusion of multipole expansion corrections, which mitigate the impact of image charges generated by the periodic
    boundary conditions, and also allows use of a larger supercell than that used to compute the electronic structure, to allow accurate analysis of isolated single-QD structures :footcite:`Vukmirovc_JAP_2008`.


----

Regarding the calculational supercell, both the linescan and electronic structure codes support the following:

1.  A wide variety of QD shapes (see :ref:`Geometry`) are supported by default. Users can construct more complex single-QD morphologies and multi-QD geometries by combining
    multiple shapes in a single supercell. Additionally, other QD shapes can be treated by  determining the Fourier transform :math:`\tilde{\chi}({\mathbf{G}})` of the corresponding
    characteristic function :math:`\chi({\mathbf{r}})`, and adding it as a supported shape to the source code.

2.  Three-dimensional rotation of the supported QD shapes, to facilitate analysis of QDs grown along an arbitrary direction. The code is best suited to analyse heterostructures grown
    along the [001] direction, as it implements the 8-band zinc-blende Hamiltonian in an [001]-oriented total angular momentum basis. To analyse structures grown along other directions,
    the QD shape can be rotated inside the [001]-oriented supercell to align with the intended growth axis. However, it should be noted that the computational efficiency of the calculation
    may not be optimal for structures with different dimensions along the growth and in-plane directions, since the plane wave basis set remains defined with respect to the [100], [010],
    and [001] axes.

----

The electronic structure code supports parallelisation using MPI and/or OpenMP. The optimal hybrid MPI/OpenMP parallelisation scheme for a given calculation will depend on the memory and
processing capabilities of an individual system, and the specifics of a given calculation (e.g. number of shapes included in the calculation supercell, and number of plane waves employed
along each direction). It is generally recommended to employ a maximum of 4 OpenMP threads per MPI process.

The linescan code typically supports larger plane wave basis sets than the electronic structure code, due to its significantly lower memory requirement. The linescan code only supports
OpenMP parallelisation due to the lower computation time requirements.

Typical QUDOS calculations have computational cost exceeding the resources available on standard personal computers (PCs), especially memory usage which can reach 100s of GBs depending
on the number of employed plane waves along each direction, number of shapes, and required speed of calculation. The code is therefore designed to be run on high-performance computing
systems, but it should be possible to run calculations of some structures on high-spec PCs. To assist in determining the amount of random access memory (RAM) required for an electronic
structure calculation, QUDOS also includes a utility executable, **8_band_pwem_memory.ex**, which is intended to be run prior to the electronic structure code using the same input parameters.
It calculates the approximate peak memory that will be used by each MPI process at various stages of the electronic structure code's execution, to provide quantitative guidance on the amount
of RAM required to run a given QUDOS calculation.


**References**

.. footbibliography::

