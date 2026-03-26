.. QUDOS documentation master file, created by
   sphinx-quickstart on Mon Feb  9 15:40:01 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Inputs
======

**Linescan code**

The linescan code requires two input files: **supercell.in** and **parameters_linescan.in**, which are described respectively in the :ref:`Geometry` and :ref:`Numerical` sections.

**Electronic structure code**

The electronic structure code requires two input files: **supercell.in** and **parameters.in**, which are described respectively in the :ref:`Geometry` and :ref:`Numerical` sections.
The **supercell.in** file is the same as that required for the linescan code.



.. _Geometry:

Geometry
--------



The geometry, i.e. the supercell and shape parameters, is described via the **supercell.in** file. 

QUDOS uses three differently sized supercells. The first is the electronic structure supercell,
for which the electronic structure is calculated. The second is the larger strain supercell, within which the strain is initially calculated before being downsampled for inclusion in the electronic
structure calculation. This lessons the impact of the inherent periodic boundary conditions for more accurate single QD calculations. The third is the larger Coulomb supercell within which the Coulomb
matrix elements are calculated. This lessens the impact of image charges due to periodic boundary conditions that can be considerable due to the long-range nature of the Coulomb potential. 


The supercell parameters must: (1) be defined prior to the shape parameters, (2) all be included, and (3) be defined in the order given in the following example::

   QD shape: 0 
   Material: 0 
   Compos x: 0.0
   Compos y: 0.0
   Length 1: 60.0
   Length 2: 60.0
   Length 3: 35.0
   Strain length 1: 3
   Strain length 2: 3
   Strain length 3: 4
   Coulomb length 1: 1
   Coulomb length 2: 1
   Coulomb length 3: 2

The descriptions of the supercell parameters are as follows:

**QD shape**: Indicates the shape of the supercell. Only a cuboidal (shape type index 0) supercell shape is currently supported by QUDOS.

**Material**: Indicates the quaternary material of the matrix (or underlying supercell). In the above example, the quaternary material employed is Al\ :sub:`x`\ Ga\ :sub:`1-x`\ P\ :sub:`y`\ As :sub:`1-y`\, which is given by material index 0.
A list of the available quaternary materials, and their identifying indices, is given in the :ref:`Material` section.

**Compos x**: x composition value of the selected quaternary material.

**Compos y**: y composition value of the selected quaternary material.

**Length 1**: Length of the electronic structure supercell along [100] (nm).

**Length 2**: Length of the electronic structure supercell along [010] (nm).

**Length 3**: Length of the electronic structure supercell along [001] (nm).

**Strain length 1**: Length / relative length of the strain supercell along [100]. If using free strain supercell lengths (Strain_supercell_type=0, see :ref:`Numerical` ), this value is the actual length in nm.
If using relative supercell lengths (Strain_supercell_type=1, see :ref:`Numerical` ), this value is an integer that indicates multiples of the electronic structure supercell length along [100] (Length 1).
The above example assumes the use of relative supercell lengths, thus the actual length of the strain supercell along [100] is 3*(Length 1), giving a length of 180.0 nm.

**Strain length 2**: See 'Strain length 1', but for direction [010].

**Strain length 3**: See 'Strain length 1', but for direction [001].

**Coulomb length 1**: Length / relative length of the Coulomb supercell along [100]. If using free Coulomb supercell lengths (Coulomb_supercell_type=0, see :ref:`Numerical` ), this value is the actual length in nm.
If using relative supercell lengths (Coulomb_supercell_type=1, see :ref:`Numerical` ), this value is an integer that indicates multiples of the electronic structure supercell length along [100] (Length 1).
The above example assumes the use of relative supercell lengths, thus the actual length of the Coulomb supercell along [100] is 1*(Length 1), giving a length of 60.0 nm.

**Coulomb length 2**: See 'Coulomb length 1', but for direction [010].

**Coulomb length 3**: See 'Coulomb length 1', but for direction [001].


The shape parameters must: (1) be defined after the supercell parameters, (2) be defined for each shape and therefore the total number must equal to one less than the 'N_shape' numerical input parameter (see :ref:`Numerical`), (3) all be included, and (4) be defined in the order given in the following example::

   Embedded: 0 
   QD shape: 6 
   Material: 3 
   Compos x: 0.0
   Compos y: 0.0
   Centre x: 0.0
   Centre y: 0.0
   Centre z: -3.1
   Length 1: 15.5
   Length 2: 15.5
   Length 3: 7.75
   Length 4: 6.2
   Length 5: 0.0
   Euler alpha angle: 0.0
   Euler beta angle: 0.0
   Euler gamma angle: 0.0


The descriptions of the shape parameters are as follows:

**Embedded**: Index of the shape the current shape is embedded in. The indices are assigned to the shapes in the order in which they are defined in **supercell.in**. Index 0 refers to the supercell.

**QD shape**: Indicates the type of shape. The supported shapes and their identifying indices are given in the :ref:`table-shapes`.

**Material**: Indicates the quaternary material of the shape. In the above example, the quaternary material employed is Ga\ :sub:`x`\ In\ :sub:`1-x`\ P\ :sub:`y`\ As :sub:`1-y`\, which is given by material index 3.
A list of the available quaternary materials, and their identifying indices, is given in the :ref:`Material` section.

**Compos x**: x composition value of the selected quaternary material.

**Compos y**: y composition value of the selected quaternary material.

**Centre x**: Location of the shape centre along [100] (origin at centre of supercell). Each shape type has a defined centre, which are given in the :ref:`table-shapes`.

**Centre y**: Location of the shape centre along [010] (origin at centre of supercell). Each shape type has a defined centre, which are given in the :ref:`table-shapes`.

**Centre z**: Location of the shape centre along [001] (origin at centre of supercell). Each shape type has a defined centre, which are given in the :ref:`table-shapes`.

**Length 1**: Shape parameter 1. See the :ref:`table-shapes`.

**Length 2**: Shape parameter 2. See the :ref:`table-shapes`.

**Length 3**: Shape parameter 3. See the :ref:`table-shapes`.

**Length 4**: Shape parameter 4. See the :ref:`table-shapes`.

**Length 5**: Shape parameter 5. See the :ref:`table-shapes`.

**Euler alpha angle**: Rotation angle about *z* (degrees).*

**Euler beta angle**: Rotation angle about *y'* (degrees).*

**Euler gamma angle**: Rotation angle about *z''* (degrees).*

\*Each shape can be rotated in the supercell to point along an arbitrary direction. The shapes are rotated intrinsically about *z*, *y'*, *z''*, described by proper Euler angles :math:`\alpha, \beta, \gamma`.
For example, to point a shape along [111], :math:`\alpha = 45.0^\circ, \beta = 54.7356^\circ (\text{atan}(\sqrt{2}/1))`. Rotation about [111] is then controlled by :math:`\gamma`.

The following is a table listing the parameters of the supported shapes. "Shape" is the shape type. "Index" is the identifying integer that defines how a shape is selected in the **supercell.in** file via "QD shape".
"L1"..."L5" are the shape parameters, the values of which are set in the **supercell.in** file via "Length 1"..."Length 5". "0x", "0y", "0z" define the location of the shape "centre", where the coordinate system origin is at the centre of the base of the shape.
The shape centre is then positioned in the supercell via the "Centre x", "Centre y", "Centre z" parameters in the **supercell.in file**.


.. _table-shapes:

.. list-table:: Table of shape parameters
   :header-rows: 1
   :widths: 1 1 1 1 1 1 1 1 1 1

   * - Shape
     - Index
     - L1
     - L2
     - L3
     - L4
     - L5
     - 0x
     - 0y
     - 0z
   * - Cuboid
     - 0
     - x length (nm)
     - y length (nm)
     - z length (nm)
     - N/A
     - N/A
     - 0
     - 0
     - L3/2
   * - Sphere
     - 1
     - Radius (nm)
     - N/A
     - N/A
     - N/A
     - N/A
     - 0
     - 0
     - L1
   * - Ellipsoid
     - 2
     - x semi-axis (nm)
     - y semi-axis (nm)
     - z semi-axis (nm)
     - N/A
     - N/A
     - 0
     - 0
     - L3
   * - Hemisphere
     - 3
     - Radius (nm)
     - N/A
     - N/A
     - N/A
     - N/A
     - 0
     - 0
     - 0
   * - Hemiellipsoid (lens)
     - 4
     - x semi-axis
     - y semi-axis
     - z semi-axis
     - N/A
     - N/A
     - 0
     - 0
     - 0
   * - Pyramid (4 sides)
     - 5
     - Half x length (nm)
     - Half y length (nm)
     - Height (nm)
     - N/A
     - N/A
     - 0
     - 0
     - 0
   * - Truncated pyramid (4 sides)
     - 6
     - Base half x length (nm)
     - Base half y length (nm)
     - Top half x length (nm)
     - Height (nm)
     - N/A
     - 0
     - 0
     - 0
   * - N-sided pyramid
     - 7
     - Pyramid base centre to side length (nm)
     - Pyramid height (nm)
     - Number of sides
     - N/A
     - N/A
     - 0
     - 0
     - 0
   * - N-sided truncated pyramid
     - 8
     - Pyramid base centre to side length (nm)
     - Pyramid top centre to side length (nm)
     - Height (nm)
     - Number of sides
     - N/A
     - 0
     - 0
     - 0
   * - Cylinder
     - 9
     - Radius (nm)
     - Height (nm)
     - N/A
     - N/A
     - N/A
     - 0
     - 0
     - L2/2
   * - Cylindrical ring
     - 10
     - Outer radius (nm)
     - Inner radius (nm)
     - Height (nm)
     - N/A
     - N/A
     - 0
     - 0
     - L3/2
   * - Cone
     - 11
     - Radius (nm)
     - Height (nm)
     - N/A
     - N/A
     - N/A
     - 0
     - 0
     - 0
   * - Truncated cone
     - 12
     - Base radius (nm)
     - Top radius (nm)
     - Height (nm)
     - N/A
     - N/A
     - 0
     - 0
     - 0
   * - Truncated conic ring
     - 13
     - Base outer radius (nm)
     - Top outer radius (nm)
     - Inner radius (nm) (:math:`\leq` L2)
     - Height (nm)
     - N/A
     - 0
     - 0
     - 0
   * - Elliptical cone
     - 14
     - x semi-axis (nm)
     - y semi-axis (nm)
     - Height (nm)
     - N/A
     - N/A
     - 0
     - 0
     - 0
   * - Truncated elliptical cone
     - 15
     - Base x semi-axis (nm)
     - Base y semi-axis (nm)
     - Top x semi-axis (nm)
     - Height (nm)
     - N/A
     - 0
     - 0
     - 0
   * - N-sided prism
     - 16
     - Prism centre to side length (nm)
     - Height (nm)
     - Number of sides
     - N/A
     - N/A
     - 0
     - 0
     - L2/2
   * - Quantum well
     - 17
     - Quantum confinement direction (0: [100],  1: [010],  2: [001])
     - Quantum well size (nm)
     - N/A
     - N/A
     - N/A
     - 0
     - 0
     - 0
   * - Quantum wire
     - 18
     - Quantum confinement plane (0: (100),  1: (010),  2: (001) )
     - Length along ([010] if L1 = 0, [100] if L1 = 1, [100] if L1 = 2)
     - Length along ([001] if L1 = 0, [001] if L1 = 1, [010] if L1 = 2)
     - N/A
     - N/A
     - 0
     - 0
     - 0


**Note:**

* The strain and Coulomb supercells must have dimension lengths equal to or larger than those of the electronic structure supercell.

* Shapes must be placed so that they are fully contained within the electronic structure supercell.

* Shapes can be placed fully inside other shapes via the "Embedded" parameter, but not partially. However, note that position-dependent elastic constants are not supported if a shape is placed within another shape.

* A defined shape must not must not touch only one of the supercell boundaries along a Cartesian direction, i.e. if the boundary of a shape is at :math:`-L_{\alpha}`, the other
  boundary of the shape along :math:`\alpha` must also be at :math:`L_{\alpha}`. This condition does not limit calculations in practice due to the periodic boundary conditions - if a shape
  only touches one supercell boundary along a direction, everything in the supercell can be shifted the same distance away from that boundary so that the shape no longer overlaps with either
  of the two boundaries along that direction.

* The first side of the N-sided shapes (pyramid, truncated pyramid, and prism) is always parallel to the [010] direction.

* The quantum well and wire shapes are distinct from a cuboid with relevant side lengths that span the corresponding lengths of the electronic structure supercell, as the quantum well and wire shapes account for the
  use of a larger strain supercell.

* Rotations are not supported for the quantum well and wire shapes.



.. _Numerical:

Numerical parameters
--------------------

The numerical parameters are input to QUDOS via an input file or via command line arguments. If the same parameter is set via the input file and the command line argument,
the command line argument overrides the input file argument. Parameters are set via the command line argument according to "-parameter=value", where "parameter" is the parameter and "value" is the value.
Different input parameters are required for the linescan code and the electronic structure code. In each case, all input parameters must be passed to the code.


----


**Linescan code**

Parameters are read via the input file **parameters_linescan.in**. The input parameters and their descriptions are as follows:

**N_shape**: Number of shapes in the supercell, including the supercell shape.

**M_Gx**: :math:`\text{N}\_\text{Gx}=2\times \text{M}\_\text{Gx}+1=` No. of plane waves used to define the envelope functions in the electronic structure supercell for the [100] dimension. 

**M_Gy**: See **M_Gx**, but for [010]. 

**M_Gz**: See **M_Gx**, but for [001].

**Strain_supercell_type**: Indicates whether to use free strain supercell lengths (0), or to set the lengths to be integer multiples of those of the electronic structure supercell (1) enabling much
faster strain and piezoelectric potential downsampling from the strain supercell to the electronic structure supercell.

**M_S_Gx**: :math:`\text{N}\_\text{S}\_\text{Gx}=2\times \text{M}\_\text{S}\_\text{Gx}+1=` No. of plane waves used to define the strain fields in the strain supercell for the [100] dimension. 

**M_S_Gy**: See **M_S_Gx**, but for [010]. 

**M_S_Gz**: See **M_S_Gx**, but for [001].

**T**: System temperature (K).

**Fx**: Applied electric field [100] component (kV/cm).

**Fy**: Applied electric field [010] component (kV/cm).

**Fz**: Applied electric field [001] component (kV/cm).

**N_Line**: Number of real space grid points for linescans.

**Isotropic**: Indicates whether to use isotropic approximation for elastic tensor. 0 for full anisotropic elastic tensor, 1 for isotropic approximation to elastic tensor.

**Elastic_const**: Indicates whether to use position-dependent elastic constants (approximated via a perturbative technique). 0 for matrix material elastic constants, 1 for position-dependent elastic constants.

**Elastic_const_order**: Correction order in case of position-dependent elastic constants. 

**Piezo**: Order of piezoelectric potential to include in calculation. 0 for no piezoelectricity, 1 for 1st order only, 2 for 1st and 2nd order, 3 for 2nd order only.

**Softening**: Indicates whether to use interface softening via Gaussian convolution. 0 for "hard" interface calculation (prone to Gibb's phenomenon), 1 for interface softening.

**Sigma_softening**: Standard deviation of Gaussian used for interface softening (nm).

**Linescan_100**: Indicates whether to compute linescans along [100]. 0 to disable, 1 to enable.

**Linescan_010**: Indicates whether to compute linescans along [010]. 0 to disable, 1 to enable.

**Linescan_001**: Indicates whether to compute linescans along [001]. 0 to disable, 1 to enable.

**Linescan_110**: Indicates whether to compute linescans along [110]. 0 to disable, 1 to enable.

**Linescan_101**: Indicates whether to compute linescans along [101]. 0 to disable, 1 to enable.

**Linescan_011**: Indicates whether to compute linescans along [011]. 0 to disable, 1 to enable.

**Linescan_111**: Indicates whether to compute linescans along [111]. 0 to disable, 1 to enable.

**Cross_sections**: Indicate whether to output cross-sections of strain, piezoelectric potential, and band-edge energies. 0 to disable, 1 to enable.

**N_cross**: Number of real-space grid points for each of the two cross-section dimensions.

**Cross_section_100**: Indicates whether to compute cross-sections with normal along [100]. 0 to disable, 1 to enable.

**Cross_section_010**: Indicates whether to compute cross-sections with normal along [010]. 0 to disable, 1 to enable.

**Cross_section_001**: Indicates whether to compute cross-sections with normal along [001]. 0 to disable, 1 to enable.

**Cross_section_m110**: Indicates whether to compute cross-sections with normal along [:math:`\bar{1}10`\ ]. 0 to disable, 1 to enable.

**Cross_section_110**: Indicates whether to compute cross-sections with normal along [110]. 0 to disable, 1 to enable.

**Cross_section_m101**: Indicates whether to compute cross-sections with normal along [:math:`\bar{1}01`\ ]. 0 to disable, 1 to enable.

**Cross_section_0m11**: Indicates whether to compute cross-sections with normal along [:math:`0\bar{1}1`\ ]. 0 to disable, 1 to enable.

**Cross_section_111**: Indicates whether to compute cross-sections with normal along [111]. 0 to disable, 1 to enable. Note that for this case, these are expanded cross-sections that include contents of neighbouring supercells, enabling
the visualisation of the effects of the periodic boundary conditions on the environment for this particular direction.

----


**Electronic structure code**

Parameters are read via the input file **parameters.in**. The input parameters and their descriptions are as follows:

**N_shape**: Number of shapes in the supercell, including the supercell shape.

**M_Gx**: :math:`\text{N}\_\text{Gx}=2\times \text{M}\_\text{Gx}+1=` No. of plane waves used to define the envelope functions in the electronic structure supercell for the [100] dimension. 

**M_Gy**: See **M_Gx**, but for [010]. 

**M_Gz**: See **M_Gx**, but for [001].

**Strain_supercell_type**: Indicates whether to use free strain supercell lengths (0), or to set the lengths to be integer multiples of those of the electronic structure supercell (1) enabling much
faster strain and piezoelectric potential downsampling from the strain supercell to the electronic structure supercell.

**M_S_Gx**: :math:`\text{N}\_\text{S}\_\text{Gx}=2\times \text{M}\_\text{S}\_\text{Gx}+1=` No. of plane waves used to define the strain fields in the strain supercell for the [100] dimension. 

**M_S_Gy**: See **M_S_Gx**, but for [010]. 

**M_S_Gz**: See **M_S_Gx**, but for [001].

**N_VB**: Number of valence band states to compute (excluding factor of 2 to account for spin degeneracy).

**N_CB**: Number of conduction band states to compute (excluding factor of 2 to account for spin degeneracy).

**T**: System temperature (K).

**Fx**: Applied electric field [100] component (kV/cm).

**Fy**: Applied electric field [010] component (kV/cm).

**Fz**: Applied electric field [001] component (kV/cm).

**Kx**: Supercell wave vector [100] component (for miniband dispersion) (:math:`\text{nm}^{-1}`). 

**Ky**: Supercell wave vector [010] component (for miniband dispersion) (:math:`\text{nm}^{-1}`). 

**Kz**: Supercell wave vector [001] component (for miniband dispersion) (:math:`\text{nm}^{-1}`). 

**N_Line**: Number of real space grid points for envelope function and probability density linescans.

**Strain**: Indicates whether to include strain in the calculation. 0 to disable, 1 to enable.

**Isotropic**: Indicates whether to use isotropic approximation for elastic tensor. 0 for full anisotropic elastic tensor, 1 for isotropic approximation to elastic tensor.

**Elastic_const**: Indicates whether to use position-dependent elastic constants (approximated via a perturbative technique). 0 for matrix material elastic constants, 1 for position-dependent elastic constants.

**Elastic_const_order**: Correction order in case of position-dependent elastic constants. 

**Piezo**: Order of piezoelectric potential to include in calculation. 0 for no piezoelectricity, 1 for 1st order only, 2 for 1st and 2nd order, 3 for 2nd order only.

**Eigensolver**: Indicates whether to use the minimal workspace version of the eigensolver (0), or the optimal workspace version of the eigensolver (1). The optimal workspace version is recommended unless memory (RAM) is limited.

**Overlaps**: Indicates whether to calculate overlaps between conduction and valence band states. 0 to disable, 1 to enable.

**Localisation**: Indicates whether to calculate carrier localisation of supercell shapes. 0 to disable, 1 to enable.

**MME**: Indicates whether to include calculation of interband momentum matrix elements. 0 to disable, 1 to enable.

**CME**: Indicates whether to include calculation of Coulomb matrix elements with optional multipole expansion corrections.
0 to disable, 1 to enable without corrections, 2 to enable with monopole corrections, 3 to enable with monopole and dipole corrections, 4 to enable with monopole, dipole, and quadrupole corrections.

**Coulomb_supercell_type**: Indicates whether to use free Coulomb supercell lengths (0) or to set them to be integer multiples of those of the electronic structure supercell (1) allowing much faster upsampling.

**M_C_Gx**: :math:`\text{N}\_\text{C}\_\text{Gx}=2\times \text{M}\_\text{C}\_\text{Gx}+1=` No. of plane waves used to define the envelope functions in the Coulomb supercell for the [100] dimension. 

**M_C_Gy**: See **M_C_Gx**, but for [010]. 

**M_C_Gz**: See **M_C_Gx**, but for [001].

**Linescan_100**: Indicates whether to compute linescans along [100]. 0 to disable, 1 to enable.

**Linescan_010**: Indicates whether to compute linescans along [010]. 0 to disable, 1 to enable.

**Linescan_001**: Indicates whether to compute linescans along [001]. 0 to disable, 1 to enable.

**Linescan_110**: Indicates whether to compute linescans along [110]. 0 to disable, 1 to enable.

**Linescan_101**: Indicates whether to compute linescans along [101]. 0 to disable, 1 to enable.

**Linescan_011**: Indicates whether to compute linescans along [011]. 0 to disable, 1 to enable.

**Linescan_111**: Indicates whether to compute linescans along [111]. 0 to disable, 1 to enable.

**Cross_sections**: Indicate whether to output cross-sections of real-space envelope functions and probability densities.

**N_cross**: Number of real-space grid points for each of the two cross-section dimensions.

**Cross_section_100**: Indicates whether to compute cross-sections with normal along [100]. 0 to disable, 1 to enable.

**Cross_section_010**: Indicates whether to compute cross-sections with normal along [010]. 0 to disable, 1 to enable.

**Cross_section_001**: Indicates whether to compute cross-sections with normal along [001]. 0 to disable, 1 to enable.

**Cross_section_m110**: Indicates whether to compute cross-sections with normal along [:math:`\bar{1}10`\ ]. 0 to disable, 1 to enable.

**Cross_section_110**: Indicates whether to compute cross-sections with normal along [110]. 0 to disable, 1 to enable.

**Cross_section_m101**: Indicates whether to compute cross-sections with normal along [:math:`\bar{1}01`\ ]. 0 to disable, 1 to enable.

**Cross_section_0m11**: Indicates whether to compute cross-sections with normal along [:math:`0\bar{1}1`\ ]. 0 to disable, 1 to enable.

**Cross_section_111**: Indicates whether to compute cross-sections with normal along [111]. 0 to disable, 1 to enable. Note that for this case, these are expanded cross-sections that include contents of neighbouring supercells, enabling
the visualisation of the effects of the periodic boundary conditions on the environment for this particular direction.

**3D_output**: Indicates whether to output 3D data of real-space envelope functions and probability densities.

**N_3D**: Number of real-space grid points in each direction for 3D outputs.




.. _Material:

Material parameters
-------------------

The material parameters are set in QUDOS by selecting the quaternary alloys of the matrix (underlying supercell) and supercell shapes from a predefined list by 
setting the shape "Material" parameter in the **supercell.in** file with the index of the material.
The "x" and "y" composition values of the quaternary materials are then set via the "Compos x" and "Compos y" parameters.

QUDOS supports quaternary alloys of the forms :math:`\text{A}_{x}\text{B}_{1-x}\text{C}_{y}\text{D}_{1-y}`,
:math:`\text{A}_{x}\text{B}_{y}\text{C}_{1-x-y}\text{D}`, and :math:`\text{A}\text{B}_{x}\text{C}_{y}\text{D}_{1-x-y}`.
In the case of :math:`\text{A}_{x}\text{B}_{1-x}\text{C}_{y}\text{D}_{1-y}`, the quaternary alloy material parameter :math:`\mathcal{P}(\text{A}_{x}\text{B}_{1-x}\text{C}_{y}\text{D}_{1-y})` is interpolated from the underlying binaries according to:

.. math::
    \mathcal{P}(\text{A}_{x}\text{B}_{1-x}\text{C}_{y}\text{D}_{1-y})&= \frac{x(1-x)\left[y\mathcal{P}(\text{A}_{x}\text{B}_{1-x}\text{C})+(1-y)\mathcal{P}(\text{A}_{x}\text{B}_{1-x}\text{D})\right]}{x(1-x)+y(1-y)} \\
    &+ \frac{y(1-y)\left[x\mathcal{P}(\text{A}_{y}\text{C}_{1-y}\text{D})+(1-x)\mathcal{P}(\text{B}_{y}\text{C}_{1-y}\text{D})\right]}{x(1-x)+y(1-y)} \\
    &-xy(1-x)(1-y)\mathcal{B}(\text{ABCD})\, ,

where

.. math::
  \mathcal{P}(\text{A}_{x}\text{B}_{1-x}\text{C}) = x\mathcal{P}(\text{AB})+(1-x)\mathcal{P}(\text{AC})+x(1-x)\mathcal{B}(\text{ABC}) \, ,

where

:math:`\mathcal{P}(\text{AB})` is the material parameter of the binary :math:`\text{AB}`, :math:`\mathcal{B}(\text{ABC})` is the bowing parameter of the ternary :math:`\text{ABC}` associated with the particular material parameter, 
and :math:`\mathcal{B}(\text{ABCD})` is the bowing parameter of the quaternary :math:`\text{ABCD}` associated with the particular material parameter.
QUDOS assumes :math:`\mathcal{B}(\text{ABC})` and :math:`\mathcal{B}(\text{ABCD})` to be non-zero only in the cases of the average valence band energies, band gap energies, and spin-orbit splitting energies. Similarly, :math:`\mathcal{P}(\text{A}_{x}\text{B}_{y}\text{C}_{1-x-y}\text{D})`
and :math:`\mathcal{P}(\text{A}\text{B}_{x}\text{C}_{y}\text{D}_{1-x-y})` are given by

.. math::
    \mathcal{P}(\text{A}_{x}\text{B}_{y}\text{C}_{1-x-y}\text{D})= x\mathcal{P}(\text{AD}) + y\mathcal{P}(\text{BD}) + (1-x-y)\mathcal{P}(\text{CD}) - xy\mathcal{B}(\text{ABD}) \\
    -x(1-x-y)\mathcal{B}(\text{ACD})-y(1-x-y)\mathcal{B}(\text{BCD})-xy(1-x-y)\mathcal{B}(\text{ABCD}) \, ,

and

.. math::
    \mathcal{P}(\text{A}\text{B}_{x}\text{C}_{y}\text{D}_{1-x-y})= x\mathcal{P}(\text{AB}) + y\mathcal{P}(\text{AC}) + (1-x-y)\mathcal{P}(\text{AD}) - xy\mathcal{B}(\text{ABC}) \\
    -x(1-x-y)\mathcal{B}(\text{ABD})-y(1-x-y)\mathcal{B}(\text{ACD})-xy(1-x-y)\mathcal{B}(\text{ABCD}) \, .

The following is a table of the quaternary materials supported by QUDOS, and their associated indices:

.. _table-quaternaries:

.. list-table:: Table of quaternary materials
   :widths: 50 50
   :header-rows: 1

   * - Material
     - Index
   * - :math:`\text{Al}_{x}\text{Ga}_{1-x}\text{P}_{y}\text{As}_{1-y}`
     - 0
   * - :math:`\text{Al}_{x}\text{Ga}_{1-x}\text{P}_{y}\text{Sb}_{1-y}`
     - 1
   * - :math:`\text{Al}_{x}\text{Ga}_{1-x}\text{As}_{y}\text{Sb}_{1-y}`
     - 2
   * - :math:`\text{Ga}_{x}\text{In}_{1-x}\text{P}_{y}\text{As}_{1-y}`
     - 3
   * - :math:`\text{Ga}_{x}\text{In}_{1-x}\text{P}_{y}\text{Sb}_{1-y}`
     - 4
   * - :math:`\text{Ga}_{x}\text{In}_{1-x}\text{As}_{y}\text{Sb}_{1-y}`
     - 5
   * - :math:`\text{Al}_{x}\text{In}_{1-x}\text{P}_{y}\text{As}_{1-y}`
     - 6
   * - :math:`\text{Al}_{x}\text{In}_{1-x}\text{P}_{y}\text{Sb}_{1-y}`
     - 7
   * - :math:`\text{Al}_{x}\text{In}_{1-x}\text{As}_{y}\text{Sb}_{1-y}`
     - 8
   * - :math:`\text{Al}_{x}\text{Ga}_{y}\text{In}_{1-x-y}\text{P}`
     - 9
   * - :math:`\text{Al}_{x}\text{Ga}_{y}\text{In}_{1-x-y}\text{As}`
     - 10
   * - :math:`\text{Al}_{x}\text{Ga}_{y}\text{In}_{1-x-y}\text{Sb}`
     - 11
   * - :math:`\text{Al}\text{P}_{x}\text{As}_{y}\text{Sb}_{1-x-y}`
     - 12
   * - :math:`\text{Ga}\text{P}_{x}\text{As}_{y}\text{Sb}_{1-x-y}`
     - 13
   * - :math:`\text{In}\text{P}_{x}\text{As}_{y}\text{Sb}_{1-x-y}`
     - 14


The material parameters of the binary materials :math:`\mathcal{P}(\text{AB})`, the ternary bowing parameters :math:`\mathcal{B}(\text{ABC})`, and the quaternary bowing parameters :math:`\mathcal{B}(\text{ABCD})`
are set via the **material_parameters.h** file located in the code/source/parameters directory. For :math:`\mathcal{B}(\text{ABCD})`, the arrays are indexed according to the :ref:`table-quaternaries`.
For :math:`\mathcal{B}(\text{ABC})`, the arrays are indexed according to:

.. _table-ternaries:

.. list-table:: Table of ternary materials
   :widths: 50 50
   :header-rows: 1

   * - Material
     - Index
   * - :math:`\text{Al}_{x}\text{Ga}_{1-x}\text{P}`
     - 0
   * - :math:`\text{Al}_{x}\text{Ga}_{1-x}\text{As}`
     - 1
   * - :math:`\text{Al}_{x}\text{Ga}_{1-x}\text{Sb}`
     - 2
   * - :math:`\text{Ga}_{x}\text{In}_{1-x}\text{P}`
     - 3
   * - :math:`\text{Ga}_{x}\text{In}_{1-x}\text{As}`
     - 4
   * - :math:`\text{Ga}_{x}\text{In}_{1-x}\text{Sb}`
     - 5
   * - :math:`\text{Al}_{x}\text{In}_{1-x}\text{P}`
     - 6
   * - :math:`\text{Al}_{x}\text{In}_{1-x}\text{As}`
     - 7
   * - :math:`\text{Al}_{x}\text{In}_{1-x}\text{Sb}`
     - 8
   * - :math:`\text{Al}\text{P}_{x}\text{As}_{1-x}`
     - 9
   * - :math:`\text{Al}\text{P}_{x}\text{Sb}_{1-x}`
     - 10
   * - :math:`\text{Al}\text{As}_{x}\text{Sb}_{1-x}`
     - 11
   * - :math:`\text{Ga}\text{P}_{x}\text{As}_{1-x}`
     - 12
   * - :math:`\text{Ga}\text{P}_{x}\text{Sb}_{1-x}`
     - 13
   * - :math:`\text{Ga}\text{As}_{x}\text{Sb}_{1-x}`
     - 14
   * - :math:`\text{In}\text{P}_{x}\text{As}_{1-x}`
     - 15
   * - :math:`\text{In}\text{P}_{x}\text{Sb}_{1-x}`
     - 16
   * - :math:`\text{In}\text{As}_{x}\text{Sb}_{1-x}`
     - 17

For :math:`\mathcal{P}(\text{AB})`, the arrays are indexed according to:

.. _table-binaries:

.. list-table:: Table of binary materials
   :widths: 50 50
   :header-rows: 1

   * - Material
     - Index
   * - :math:`\text{AlN}`
     - 0
   * - :math:`\text{AlP}`
     - 1
   * - :math:`\text{AlAs}`
     - 2
   * - :math:`\text{AlSb}`
     - 3
   * - :math:`\text{AlBi}`
     - 4
   * - :math:`\text{GaN}`
     - 5
   * - :math:`\text{GaP}`
     - 6
   * - :math:`\text{GaAs}`
     - 7
   * - :math:`\text{GaSb}`
     - 8
   * - :math:`\text{GaBi}`
     - 9
   * - :math:`\text{InN}`
     - 10
   * - :math:`\text{InP}`
     - 11
   * - :math:`\text{InAs}`
     - 12
   * - :math:`\text{InSb}`
     - 13
   * - :math:`\text{InBi}`
     - 14


It is recommended to modify the material parameters to values that are most applicable to your application of interest. Note that when the **material_parameters.h** file is modified, QUDOS
must be recompiled (see :ref:`Compilation`) for the changes to take effect. Make sure that **8_band_line.ex** and **8_band_pwem.ex** have been deleted prior to recompiling to ensure that
the compilation proceeds correctly.