
Outputs
=======
The following sections describe the output files of QUDOS, for the linescan code and the electronic structure code respectively.



Linescan code
-------------

| **basis_set.dat**
| Outputs information regarding the electronic structure and strain supercell plane wave basis sets, i.e. the number of plane waves used for each dimension.

| **supercell.dat**
| Outputs information regarding the calculation parameters, and the electronic structure and strain supercells, including the shapes within them.

The following calculation parameters are provided: the temperature, the applied electric field, whether the isotropic or anisotropic elastic tensor was used for the strain calculation,
whether position-independent or position-dependent elastic constants were used in the strain calculation, and the order of pieozoelectric potential included.

The following supercell-specific information is provided: material, electronic structure supercell lengths, and strain supercell lengths.

The following information is provided for each shape in the electronic structure supercell: material; shape type; location of the shape centre within the electronic structure supercell;
shape parameters; lattice constant; band gap; spin-orbit splitting energy; modified inverse electron effective mass; kane parameter; first, second, and third modified Luttinger parameters;
unstrained spin-split-off, light hole, heavy hole, and conduction band-edge energies.


| **strain_linescan_###.dat**
| Outputs strain calculation results for the linescan direction ###. The columns are as follows:

| Column 1 - :math:`x`  component of position vector (nm).
| Column 2 - :math:`y`  component of position vector (nm).
| Column 3 - :math:`z`  component of position vector (nm).
| Column 4 - :math:`xx` component of strain tensor   (%).
| Column 5 - :math:`yy` component of strain tensor   (%).
| Column 6 - :math:`zz` component of strain tensor   (%).
| Column 7 - :math:`yz` component of strain tensor   (%).
| Column 8 - :math:`xz` component of strain tensor   (%).
| Column 9 - :math:`xy` component of strain tensor   (%).


| **piezo_linescan_1st_and_2nd_order_###.dat**
| Outputs piezoelectric potential calculation results for the linescan direction ###. Note that "1st_and_2nd_order" assumes that both the 1st and 2nd order were enabled in the piezoelectric potential calculation. The title will include just "1st_order" if only the 1st order is enabled. The columns are as follows:

| Column 1 - :math:`x`    component of position vector         (nm).
| Column 2 - :math:`y`    component of position vector         (nm).
| Column 3 - :math:`z`    component of position vector         (nm).
| Column 4 - Real      part of piezoelectric potential (eV).
| Column 5 - Imaginary part of piezoelectric potential (should be negligible).


| **band_edge_energies_linescan_###.dat**
| Outputs band edge energies calculation results for the linescan direction ###. The columns are as follows:

| Column 1 - :math:`x` component of position vector  (nm).
| Column 2 - :math:`y` component of position vector  (nm).
| Column 3 - :math:`z` component of position vector  (nm).
| Column 4 - spin-split-off band edge energy (eV).
| Column 5 - light-hole     band edge energy (eV).
| Column 6 - heavy-hole     band edge energy (eV).
| Column 7 - conduction     band edge energy (eV).


| **strain_cross_section_###.dat**
| Outputs strain calculation results for the cross-section normal direction ###. The columns are as follows:

| Column 1 - :math:`x`  component of position vector (nm).
| Column 2 - :math:`y`  component of position vector (nm).
| Column 3 - :math:`z`  component of position vector (nm).
| Column 4 - :math:`xx` component of strain tensor   (%).
| Column 5 - :math:`yy` component of strain tensor   (%).
| Column 6 - :math:`zz` component of strain tensor   (%).
| Column 7 - :math:`yz` component of strain tensor   (%).
| Column 8 - :math:`xz` component of strain tensor   (%).
| Column 9 - :math:`xy` component of strain tensor   (%).


| **piezo_cross_section_1st_and_2nd_order_###.dat**
| Outputs piezoelectric potential calculation results for the cross-section normal direction ###. Note that "1st_and_2nd_order" assumes that both the 1st and 2nd order were enabled in the piezoelectric potential calculation. The title will include just "1st_order" if only the 1st order is enabled. The columns are as follows:

| Column 1 - :math:`x`    component of position vector         (nm).
| Column 2 - :math:`y`    component of position vector         (nm).
| Column 3 - :math:`z`    component of position vector         (nm).
| Column 4 - Real      part of piezoelectric potential (eV).
| Column 5 - Imaginary part of piezoelectric potential (should be negligible).


| **band_edge_energies_cross_section_###.dat**
| Outputs band edge energies calculation results for the cross-section normal direction ###. The columns are as follows:

| Column 1 - :math:`x` component of position vector  (nm).
| Column 2 - :math:`y` component of position vector  (nm).
| Column 3 - :math:`z` component of position vector  (nm).
| Column 4 - spin-split-off band edge energy (eV).
| Column 5 - light-hole     band edge energy (eV).
| Column 6 - heavy-hole     band edge energy (eV).
| Column 7 - conduction     band edge energy (eV).



Electronic structure code
-------------------------

| **basis_set.dat**
| Outputs information regarding the electronic structure and strain supercell plane wave basis sets, i.e. the number of plane waves used for each dimension. Also outputs the set supercell miniband dispersion wavevector.

| **supercell.dat**
| Outputs information regarding the calculation parameters, and the electronic structure, strain, and Coulomb supercells, including the shapes within them.

The following calculation parameters are provided: the temperature, the applied electric field, whether the isotropic or anisotropic elastic tensor was used for the strain calculation,
whether position-independent or position-dependent elastic constants were used in the strain calculation, and the order of pieozoelectric potential included.

The following supercell-specific information is provided: material, electronic structure supercell lengths, strain supercell lengths, and Coulomb supercell lengths.

The following information is provided for each shape in the electronic structure supercell: material; shape type; location of the shape centre within the electronic structure supercell;
shape parameters; lattice constant; band gap; spin-orbit splitting energy; modified inverse electron effective mass; kane parameter; first, second, and third modified Luttinger parameters;
unstrained spin-split-off, light hole, heavy hole, and conduction band-edge energies.

| **energies_cb.dat**
| Outputs calculated conduction band state energies, and the corresponding state Bloch character associated with each band. The columns are as follows:

| Column 1  - Energy                                              (eV).
| Column 2  - Bloch character - SO       			     (%).                         
| Column 3  - Bloch character - LH                                (%).
| Column 4  - Bloch character - HH                                (%).
| Column 5  - Bloch character - CB                                (%).

| **energies_vb.dat**
| Outputs calculated valence band state energies, and the corresponding state Bloch character associated with each band. The columns are as follows:

| Column 1  - Energy                                              (eV).
| Column 2  - Bloch character - SO       			     (%).                         
| Column 3  - Bloch character - LH                                (%).
| Column 4  - Bloch character - HH                                (%).
| Column 5  - Bloch character - CB                                (%).

| **envelope_function_cb_x_###.dat**
| Outputs calculated real-space envelope function linescan for conduction band state x along linescan direction ###, where x is indexed according to increasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5  - Real      part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 6  - Imaginary part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 7  - Real      part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 8  - Imaginary part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 9  - Real      part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 10  - Imaginary part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 11  - Real      part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 12  - Imaginary part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 13  - Real      part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 14  - Imaginary part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 15  - Real      part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 16  - Imaginary part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 17  - Real      part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 18 - Imaginary part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 19  - Real      part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 20 - Imaginary part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).


| **envelope_function_vb_x_###.dat**
| Outputs calculated real-space envelope function linescan for valence band state x along linescan direction ###, where x is indexed according to decreasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5  - Real      part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 6  - Imaginary part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 7  - Real      part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 8  - Imaginary part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 9  - Real      part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 10  - Imaginary part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 11  - Real      part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 12  - Imaginary part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 13  - Real      part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 14  - Imaginary part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 15  - Real      part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 16  - Imaginary part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 17  - Real      part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 18 - Imaginary part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 19  - Real      part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 20 - Imaginary part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).


| **probability_density_cb_x_###.dat**
| Outputs calculated real-space probability density linescan for conduction band state x along linescan direction ###, where x is indexed according to increasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5 - Probability density - Total   (:math:`\text{nm}^{-3}`).
| Column 6 - Probability density - SO component only (:math:`\text{nm}^{-3}`).
| Column 7 - Probability density - LH component only (:math:`\text{nm}^{-3}`).
| Column 8 - Probability density - HH component only (:math:`\text{nm}^{-3}`).
| Column 9 - Probability density - CB component only (:math:`\text{nm}^{-3}`).


| **probability_density_vb_x_###.dat**
| Outputs calculated real-space probability density linescan for valence band state x along linescan direction ###, where x is indexed according to decreasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5 - Probability density - Total   (:math:`\text{nm}^{-3}`).
| Column 6 - Probability density - SO component only (:math:`\text{nm}^{-3}`).
| Column 7 - Probability density - LH component only (:math:`\text{nm}^{-3}`).
| Column 8 - Probability density - HH component only (:math:`\text{nm}^{-3}`).
| Column 9 - Probability density - CB component only (:math:`\text{nm}^{-3}`).


| **envelope_function_cross_section_cb_x_###.dat**
| Outputs calculated real-space envelope function cross-section for conduction band state x for cross-section with normal direction ###, where x is indexed according to increasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5  - Real      part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 6  - Imaginary part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 7  - Real      part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 8  - Imaginary part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 9  - Real      part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 10  - Imaginary part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 11  - Real      part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 12  - Imaginary part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 13  - Real      part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 14  - Imaginary part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 15  - Real      part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 16  - Imaginary part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 17  - Real      part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 18 - Imaginary part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 19  - Real      part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 20 - Imaginary part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).


| **envelope_function_cross_section_vb_x_###.dat**
| Outputs calculated real-space envelope function cross-section for valence band state x for cross-section with normal direction ###, where x is indexed according to decreasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5  - Real      part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 6  - Imaginary part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 7  - Real      part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 8  - Imaginary part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 9  - Real      part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 10  - Imaginary part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 11  - Real      part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 12  - Imaginary part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 13  - Real      part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 14  - Imaginary part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 15  - Real      part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 16  - Imaginary part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 17  - Real      part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 18 - Imaginary part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 19  - Real      part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 20 - Imaginary part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).


| **probability_density_cross_section_cb_x_###.dat**
| Outputs calculated real-space probability density cross-section for conduction band state x for cross-section with normal direction ###, where x is indexed according to increasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5 - Probability density - Total   (:math:`\text{nm}^{-3}`).
| Column 6 - Probability density - SO component only (:math:`\text{nm}^{-3}`).
| Column 7 - Probability density - LH component only (:math:`\text{nm}^{-3}`).
| Column 8 - Probability density - HH component only (:math:`\text{nm}^{-3}`).
| Column 9 - Probability density - CB component only (:math:`\text{nm}^{-3}`).


| **probability_density_cross_section_vb_x_###.dat**
| Outputs calculated real-space probability density cross-section for valence band state x for cross-section with normal direction ###, where x is indexed according to decreasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5 - Probability density - Total   (:math:`\text{nm}^{-3}`).
| Column 6 - Probability density - SO component only (:math:`\text{nm}^{-3}`).
| Column 7 - Probability density - LH component only (:math:`\text{nm}^{-3}`).
| Column 8 - Probability density - HH component only (:math:`\text{nm}^{-3}`).
| Column 9 - Probability density - CB component only (:math:`\text{nm}^{-3}`).


| **envelope_function_3D_cb_x.dat**
| Outputs calculated real-space envelope function in three-dimensional space for conduction band state x, where x is indexed according to increasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5  - Real      part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 6  - Imaginary part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 7  - Real      part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 8  - Imaginary part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 9  - Real      part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 10  - Imaginary part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 11  - Real      part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 12  - Imaginary part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 13  - Real      part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 14  - Imaginary part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 15  - Real      part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 16  - Imaginary part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 17  - Real      part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 18 - Imaginary part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 19  - Real      part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 20 - Imaginary part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).


| **envelope_function_3D_vb_x.dat**
| Outputs calculated real-space envelope function in three-dimensional space for valence band state x, where x is indexed according to decreasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5  - Real      part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 6  - Imaginary part of SO S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 7  - Real      part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 8  - Imaginary part of SO S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 9  - Real      part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 10  - Imaginary part of LH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 11  - Real      part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 12  - Imaginary part of LH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 13  - Real      part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 14  - Imaginary part of HH S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 15  - Real      part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 16  - Imaginary part of HH S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 17  - Real      part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 18 - Imaginary part of CB S1 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 19  - Real      part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).
| Column 20 - Imaginary part of CB S2 component of envelope function (:math:`\text{nm}^{-3/2}`).


| **probability_density_3D_cb_x.dat**
| Outputs calculated real-space probability density in three-dimensional space for conduction band state x, where x is indexed according to increasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5 - Probability density - Total   (:math:`\text{nm}^{-3}`).
| Column 6 - Probability density - SO component only (:math:`\text{nm}^{-3}`).
| Column 7 - Probability density - LH component only (:math:`\text{nm}^{-3}`).
| Column 8 - Probability density - HH component only (:math:`\text{nm}^{-3}`).
| Column 9 - Probability density - CB component only (:math:`\text{nm}^{-3}`).


| **probability_density_3D_vb_x.dat**
| Outputs calculated real-space probability density in three-dimensional space for valence band state x, where x is indexed according to decreasing energy. The columns are as follows:

| Column 1  - x position                                          (nm).
| Column 2  - y position                                          (nm).
| Column 3  - z position                                          (nm).
| Column 4  - Energy                                              (eV).
| Column 5 - Probability density - Total   (:math:`\text{nm}^{-3}`).
| Column 6 - Probability density - SO component only (:math:`\text{nm}^{-3}`).
| Column 7 - Probability density - LH component only (:math:`\text{nm}^{-3}`).
| Column 8 - Probability density - HH component only (:math:`\text{nm}^{-3}`).
| Column 9 - Probability density - CB component only (:math:`\text{nm}^{-3}`).


| **overlaps.dat**
| Outputs the calculated spatial overlaps between conduction and valence band states :math:`\left(\langle\psi_{\text{CB}}|\psi_{\text{VB}}\rangle\right)`. The columns are as follows:

| Column 1  - CB state index (in order of increasing energy)
| Column 2  - VB state index (in order of decreasing energy)               
| Column 3  - Overlap real component
| Column 4  - Overlap complex component


| **carrier_localisation_cb.dat**
| Outputs the calculated percentage of the probability density associated with each conduction band eigenstate that resides within a shape in the supercell, i.e. :math:`\int_{\Omega_{s}}\langle\psi_{n}|\psi_{n}\rangle \, \text{d}^{3}\mathbf{r}`, where :math:`\Omega_{s}` is the volume of shape :math:`s`. The rows represent the states in order of increasing energy, and the columns represent the shapes in the order input via supercell.in. The first column gives the percentage inside the electronic structure supercell, but outside any other shape in the supercell. Assuming that shapes are not placed within shapes other than the supercell, the entries in a given row of the output file should sum to 100.0%.


| **carrier_localisation_vb.dat**
| Outputs the calculated percentage of the probability density associated with each valence band eigenstate that resides within a shape in the supercell, i.e. :math:`\int_{\Omega_{s}}\langle\psi_{n}|\psi_{n}\rangle \, \text{d}^{3}\mathbf{r}`, where :math:`\Omega_{s}` is the volume of shape :math:`s`. The rows represent the states in order of decreasing energy, and the columns represent the shapes in the order input via supercell.in. The first column gives the percentage inside the electronic structure supercell, but outside any other shape in the supercell. Assuming that shapes are not placed within shapes other than the supercell, the entries in a given row of the output file should sum to 100.0%.


| **momentum_matrix_elements.dat**
| Outputs the calculated interband momentum matrix elements, i.e. :math:`\widehat{e} \cdot \textbf{p}_{mn} = \frac{ m_{0} }{ \hbar } \langle \psi_{m} ( \textbf{r} ) \vert \widehat{e} \cdot \nabla_{ \textbf{k} } \widehat{H} \vert \psi_{n} ( \textbf{r} ) \rangle`, where :math:`\widehat{e} \cdot \nabla_{ \textbf{k} } \widehat{H}` is the directional derivative of the QD Hamiltonian :math:`\widehat{H}` with respect to wave vector :math:`\textbf{k}`, along the direction specified by the photon polarisation :math:`\widehat{e}`. The momentum matrix elements are given in units of :math:`\text{eV}^{1/2}` so that the squared absolute value of the element is given in terms of :math:`\text{eV}`, akin to the Kane energy. The columns are as follows:

| Column 1  - CB state index (in order of increasing energy)
| Column 2  - VB state index (in order of decreasing energy)			     	                       
| Column 3  - Momentum matrix element for [100] (:math:`\text{eV}^{1/2}`) - Real part
| Column 4  - Momentum matrix element for [100] (:math:`\text{eV}^{1/2}`) - Imag part
| Column 5  - Momentum matrix element for [100] (:math:`\text{eV}^{1/2}`) - Absolute value
| Column 6  - Momentum matrix element for [010] (:math:`\text{eV}^{1/2}`) - Real part
| Column 7  - Momentum matrix element for [010] (:math:`\text{eV}^{1/2}`) - Imag part
| Column 8  - Momentum matrix element for [010] (:math:`\text{eV}^{1/2}`) - Absolute value
| Column 9  - Momentum matrix element for [001] (:math:`\text{eV}^{1/2}`) - Real part
| Column 10  - Momentum matrix element for [001] (:math:`\text{eV}^{1/2}`) - Imag part
| Column 11  - Momentum matrix element for [001] (:math:`\text{eV}^{1/2}`) - Absolute value
| Column 12  - Momentum matrix element for [110] (:math:`\text{eV}^{1/2}`) - Real part
| Column 13  - Momentum matrix element for [110] (:math:`\text{eV}^{1/2}`) - Imag part
| Column 14  - Momentum matrix element for [110] (:math:`\text{eV}^{1/2}`) - Absolute value
| Column 15  - Momentum matrix element for [-110] (:math:`\text{eV}^{1/2}`) - Real part
| Column 16  - Momentum matrix element for [-110] (:math:`\text{eV}^{1/2}`) - Imag part
| Column 17  - Momentum matrix element for [-110] (:math:`\text{eV}^{1/2}`) - Absolute value


| **coulomb_matrix_elements_ee.dat**
| Outputs the calculated Coulomb matrix elements, :math:`V^{e}_{ijkl}`, between electron states, i.e. :math:`V^{e}_{ijkl}=\int_{\Omega}\int_{\Omega} \psi_{i}^{*}(\mathbf{r})\psi_{j}^{*}(\mathbf{r'})V(|\mathbf{r}-\mathbf{r'}|)\psi_{k}(\mathbf{r'})\psi_{l}(\mathbf{r})~\text{d}^{3}\mathbf{r'}~\text{d}^{3}\mathbf{r}`. The columns are as follows:

| Column 1  - Electron state index i (in order of increasing energy)
| Column 2  - Electron state index j (in order of increasing energy)
| Column 3  - Electron state index k (in order of increasing energy)
| Column 4  - Electron state index l (in order of increasing energy)               
| Column 5  - Real part of Coulomb matrix element (eV)
| Column 6  - Imag part of Coulomb matrix element (eV)
| Column 7  - Magnitude of Coulomb matrix element (eV)


| **coulomb_matrix_elements_hh.dat**
| Outputs the calculated Coulomb matrix elements, :math:`V^{h}_{ijkl}`, between hole states, i.e. :math:`V^{h}_{ijkl}=\int_{\Omega}\int_{\Omega} \psi_{i}^{*}(\mathbf{r})\psi_{j}^{*}(\mathbf{r'})V(|\mathbf{r}-\mathbf{r'}|)\psi_{k}(\mathbf{r'})\psi_{l}(\mathbf{r})~\text{d}^{3}\mathbf{r'}~\text{d}^{3}\mathbf{r}`. The columns are as follows:

| Column 1  - Hole state index i (in order of decreasing energy)
| Column 2  - Hole state index j (in order of decreasing energy)
| Column 3  - Hole state index k (in order of decreasing energy)
| Column 4  - Hole state index l (in order of decreasing energy)               
| Column 5  - Real part of Coulomb matrix element (eV)
| Column 6  - Imag part of Coulomb matrix element (eV)
| Column 7  - Magnitude of Coulomb matrix element (eV)


| **coulomb_matrix_elements_eh_direct.dat**
| Outputs the calculated direct Coulomb matrix elements, :math:`V^{\text{dir}}_{ijkl}`, between electron and hole states, i.e. :math:`V^{\text{dir}}_{ijkl}=\int_{\Omega}\int_{\Omega} \psi_{i}^{*}(\mathbf{r})\psi_{j}^{*}(\mathbf{r'})V(|\mathbf{r}-\mathbf{r'}|)\psi_{k}(\mathbf{r'})\psi_{l}(\mathbf{r})~\text{d}^{3}\mathbf{r'}~\text{d}^{3}\mathbf{r}`, where :math:`i, l` are electron states and :math:`j, k` are hole states. The columns are as follows:

| Column 1  - Electron state index i (in order of increasing energy)
| Column 2  - Hole state index j (in order of decreasing energy)
| Column 3  - Hole state index k (in order of decreasing energy)
| Column 4  - Electron state index l (in order of increasing energy)               
| Column 5  - Real part of Coulomb matrix element (eV)
| Column 6  - Imag part of Coulomb matrix element (eV)
| Column 7  - Magnitude of Coulomb matrix element (eV)


| **coulomb_matrix_elements_eh_exchange.dat**
| Outputs the calculated exchange Coulomb matrix elements, :math:`V^{\text{ex}}_{ijkl}`, between electron and hole states, i.e. :math:`V^{\text{ex}}_{ijkl}=\int_{\Omega}\int_{\Omega} \psi_{i}^{*}(\mathbf{r})\psi_{j}^{*}(\mathbf{r'})V(|\mathbf{r}-\mathbf{r'}|)\psi_{k}(\mathbf{r'})\psi_{l}(\mathbf{r})~\text{d}^{3}\mathbf{r'}~\text{d}^{3}\mathbf{r}`, where :math:`i, k` are electron states and :math:`j, l` are hole states. The columns are as follows:

| Column 1  - Electron state index i (in order of increasing energy)
| Column 2  - Hole state index j (in order of decreasing energy)
| Column 3  - Electron state index k (in order of increasing energy)
| Column 4  - Hole state index l (in order of decreasing energy)               
| Column 5  - Real part of Coulomb matrix element (eV)
| Column 6  - Imag part of Coulomb matrix element (eV)
| Column 7  - Magnitude of Coulomb matrix element (eV)

