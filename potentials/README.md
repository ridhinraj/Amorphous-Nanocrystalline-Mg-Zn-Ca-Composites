# Interatomic potential

The simulations use the 2NN MEAM Mg-Zn-Ca potential of Jang, Seol and Lee (2019), overlaid with ZBL
for short-range repulsion. The potential files are distributed by NIST and are not redistributed
here.

Download record **2019--Jang-H-S-Seol-D-Lee-B-J--Mg-Zn-Ca** from the NIST Interatomic Potentials
Repository (https://www.ctcms.nist.gov/potentials/) and place the two files in this folder:

```
potentials/library.meam
potentials/CaZnMg.meam
```

`potential.mod` in the project root refers to them as `../potentials/...`, because every case runs in
its own subdirectory.

The potential as set up here was checked against the HCP Mg lattice parameters, the cohesive energy,
the six independent elastic constants, the basal generalized stacking-fault energy and the liquid
mixing enthalpy. The results are in Table 2 and Section 3.1 of the paper; the inputs are in
`md/lattice_relax.in`, `md/elastic/`, `md/gsfe_basal.in` and `md/liquid_*.in`.
