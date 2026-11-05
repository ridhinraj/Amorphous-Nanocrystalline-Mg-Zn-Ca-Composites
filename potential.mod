# ============================================================
# potential.mod -- shared interatomic potential definition.
#
# ORDERING REQUIREMENT: include this AFTER create_box or
# read_data, never before -- pair_coeff needs atom types 1/2/3
# to already exist. (This is why it's a separate file from
# common.mod, which is included at the very top of every script,
# before types exist in the melt-quench script.)
#
# Element ordering: the first list is the order in which library.meam
# indexes the ternary parameters; the second matches the mass
# assignments below (Mg=1, Zn=2, Ca=3). This pairing was verified
# against the HCP Mg lattice parameters, cohesive energy, elastic
# constants, basal GSFE and liquid mixing enthalpy; the results are in
# Table 2 and Section 3.1 of the paper, and the inputs are in
# md/lattice_relax.in, md/elastic/, md/gsfe_basal.in and md/liquid_*.in.
# ============================================================

mass            1 24.305      # Mg
mass            2 65.38       # Zn
mass            3 40.078      # Ca

pair_style      hybrid/overlay meam zbl 1.0 2.0
pair_coeff      * * meam ../potentials/library.meam Ca Zn Mg ../potentials/CaZnMg.meam Mg Zn Ca
pair_coeff      1 1 zbl 12.0 12.0
pair_coeff      1 2 zbl 12.0 30.0
pair_coeff      1 3 zbl 12.0 20.0
pair_coeff      2 2 zbl 30.0 30.0
pair_coeff      2 3 zbl 30.0 20.0
pair_coeff      3 3 zbl 20.0 20.0
