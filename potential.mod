# ============================================================
# potential.mod -- shared interatomic potential definition.
#
# ORDERING REQUIREMENT: include this AFTER create_box or
# read_data, never before -- pair_coeff needs atom types 1/2/3
# to already exist. (This is why it's a separate file from
# common.mod, which is included at the very top of every script,
# before types exist in the melt-quench script.)
#
# STILL UNVERIFIED: the "Ca Zn Mg" / "Mg Zn Ca" element-order
# pairing below is what composition 1 has been running with
# throughout this project. The second list (Mg Zn Ca) correctly
# matches the mass assignments (Mg=1, Zn=2, Ca=3). The FIRST
# list's order needs to match however CaZnMg.meam internally
# indexes its ternary parameters, and that has never been
# independently confirmed -- run the Phase 1 pure-Mg cohesive
# energy / lattice constant check against this exact pair_style
# line before fully trusting downstream results.
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
