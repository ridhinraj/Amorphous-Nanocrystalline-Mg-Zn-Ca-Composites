# ============================================================
# common.mod -- core settings shared by every LAMMPS script in
# this project. Include this FIRST, before create_box/read_data.
# ============================================================

units           metal
dimension       3
boundary        p p p
atom_style      atomic
atom_modify     map array

# Centralized here so every stage (quench, relax, tension) uses
# the same neighbor settings
neighbor        4.0 bin
neigh_modify    delay 0 every 1 check yes