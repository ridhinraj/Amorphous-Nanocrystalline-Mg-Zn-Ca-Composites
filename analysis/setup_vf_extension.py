"""
setup_vf_extension.py

Generalizes the volume-fraction sweep beyond run_1: any of the 10
compositions' already-quenched glass can be reused as the source for
additional VF cases, and the seed radius is no longer fixed at 12 Ang.

This directly tests two open limitations at once:
  1. Whether the VF trend (yield stress decreasing 5% -> 15% VF)
     generalizes beyond the single composition it's currently
     confirmed on.
  2. Whether it generalizes beyond the single (12 Ang) seed geometry
     it's currently confirmed on.

Cost structure is identical to the original vf_1-3: since every
composition already completed Stage 1 (melt-quench) months ago, each
new case here only needs Stage 2-4 (build, relax, tension) --
~1 day each, not the ~2 weeks a new composition would cost.

Reuses BUILDER_TEMPLATE, RELAX_TEMPLATE, and TENSION_TEMPLATE directly
from the already-tested setup_job_array.py, rather than retyping them,
so every accumulated fix (pd.RangeIndex, cleared velocities, frozen
L0, y/z NPT coupling through the pulling run) carries over automatically.

Run from the project root:
    module load apps/anaconda3/2024.06   # or shap_env, either has what this needs
    python3 setup_vf_extension.py
"""

import os
import sys

sys.path.insert(0, os.getcwd())
from setup_job_array import BUILDER_TEMPLATE, RELAX_TEMPLATE, TENSION_TEMPLATE

# ============================================================
# Extension cases -- each is (case_name, source_run, seed_radius, target_vf).
# source_run is whichever composition's already-quenched glass to
# reuse; target_vf is the fraction (0.05 = 5%).
#
# Prioritized default set, given the cost of each (~1 day):
#   - run_6 at the SAME 12 Ang radius as the original vf_1-3, at all
#     three original VF targets -- tests composition-generality
#     directly comparable to the existing baseline-composition result.
#   - Two new seed-radius points (15, 20 Ang) at the baseline
#     composition (run_1) and a single, shared 10% VF target --
#     tests seed-geometry-generality without a full second 3-point
#     sweep at each radius.
#
# This is a deliberately small, prioritized extension (5 new cases,
# ~5 days total), not a full 9-composition x N-radius cross product
# (which would be 25+ cases and likely doesn't fit before submission).
# Add more tuples below if there's schedule room for more.
# ============================================================
EXTENSION_CASES = [
    ("vf_r6_1", "run_6", 12.0, 0.05),
    ("vf_r6_2", "run_6", 12.0, 0.10),
    ("vf_r6_3", "run_6", 12.0, 0.15),
    ("vf_rad15_1", "run_1", 15.0, 0.10),
    ("vf_rad20_1", "run_1", 20.0, 0.10),
]


def main():
    for case_name, source_run, seed_radius, target_vf in EXTENSION_CASES:
        os.makedirs(case_name, exist_ok=True)

        amorphous_file = f"../{source_run}/quenched_matrix_300K.data"
        rng_seed = 3000 + abs(hash((source_run, seed_radius, target_vf))) % 1000

        builder_script = (
            BUILDER_TEMPLATE
            .replace("__AMORPHOUS_FILE__", amorphous_file)
            .replace("__SEED_RADIUS__", str(seed_radius))
            .replace("__TARGET_VF__", str(target_vf))
            .replace("__RNG_SEED__", str(rng_seed))
        )
        with open(f"{case_name}/build_an_composite.py", "w") as f:
            f.write(builder_script)

        with open(f"{case_name}/in.relax_composite", "w") as f:
            f.write(RELAX_TEMPLATE)

        with open(f"{case_name}/week4_tension.in", "w") as f:
            f.write(TENSION_TEMPLATE)

        # Written so submit_vf_extension_array.sh can look up the
        # correct dependency file per case without hardcoding a
        # single source composition (unlike the original
        # submit_vf_array.sh, which only ever checked run_1).
        with open(f"{case_name}/source_composition.txt", "w") as f:
            f.write(source_run)

        print(f"{case_name}: source={source_run}, seed_radius={seed_radius} Ang, "
              f"target_vf={target_vf:.0%}")

    print(f"\nGenerated {len(EXTENSION_CASES)} extension cases: "
          f"{', '.join(c[0] for c in EXTENSION_CASES)}")
    print("Each depends on its source composition's quenched_matrix_300K.data "
          "already existing -- all 10 compositions confirmed complete, so this "
          "should be immediately runnable.")


if __name__ == "__main__":
    main()
