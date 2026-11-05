#!/usr/bin/env python3
"""
analyze_crystallinity.py

Measure the crystalline (HCP alpha-Mg) fraction of a melt-quench
or devitrification-anneal dump file, using OVITO's Polyhedral
Template Matching modifier. Used to check Week 3 progress toward
target nanocrystal volume fractions (0/10/20/30/40%).

Requires: pip install ovito --break-system-packages

Usage:
    python analyze_crystallinity.py anneal_T450.0_100000.dump
    python analyze_crystallinity.py "anneal_T450.0_*.dump"   # batch, quoted glob
"""
import sys
import glob

from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

# PTM structure type codes (OVITO convention)
STRUCT_NAMES = {0: "Other/amorphous", 1: "FCC", 2: "HCP", 3: "BCC", 4: "ICO"}


def analyze_one(path):
    pipeline = import_file(path)
    pipeline.modifiers.append(PolyhedralTemplateMatchingModifier())
    data = pipeline.compute()

    struct_types = data.particles["Structure Type"]
    n_total = len(struct_types)
    counts = {code: int((struct_types == code).sum()) for code in STRUCT_NAMES}

    print(f"\nFile: {path}   (N = {n_total} atoms)")
    for code, name in STRUCT_NAMES.items():
        frac = counts[code] / n_total
        print(f"  {name:16s}: {frac:7.3%}")

    return counts[2] / n_total  # HCP fraction = alpha-Mg nanocrystal fraction


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    paths = sorted(glob.glob(sys.argv[1])) if "*" in sys.argv[1] else [sys.argv[1]]
    if not paths:
        print(f"No files matched: {sys.argv[1]}")
        sys.exit(1)

    for p in paths:
        analyze_one(p)


if __name__ == "__main__":
    main()
