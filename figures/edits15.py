# -*- coding: utf-8 -*-
"""Style edits applied to the joined paragraph text of the t15_* modules at load time."""
R = [
("Compositions 1 and 9 were resumed from a checkpoint during the first quench, so their effective quench rates were",
 "Compositions 1 and 9 were resumed from a checkpoint during the first quench, which gave effective quench rates of"),
("All packings in a matrix share that matrix’s quench, so the replicates sample seed placement and not independent quenches.",
 "Because all packings in a matrix share that matrix’s quench, the replicates sample seed placement and not independent quenches."),
("The barostat held both lateral stresses at zero, so the stress state approximates free uniaxial tension, and all boundaries remained periodic.",
 "With both lateral stresses held at zero by the barostat, the stress state approximates free uniaxial tension. All boundaries remained periodic."),
("LAMMPS does not divide by any volume, so these raw values have units of stress × volume.",
 "LAMMPS does not divide by any volume, and these raw values therefore have units of stress × volume."),
("Per-atom stress cannot distinguish elastic stress concentration from irreversible rearrangement, so the non-affine squared displacement D²min [@falk] was also computed",
 "Because per-atom stress cannot distinguish elastic stress concentration from irreversible rearrangement, the non-affine squared displacement D²min [@falk] was also computed"),
("The 4.5 Å gap between shell and bulk equals the D²min cutoff, so an atom’s score is not influenced by neighbors in the adjacent region.",
 "The 4.5 Å gap between shell and bulk equals the D²min cutoff, which keeps an atom’s score free of neighbors in the adjacent region."),
("fail the selection criteria, so shells were placed at both radii.", "fail the selection criteria. Shells were therefore placed at both radii."),
("For uncorrelated points ql scales roughly as N_{b}^{−1/2}, so non-averaged values depend on the coordination number.",
 "For uncorrelated points ql scales roughly as N_{b}^{−1/2}, which makes non-averaged values depend on the coordination number."),
("The second-campaign predictions were computed from these frozen models after the runs had finished, so they are retrospective tests.",
 "The second-campaign predictions were computed from these frozen models after the runs had finished and are therefore retrospective tests."),
("The Ca–Ca partial is effectively zero below 3.25 Å, so calcium is dispersed rather than clustered.",
 "The Ca–Ca partial is effectively zero below 3.25 Å, which shows that calcium is dispersed rather than clustered."),
("No steady flow stress is reached, so σ15 is specific to that strain.", "No steady flow stress is reached, and σ15 is specific to that strain."),
("The increase is more than ten times the packing standard deviation of this matrix, so it lies well outside seed-placement scatter.",
 "At more than ten times the packing standard deviation of this matrix, the increase lies well outside seed-placement scatter."),
("The analysis counts clusters interval by interval and does not track individual zones, so zone lifetimes are not known.",
 "The analysis counts clusters interval by interval and does not track individual zones, and zone lifetimes are therefore unknown."),
("(Fig. 12b), so Ca atoms are more strongly separated from one another than Zn atoms are.",
 "(Fig. 12b). Ca atoms are thus more strongly separated from one another than Zn atoms are."),
("shifts the peak strain by up to 0.26%, so the strengths are unaffected, whereas", "shifts the peak strain by up to 0.26%. The strengths are unaffected, whereas"),
("between successive pairs, so the dependence is not log-linear over this range.", "between successive pairs, and the dependence is not log-linear over this range."),
("(errors +0.6% and −1.2%), so on these two tests the surrogate did no better", "(errors +0.6% and −1.2%). On these two tests the surrogate therefore did no better"),
("Simulated rates exceed experimental ones by many orders of magnitude, so only the direction of the rate effect is comparable.",
 "Because simulated rates exceed experimental ones by many orders of magnitude, only the direction of the rate effect is comparable."),
("Each matrix was quenched once, so the ten compositions and the three pure glasses carry no estimate",
 "Each matrix was quenched once, and the ten compositions and the three pure glasses carry no estimate"),
("parallel to the tensile axis, so the composites are effectively textured and", "parallel to the tensile axis. The composites are therefore effectively textured, and"),
("and DBSCAN uses no periodic images, so clusters crossing the cell boundary are split.", "DBSCAN uses no periodic images and splits clusters that cross the cell boundary."),
("The strain rates lie many orders of magnitude above experimental ones, so only the direction of rate and temperature effects is comparable, and the cells are defect-free, which excludes flaw-controlled failure.",
 "The strain rates lie many orders of magnitude above experimental ones, and only the direction of rate and temperature effects is comparable. The cells are defect-free, which excludes flaw-controlled failure."),
]

def apply(mods):
    used = set()
    for m in mods:
        for i, (k, v) in enumerate(m.BLOCKS):
            if k != 'P':
                continue
            for a, b in R:
                if a in v:
                    v = v.replace(a, b); used.add(a)
            m.BLOCKS[i] = (k, v)
    return [a[:60] for a, _ in R if a not in used]
