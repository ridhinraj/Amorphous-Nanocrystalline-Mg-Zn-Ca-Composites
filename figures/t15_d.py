# -*- coding: utf-8 -*-
# Draft 15, part D: discussion, conclusions, declarations

BLOCKS = [

("H1", "4. Discussion"),

("H2", "4.1. The calcium effect"),

("P", "Peak strength rose with Ca content at each of three Zn levels, and the pure-glass controls show the same +17.9% "
      "rise between xCa = 0.03 and 0.07 without any seed. The calcium effect is therefore a property of the glass "
      "matrix and does not come from seed geometry. At the atomic scale, Ca and Zn atoms sit in measurably different "
      "environments. Ca atoms have larger Voronoi volumes, longer bonds, lower non-averaged q6 and a stronger "
      "same-type deficit. After normalization by their own volume, they carry 1.5–1.6 times the von Mises stress of "
      "Zn atoms from the unloaded glass onward, and they take 1.7–1.9 times their share of the load increment."),

("P", "The regression narrows down which features matter. A structure-only model reproduces 81–88% of the "
      "volume-normalized Ca–Zn difference in mean log stress while explaining only 15–18% of the per-atom variance. "
      "The descriptors that carry it are the Voronoi volume, bond length and coordination, all measures of local "
      "spacing, whereas bond-orientational and chemical order add almost nothing. The excess is thus more consistent "
      "with a local-size effect than with an effect of packing order or chemical order. This remains a correlation. "
      "Structure and identity are strongly collinear, since a larger atom necessarily has a larger cell and longer "
      "bonds, and the data do not show that the higher stress on Ca atoms causes the higher strength of Ca-rich "
      "glasses. That simple static descriptors predict per-atom stress weakly also says nothing about more expressive "
      "structural fields, such as softness or graph-network predictors, which anticipate where rearrangements occur in "
      "glasses [@cubuk; @bapst] and could be applied to this system."),

("P", "The effect found here is distributed through the amorphous matrix. It differs from the role of calcium in "
      "crystalline Mg-Zn-Ca alloys, where strengthening has been attributed to solute segregation at twin and grain "
      "boundaries [@basu; @jang22]. The simulations do not exclude further contributions from crystalline boundaries "
      "or from the interphase, and both kinds of effect may coexist in real microstructures."),

("H2", "4.2. Zinc content, volume fraction and the interphase"),

("P", "Zinc raised the peak strength much less than calcium per unit composition (+6.2% from xZn = 0.25 to 0.38) and "
      "also raised the density. Two of the four compositions in the Zn series (1 and 9) were quenched slightly "
      "faster. Their peak strengths lie 0.6–0.7% above linear interpolation between their nominal-rate neighbors, "
      "comparable to the packing scatter, so the Zn trend is not attributed to quench rate, although a contribution "
      "below 1% cannot be excluded. No per-atom stress or structural signature comparable to the calcium effect was "
      "found for zinc."),

("P", "Peak strength and modulus both fell as the volume fraction rose, although the defect-free crystal is far "
      "stiffer and stronger than the glass. As discussed in Section 3.4, this points to a compliant interphase produced "
      "by the construction protocol rather than to a weak crystal, and a connected crystalline path is excluded by "
      "construction. Peak strength also increased with seed radius at fixed volume fraction and follows VF/R, which "
      "scales with interphase area per unit volume. This again points to the interphase rather than the crystal "
      "volume as the quantity that lowers the strength, although it rests on three radii in one matrix. Within this "
      "protocol, nanocrystals lowered σmax while reducing post-peak softening and the number of simultaneously active "
      "rearrangement clusters, and the Ca advantage persisted at every volume fraction."),

("P", "Hypothesis (ii) is only partly supported. Atoms in the 3 Å shell showed a modest but statistically clear "
      "excess of non-affine displacement over the bulk at every strain tested. The excess was largest before the peak "
      "and declined as the material yielded, which is the opposite of what nucleation at yielding would predict, and "
      "its size of 20–40% is far below what a dominant nucleation site would produce. A larger interphase area per "
      "unit volume lowers the strength while the non-affine excess stays modest, which suggests a distributed rather "
      "than a concentrated interphase contribution. Whether the volume-fraction trends survive a better-relaxed "
      "interphase is the main open question, and a below-Tg anneal of the composites would address it."),

("H2", "4.3. Character of plastic deformation"),

("P", "Hypothesis (i) is not supported in these cells. No discrete shear band was detected in any of the 46 glass and "
      "composite simulations, and three features of the simulations favor homogeneous flow. The cell edge of about "
      "20 nm limits the space available to a band and its process zone. The glasses were quenched at about 10¹² K/s, "
      "and rapidly quenched model glasses deform more homogeneously than well-relaxed ones [@shi05]. The cells are "
      "fully periodic, with no free surface or notch to nucleate a band. In addition, the stress-based detector was "
      "not checked against a banding configuration (Section 2.7). The absence of a band is therefore a statement "
      "about these simulation conditions and does not show that the composites resist shear banding."),

("P", "Within these conditions, plasticity proceeded through a bounded population of small, spatially separate "
      "rearrangements, consistent with shear-transformation-zone descriptions [@falk; @argon]. This population was "
      "insensitive to Ca content, despite a 16% change in peak strength, but decreased with volume fraction in the "
      "same direction as strength and softening. A direct test of hypothesis (i) requires larger cells, better-relaxed "
      "glasses made by slower quenching or longer annealing, and a strain-based localization measure [@shimizu]."),

("H2", "4.4. Strain rate, temperature and the surrogate model"),

("P", "Peak strength rose with strain rate and fell with temperature, consistent with thermally activated "
      "rearrangement, since both change the thermal energy available within the loading time. The rate effect was the "
      "same at four compositions that differ by 13% in strength (+10.3% per doubling, standard deviation 0.4 "
      "percentage points). MD of Cu₆₄Zr₃₆ glass between 10⁷ and 10⁹ s⁻¹ likewise gives higher strength at higher "
      "rates, whereas experiments on metallic glasses report increasing, decreasing or unchanged strength depending on "
      "composition and loading mode [@fan24]. Simulated rates exceed experimental ones by many orders of magnitude, so "
      "only the direction of the rate effect is comparable."),

("P", "The surrogate identified temperature and strain rate, both fixed across the first thirteen composition and "
      "volume-fraction points, as the largest gaps in the design space, and the recommended follow-up runs were "
      "carried out. Where the inputs lay within the training range, 16 of 17 simulated values fell inside the 95% "
      "intervals, with a mean absolute error of 0.99%. On the two out-of-range volume fractions, however, it did no "
      "better than linear extrapolation, and it could not anticipate the temperature dependence, which was never varied "
      "in its training data. Its proposal to vary temperature and rate was also trivial, because an input that has not "
      "been varied always has the largest uncertainty. With a design this small, and nearly one-dimensional within "
      "each series, simple trend models are adequate. The surrogate would become useful only for denser, "
      "multi-variable designs with an anisotropic kernel. No computational-savings claim is made, and we regard the "
      "surrogate as a demonstration of uncertainty-ranked follow-up selection and held-out checking, not as a "
      "validated design tool."),

("H2", "4.5. Implications for alloy design"),

("P", "Within the limits listed in Section 4.6, modest increases in Ca content are a candidate route to higher peak "
      "strength in Mg-Zn-Ca glasses and AN composites, and Ca is already a common biocompatible addition in this "
      "system [@vinogradov]. The effect belongs to the glass matrix, persists at every volume fraction tested and is "
      "not an artifact of seed geometry, although it rests on three quenches and the atomic-scale analysis is "
      "correlational. Nanocrystals made by the present construction lowered the peak strength while reducing post-peak "
      "softening and the number of active rearrangement clusters. Any benefit of the AN microstructure is therefore "
      "more likely to lie in deformation stability and damage tolerance than in peak strength [@li26]. At fixed volume "
      "fraction, fewer and larger seeds were stronger, which suggests keeping the interphase area small. These "
      "conclusions depend on the interphase structure and need confirmation with relaxed interphases. The study "
      "examined mechanical response only, in a defect-free 20 nm cell, and does not demonstrate a gain in damage "
      "tolerance. Corrosion, degradation kinetics and biocompatibility were outside its scope."),

("P", "The simulated peak strengths (0.86–1.38 GPa across the glass and composite simulations) and moduli (26–33 GPa) "
      "cannot be compared directly with measurements on cast Mg-Zn-Ca glasses. The strain rate is many orders of "
      "magnitude higher than in laboratory tests. The 20 nm defect-free cell cannot contain the casting flaws that "
      "control the measured strength of bulk Mg-Zn-Ca glasses [@zhao08]. The experimental materials are also "
      "monolithic glasses, not AN composites. An order-of-magnitude comparison is still useful. Glassy Mg₆₇Zn₂₈Ca₅ "
      "wires tested in tension fracture at 675–894 MPa, with a Weibull characteristic strength of 817 MPa and a "
      "Weibull modulus of 20.6 [@zberg09w]. Cast bulk glasses fail in compression at 716–854 MPa (Mg₆₆Zn₃₀Ca₄) and "
      "672–752 MPa (Mg₇₁Zn₂₅Ca₄), with Weibull moduli of 26 and 41 [@zhao08], and a monolithic glass of our baseline "
      "composition, Mg₆₀Zn₃₅Ca₅, fails in compression at 571 MPa [@wong17]. Our pure-glass peak strengths "
      "(1.07–1.26 GPa) are therefore 1.3–1.5 times the characteristic tensile strength of the wires, 1.3–1.9 times "
      "the flaw-limited compressive strengths and 2.1 times the strength of a cast glass of our baseline composition. "
      "This is the margin expected of a defect-free cell loaded at 10⁹ s⁻¹."),

("P", "The moduli differ more. Bulk Mg₆₅₋₇₂Zn₂₃₋₃₀Ca₄₋₅ glasses give 48.8–50.4 GPa [@jin22], whereas the pure glasses "
      "here give 29.2–32.5 GPa, only 0.58–0.67 times the measured stiffness. The crystalline elastic benchmarks of "
      "Section 3.1 cannot detect this shortfall, since they compare the potential only with tabulated values for HCP "
      "Mg."),

("H2", "4.6. Limitations"),

("P", "Cell size and glass state. The cell edge is about 20 nm, the glasses were quenched at about 10¹² K/s and the "
      "cells have no free surface. All three favor homogeneous flow, so the absence of a band in all 46 glass and "
      "composite cells does not show that bands cannot form in larger or better-relaxed cells. Whether the bounded "
      "population of rearrangement clusters would coalesce into a band remains untested."),

("P", "Replication. Each matrix was quenched once, so the ten compositions and the three pure glasses carry no "
      "estimate of quench-to-quench variability. Seed-placement scatter was measured with four packings at 10% volume "
      "fraction in each of three matrices. These packings share their matrix’s quench and loading velocity seed, and "
      "scatter at saturated volume fractions, at other radii and for the rate and temperature series was not "
      "measured. The rate and temperature series reuse one relaxed structure per composition, and all p-values are "
      "conditional on one realization."),

("P", "Construction. In the ten-composition sweep the seed radius rises with calcium content, a confound that the "
      "pure-glass controls resolve. The seed-radius dependence rests on three radii in one matrix, with one run each. "
      "The placement algorithm saturates near 18–19% volume fraction. The seeds are perfect HCP spheres inserted into "
      "the glass, so the interphase is a product of the construction and not of crystallization; it was created by "
      "deleting matrix atoms within 2.0 Å of the seed and relaxed for 50 ps, and the falling modulus suggests that it "
      "is compliant. Interphase energy, free volume and seed crystallinity during loading were not quantified. All "
      "seeds share one orientation, with the basal plane parallel to the tensile axis, so the composites are "
      "effectively textured and the effect of seed orientation was not tested. Compositions 1 and 9 were quenched "
      "faster than nominal, and an effect below about 1% on their peak strengths cannot be excluded."),

("P", "Interatomic potential. The potential was checked for lattice, cohesive and elastic properties, for the basal "
      "stacking-fault energies along both paths and for the sign of the liquid mixing enthalpy, but not against "
      "experiment for plasticity or for glass mechanics. Against first-principles values for the I₂ fault it places "
      "the intrinsic stacking-fault energy about 20% low and the unstable stacking-fault energy about 40% high. With "
      "the present seed orientation basal slip is not driven, so this error matters less here than it would for other "
      "orientations, but prismatic and pyramidal slip were not benchmarked. Pure Zn is weakly bound in this potential, "
      "which is the main uncertainty in the mixing enthalpy, and the simulated glass moduli lie well below measured "
      "values (Section 4.5)."),

("P", "Analysis. Shear-band detection used stress rather than atomic strain and had no positive control. The "
      "absolute D²min cluster count depends several-fold on the analysis settings and is treated as an index, and "
      "DBSCAN uses no periodic images, so clusters crossing the cell boundary are split. Interphase ratios depend on "
      "where the shell is placed and were evaluated with one shell width on two composites, and before the peak D²min "
      "includes elastic and thermal non-affine motion. The structural analyses used three of the ten compositions. "
      "Non-averaged q6 depends on coordination, and one neighbor cutoff was used for all species pairs. Structure and "
      "identity are collinear, per-atom R² did not exceed 0.275, and the random atom-wise split is optimistic, so these "
      "findings are correlational. The surrogate follow-up was a single round."),

("P", "Scope. Only mechanical response was studied. The strain rates lie many orders of magnitude above experimental "
      "ones, so only the direction of rate and temperature effects is comparable, and the cells are defect-free, which "
      "excludes flaw-controlled failure."),

("H1", "5. Conclusions"),

("P", "Forty-eight MD tensile simulations of Mg-Zn-Ca glasses and AN composites lead to the following conclusions, "
      "which hold for 20 nm periodic cells of glasses quenched at about 10¹² K/s and loaded at 10⁸–2×10⁹ s⁻¹."),

("P", "(1) Peak strength increased with Ca content, by 17.9% from xCa = 0.03 to 0.07 in pure glasses at xZn = 0.35, "
      "so the calcium effect is a property of the glass matrix. Zinc had a weaker effect (+6.2% from xZn = 0.25 to "
      "0.38)."),

("P", "(2) Nanocrystals lowered the peak strength by 1.16% per percentage point of volume fraction in three matrices, "
      "against a packing scatter of 0.3–0.5%, and reduced post-peak softening from 23–24% to 13%. Young’s modulus "
      "also fell although the crystal is stiffer than the glass, which points to compliant interphases left by the "
      "construction. At fixed volume fraction a larger seed gave a stronger composite (+6.0% from 12 to 20 Å)."),

("P", "(3) Peak strength fell by 0.245% per K between 250 and 400 K and rose by 10.3% per doubling of strain rate at "
      "four compositions."),

("P", "(4) No discrete shear band was detected, and hypothesis (i) is not supported in these cells. Incremental "
      "non-affine displacement showed a steady population of small rearrangement clusters, similar across Ca content "
      "and 26–56% smaller in the composites than in the glass. Cumulative D²min integrates thermal and anelastic "
      "motion and would wrongly suggest a single system-spanning cluster."),

("P", "(5) The interphase showed a 20–40% excess of non-affine motion over the bulk, largest before the peak. It is a "
      "site of modestly enhanced response rather than the dominant site of rearrangement, so hypothesis (ii) is only "
      "partly supported."),

("P", "(6) Ca atoms carried about twice the raw per-atom von Mises stress of Zn atoms from the unloaded glass onward. "
      "One third of this is the larger atomic volume of Ca. The remaining 1.5–1.6-fold excess in volume-normalized "
      "stress is predicted to 81–88% by local structure, mainly volume, bond length and coordination, and the result "
      "is correlational."),

("P", "(7) A Gaussian process surrogate placed 16 of 17 informative held-out simulations inside its 95% intervals, "
      "but it did no better than linear extrapolation on two out-of-range volume fractions and could not anticipate "
      "the temperature dependence."),

("P", "Future work should use larger cells, slower quenches and independent quenches of each matrix, relaxed or "
      "crystallization-derived interphases, several seed orientations and a strain-based localization measure."),

]

DECL = [
("H1", "CRediT authorship contribution statement"),
("P", "Ridhin V Raj: Conceptualization, Methodology, Software, Validation, Writing – original draft. "
      "Aparna Zagabathuni: Conceptualization, Resources, Supervision, Writing – review & editing."),
("H1", "Declaration of competing interest"),
("P", "The authors declare that they have no known competing financial interests or personal relationships that "
      "could have appeared to influence the work reported in this paper."),
("H1", "Acknowledgements"),
("P", "The authors acknowledge the High Performance Computing Cluster of the Central Computer Center, National "
      "Institute of Technology Calicut, for the computational resources used in this work."),
("H1", "Funding"),
("P", "This research did not receive any specific grant from funding agencies in the public, commercial, or "
      "not-for-profit sectors."),
("H1", "Data availability"),
("P", "The LAMMPS input files, the analysis and figure scripts (Python), the stress–strain curves and per-case "
      "property table of all 48 simulations, the per-atom Voronoi volumes used in Section 3.6 and the held-out "
      "prediction tables are openly available in a public repository [@repo]. The raw trajectories (dump files) are "
      "large and are available from the corresponding author on reasonable request. The interatomic potential is "
      "available from the NIST Interatomic Potentials Repository (record 2019--Jang-H-S-Seol-D-Lee-B-J--Mg-Zn-Ca) "
      "[@becker; @hale]."),
("H1", "Declaration of generative AI and AI-assisted technologies in the manuscript preparation process"),
("P", "During the preparation of this work the authors used Claude (Anthropic) in order to help write and check the "
      "analysis and figure scripts, check arithmetic and cross-references, and draft and edit the text. After using "
      "this tool, the authors reviewed and edited the content as needed and take full responsibility for the content "
      "of the published article."),
]
