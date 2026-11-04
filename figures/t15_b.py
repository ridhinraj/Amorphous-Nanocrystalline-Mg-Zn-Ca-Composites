# -*- coding: utf-8 -*-
# Draft 15, part B: results
F = "/tmp/w"

BLOCKS = [

("H1", "3. Results"),

("H2", "3.1. Validation of the interatomic potential"),

("P", "Energy minimization of HCP Mg gave a = 3.209 Å and c/a = 1.619, within 0.03% and 0.25% of experiment, and a "
      "cohesive energy of −1.549 eV/atom, 2.6% larger in magnitude than the experimental −1.51 eV/atom [@kittel] "
      "(Table 2). The relaxed-ion elastic constants satisfy hexagonal symmetry without it being imposed "
      "(C66 = (C11 − C12)/2 within 0.005 GPa) and agree within a few percent with the values tabulated by NIST for "
      "this potential [@becker; @hale]. This checks the workflow rather than the accuracy of the potential against "
      "experiment. The Hill average [@hill] of these constants gives a polycrystalline Young’s modulus of 48.3 GPa "
      "for α-Mg."),

("P", "On the basal plane the GSFE maximum along the straight ⟨112̅0⟩ path is 320.7 mJ/m², and the curve returns to "
      "the undisplaced energy at a full Burgers vector to within 10⁻⁹ eV. This maximum cannot be compared with "
      "first-principles values, which refer to the path along which the basal dislocation splits into two partials. "
      "Along that partial path the same procedure gives an intrinsic stacking-fault energy of 30.0 mJ/m², at one third "
      "of the lattice period, and an unstable stacking-fault energy of 125.6 mJ/m². First-principles calculations for "
      "the I₂ fault in pure Mg give 37.6 and 89.6 mJ/m² [@dong]. The potential therefore places the stable fault about 20% "
      "low and the barrier about 40% high. Basal slip is not driven in the present seeds, however. The tensile axis "
      "lies in the basal plane of every seed (Section 2.3), so the resolved shear stress on the basal plane is zero, "
      "and crystal plasticity in the seeds would require prismatic or pyramidal slip."),

("P", "The liquid mixing enthalpy of Mg₆₀Zn₃₅Ca₅ at 1000 K, referenced to liquid Mg, Zn and Ca prepared with the "
      "same protocol, is −7.5 kJ/mol (three seeds; standard deviation 0.05 kJ/mol). It is negative, as expected for a "
      "glass-forming alloy, and the ternary liquid is 5.7% denser than an ideal mixture of the three. This is the only "
      "benchmark in which Ca interacts with both Zn and Mg, and its sign is consistent with the ternary parameters "
      "being read as intended. Pure Zn is weakly bound in this potential and evaporates at 2000 K and zero pressure. "
      "The Zn reference, which carries 35% of the weight, is therefore the main uncertainty, and only the sign and "
      "order of magnitude of this value are relied on. The partial g(r) distances (Section 3.2) increase in the order "
      "Zn–Zn < Mg–Mg < Ca–Ca, consistent with the metallic radii, which is a further qualitative check of the ternary "
      "parameters."),

("P", "Two defect-free HCP Mg cells loaded in tension give an upper reference for the crystal. Their peak stress, an "
      "ideal strength, was 2.51 GPa along ⟨2̅110⟩ and 4.71 GPa along [0001], with moduli of 43.8 and 54.2 GPa "
      "(Supplementary Table S1). The pure glasses reach 43–50% of the first and 23–27% of the second. Together, the "
      "benchmarks cover lattice, cohesive and elastic properties, both basal stacking faults and the sign of the "
      "mixing enthalpy. They do not establish the accuracy of the potential for plastic deformation of crystalline "
      "Mg."),

("TBL", ("new:T2", None,
         "Table 2.|Potential benchmarks. Reference values: experiment [@kittel]; first-principles (DFT) basal "
         "stacking-fault energies on the minimum-energy path [@dong]; elastic constants tabulated by NIST for this "
         "potential [@becker; @hale].")),

("H2", "3.2. Structure and thermal history of the amorphous matrices"),

("P", "The total g(r) of every quenched matrix shows a single broad first peak at 3.075–3.125 Å, a split second peak "
      "near 5.3 and 6.1 Å, which is typical of metallic glasses [@cheng], and damped oscillations that converge to "
      "g(r) ≈ 1 (Fig. 3a). Across the ten compositions g(10 Å) was 0.956–1.034, with no sharper, non-decaying "
      "correlations of the kind residual crystallinity would give (Supplementary Table S3). The 0.05 Å spread of "
      "first-peak positions equals one histogram bin and is not interpreted. The same-species partials peak at "
      "2.725 Å (Zn–Zn), 3.175 Å (Mg–Mg) and 3.825 Å (Ca–Ca) (Fig. 3b), an order that follows the metallic radii. The "
      "Ca–Ca partial is effectively zero below 3.25 Å, so calcium is dispersed rather than clustered. No "
      "structure-identification analysis was applied to the matrices themselves, so a very small crystalline fraction "
      "cannot be excluded."),

("FIG", (F + "/recol14/fig02.png", 6.2,
         "Fig. 3.|(a) Total radial distribution function of the quenched matrix of Composition 4, with a single broad "
         "first peak, a split second peak and decay toward g(r) = 1. (b) Same-species partial radial distribution "
         "functions (Mg–Mg, Zn–Zn, Ca–Ca).")),

("P", "Density increased linearly with zinc content, from 2.579 g/cm³ at xZn = 0.25 to 3.203–3.209 g/cm³ at "
      "xZn = 0.40 (R² = 0.9995; Fig. 4), and varied by at most 0.3% with Ca content at fixed xZn (Supplementary "
      "Table S3). This follows from the atomic masses, since zinc (65.38 g/mol) is 2.7 times heavier than "
      "magnesium."),

("FIG", (F + "/recol14/fig03.png", 4.2,
         "Fig. 4.|Density of the quenched matrices versus zinc content. Marker size encodes calcium content, which has "
         "no resolved effect on density (at most 0.3%). Dashed line: linear fit, R² = 0.9995.")),

("P", "On cooling at about 10¹² K/s, Tg from the potential energy was 515–535 K across the ten compositions (mean "
      "526 K), and from the specific volume 546–563 K (mean 555 K) (Supplementary Table S4). The two differ because "
      "the liquid branch is curved. For Composition 3, widening the fitting window from 600 to 1200 K raised the "
      "volume-based value from 489 to 603 K and moved the energy-based value non-monotonically between 488 and 530 K "
      "(Supplementary Table S5). Over the three windows used, the energy-based Tg of a composition varied by at most "
      "14 K and the volume-based value by up to 57 K. The energy-based values are therefore used, and differences of "
      "about 20 K between compositions are not interpreted."),

("FIG", (F + "/recol14/fig04.png", 6.2,
         "Fig. 5.|Thermal history of Composition 3. (a, b) Potential energy and specific volume against temperature for "
         "the first quench (gray), the reheating to 600 K (orange) and the second cooling (dark brown). The shaded band "
         "is the range of Tg from the first quench over fitting windows, and the dashed line marks the 600 K stage. "
         "(c) Potential energy during the 500 ps hold at 600 K; black bars are 100 ps means.")),

("P", "The 600 K stage therefore lies 65–85 K above the energy-based Tg. It re-equilibrates the glass as a "
      "supercooled liquid and quenches it a second time at the same rate (0.999–1.001×10¹² K/s), rather than "
      "annealing it below Tg (Fig. 5a, b). The potential energy reached a plateau during the hold (drift of −2.2 to "
      "+1.0 meV/atom per ns; Fig. 5c), with no step that would indicate crystallization. The cycle lowered the 300 K "
      "potential energy by 0.70–1.81 meV/atom, raised the density by 0.03–0.09% and lowered the energy-based Tf by "
      "6–17 K, to a final 506–519 K with no resolvable composition dependence. Compositions 1 and 9, quenched slightly "
      "faster, have Tg (526 and 529 K) and first-quench Tf (524 and 526 K) within the ranges of the other eight. No "
      "effect of the faster quench is resolved, although with one composition per rate a small effect cannot be "
      "excluded. The effect of Tf on strength was not tested."),

("H2", "3.3. Peak strength and composition dependence"),

("P", "All 46 glass and composite curves (Fig. 6) show a near-linear rise, a broad stress maximum at 5.8–8.0% strain "
      "(6.5–7.1% for the ten compositions) and strain softening that continues to 15% strain. No steady flow stress "
      "is reached, so σ15 is specific to that strain. Because the maximum is broad, εmax is less precisely defined "
      "than σmax. Each matrix was quenched once, and the compositions are compared at one quench each."),

("FIG", (F + "/fig15/fig05_stress_strain.png", 5.5,
         "Fig. 6.|Tensile stress–strain curves grouped by the variable changed: (a) pure glasses at three calcium "
         "contents; (b–d) nanocrystal volume-fraction series in the matrices with xCa = 0.03, 0.05 and 0.07, with the "
         "pure glass dashed; (e) seed radius at 10% volume fraction; (f) strain rate; (g) temperature; (h) defect-free "
         "HCP Mg crystals along [0001] and ⟨2̅110⟩. Smoothed curves are solid and raw data faint; dots mark the peak "
         "strength. The ten-composition sweep is shown in Fig. 7. Each series uses one color family in all figures.")),

("P", "The three pure-glass controls separate the composition effect from any influence of the seeds. Their peak "
      "strengths are 1.069, 1.177 and 1.260 GPa at xCa = 0.03, 0.05 and 0.07 (xZn = 0.35), steps of +10.2% and +7.0% "
      "and +17.9% overall (Fig. 6a, Fig. 7b; Table 3). Young’s modulus rises by 11.4% over the same range, and the "
      "softening ratio stays between 22.6% and 24.2%. The calcium effect is therefore a property of the glass matrix. "
      "Each glass is a single quench, and quench-to-quench variability is not sampled (Section 4.6)."),

("P", "The ten-composition sweep agrees with this (Table 3). In that sweep the seed radius rises with Ca content (15 "
      "to 35 Å at xZn = 0.35; Table 1), so the sweep alone cannot separate composition from seed geometry. The "
      "controls show that this does not matter here. The single-seed composites lie within 1.0% of their own glasses "
      "at seed volume fractions of 0.16–2.02% and reproduce the glass steps to within two percentage points. The step "
      "from xCa = 0.05 to 0.07 is also recovered at matched geometry, with 12 Å seeds at equal volume fraction "
      "(+7.8%, +6.5% and +6.7% at 5%, 10% and 15%; Section 3.4)."),

("P", "Within the sweep, calcium raised the peak strength at each of three zinc levels (Fig. 7b), by 8.7% from "
      "xCa = 0.04 to 0.06 at xZn = 0.30, by 8.2% at xZn = 0.40, and by 16.3% across xCa = 0.03–0.07 at xZn = 0.35. "
      "Every slice increases monotonically, and σ15 follows the same order. Young’s modulus also increased with Ca, "
      "from 29.5 to 32.4 GPa at xZn = 0.35, whereas the softening ratio showed no clear trend. Raising Zn at "
      "xCa = 0.05 increased σmax monotonically but by much less, from 1.105 to 1.174 GPa (+6.2%) between xZn = 0.25 "
      "and 0.38 (Fig. 7a)."),

("FIG", (F + "/fig15/fig06_composition.png", 6.0,
         "Fig. 7.|Peak strength (filled markers) and flow stress at 15% strain (hollow markers) versus (a) Zn content "
         "at xCa ≈ 0.05 and (b) Ca content at three Zn levels. Diamonds in (b) are the pure-glass controls at "
         "xZn = 0.35.")),

("TBL", ("t14", 2,
         "Table 3.|Mechanical properties of the ten compositions and of the three pure-glass controls (xZn = 0.35). "
         "Composition numbers follow Table 1. Rows are grouped as the Zn series, the Ca pairs and the pure glasses "
         "(VF = 0).")),

("H2", "3.4. Nanocrystal volume fraction and seed radius"),

("P", "Seeds of 12 Å radius were added at 5%, 10% and 15% volume fraction to each of the three matrices, and the "
      "xCa = 0.05 matrix was also taken to its packing limit (Fig. 6b–d; Table 4). Each composite is compared with the "
      "pure glass of its own matrix. Peak strength fell almost linearly with volume fraction in every matrix, by "
      "5.4–7.1% at 5%, 10.9–12.5% at 10% and 15.6–17.9% at 15% (Fig. 8a). Fitted through the origin, the pooled slope "
      "is −1.16% per percentage point (−1.06% to −1.21% in the three matrices), and the three matrices collapse onto "
      "one line in relative terms (Fig. 8b). The two saturated cases at 18.5% and 18.9% lie close to this line "
      "extended beyond the fitted range."),

("P", "Post-peak softening decreased with volume fraction (Fig. 8c), from 22.6–24.2% in the glasses to 12.8% at "
      "18.5–18.9%. The flow stress at 15% strain therefore changed less than the peak. In other words, the "
      "nanocrystals lower the peak strength while reducing the stress drop that follows it. Seed number, and with it "
      "interphase area, rises with volume fraction at fixed radius, so these trends do not separate crystal volume "
      "from interphase area. The calcium advantage persisted at every volume fraction, at +16.6%, +15.7% and +14.7% "
      "between xCa = 0.03 and 0.07, against +17.9% in the glass, and the two effects act almost independently."),

("P", "Young’s modulus also decreased with volume fraction, from 29.2–32.5 GPa in the glasses to 26.9–28.5 GPa at "
      "15% (a drop of 8–12%; Table 4). A stiffer second phase should do the opposite. Along the loading direction the "
      "defect-free crystal has a modulus of 43.8 GPa (Section 3.1). Even the lower (Reuss) mixture bound for 15% "
      "crystal in the xCa = 0.05 glass is 32.2 GPa, whereas the composite gives 27.8 GPa. The falling modulus "
      "therefore points to a compliant region associated with the seeds. The most likely candidate is the interphase "
      "created by construction, where matrix atoms within 2.0 Å of the seed were deleted and the structure was relaxed "
      "for only 50 ps. The loss of peak strength with volume fraction may share this origin, and the interphase excess "
      "of non-affine motion (Section 3.5) is consistent with it. The interphase energy and free volume, and the "
      "crystallinity of the seeds during loading, were not quantified."),

("P", "To measure how large a difference seed placement alone can produce, four packings per matrix were loaded "
      "identically at 10% volume fraction. Their peak strengths have coefficients of variation of 0.32–0.51% (pooled "
      "0.41%), and the highest and lowest of the four differ by 0.76–1.17%. The four packings share their matrix’s "
      "quench and loading velocity seed, so this is the scatter due to seed placement and not the full run-to-run "
      "variability. It is an order of magnitude smaller than the change over the volume-fraction range and than the "
      "calcium steps."),

("P", "At 10% volume fraction in the xCa = 0.05 matrix, peak strength increased with seed radius: 1.035, 1.056 and "
      "1.097 GPa at 12, 15 and 20 Å, or +6.0% overall (Fig. 8d). The small differences in achieved volume fraction "
      "(10.01%, 9.94% and 9.81%) account for less than 0.003 GPa of this. The increase is more than ten times the "
      "packing standard deviation of this matrix, so it lies well outside seed-placement scatter. The flow stress is "
      "not monotonic in radius (0.854, 0.843 and 0.871 GPa), however, and each radius was run once. At fixed volume "
      "fraction the interphase area per unit volume scales as VF/R, and a larger seed carries less interphase for the "
      "same crystal volume. The six points of the 12 Å series, including the pure glass, follow "
      "σmax = 1.170 − 0.160 (VF/R) GPa, with VF in percent and R in Å. Fitted without the 15 and 20 Å cases, this "
      "relation predicts them to within 0.008 GPa. This is a consistency check rather than a demonstration, but it "
      "points to interphase area, not crystal volume, as the quantity that lowers the strength."),

("P", "Targets of 20% and 25% reached the limit of the placement algorithm, and both saturated, at 18.91% and "
      "18.49%. Their peak strengths differ by 1.2% in the direction opposite to the volume-fraction trend, which is "
      "more than the packing scatter measured at 10%. Scatter at these saturated fractions was not measured, and the "
      "pair is not interpreted."),

("FIG", (F + "/fig15/fig07_volume_fraction.png", 6.3,
         "Fig. 8.|(a) Peak strength versus achieved volume fraction in the three matrices (12 Å seeds). Large open "
         "symbols at VF = 0 are the pure glasses, small open symbols the single-seed composites of Section 3.3, and bars "
         "the standard deviation of four packings at 10%. (b) Change relative to the matched pure glass; the dotted "
         "line is the pooled slope of −1.16% per percentage point. (c) Post-peak softening from the peak to 15% "
         "strain. (d) Peak strength versus seed radius at 10% volume fraction in the xCa = 0.05 matrix; the gray band "
         "is ±1 standard deviation of the 12 Å packings.")),

("TBL", ("t14", 3,
         "Table 4.|Pure-glass, volume-fraction, packing-replicate and seed-radius cases. VF is the achieved nanocrystal "
         "volume fraction from each build log, and Δσmax is the change relative to the pure glass of the same matrix. "
         "Rows marked n = 4 give the mean and sample standard deviation of four packings at 10% volume fraction.")),

("H2", "3.5. Character of plastic deformation"),

("P", "Under hypothesis (i), strength would be set by the nucleation and arrest of a discrete, localized shear band. "
      "The von Mises stress clustering of Section 2.7 found no such band in any recorded frame of any of the 46 glass "
      "and composite simulations (Supplementary Table S6). This covers a 20-fold range of strain rate, temperatures of "
      "250–400 K, volume fractions up to 19% and the pure glasses. The test uses stress rather than strain and was "
      "not checked on a cell known to band, so it shows that no strongly localized high-stress region forms, not that "
      "flow is homogeneous. Non-affine displacement was therefore computed for the three compositions of the "
      "xZn = 0.35 series, for the 5%, 10% and 15% composites of the xCa = 0.05 matrix and for the three pure "
      "glasses."),

("P", "Cumulative D²min integrates all non-affine displacement since the undeformed frame. With the threshold "
      "calibrated on the first 0.5% strain, the number of distinct above-threshold clusters rose to a peak of 273–458 "
      "at 4.5–6% strain and then fell as the clusters merged (Fig. 9a). By 15% strain, 64–70% of all atoms were above "
      "the threshold, and the largest cluster held at least 99.9% of them. This reflects thermal, anelastic and "
      "plastic non-affine motion accumulated over the whole loading history. It shows that most of the material has "
      "moved non-affinely by 15% strain, but it cannot locate the activity at a given strain. In a synthetic control in "
      "which one group of atoms was displaced once and then moved only thermally, the cumulative D²min of that group "
      "kept rising at every later frame. A zero-strain trajectory of the same length, which would isolate the thermal "
      "part in the real cells, was not run. Incremental D²min measures the activity within each 0.5% interval and is "
      "used for all cluster results below."),

("FIG", (F + "/b15/fig09_d2min.png", 6.2,
         "Fig. 9.|(a) Number of above-threshold D²min clusters in Composition 1. Cumulative D²min (each frame compared "
         "with the undeformed frame) integrates all non-affine displacement and merges into one cluster; incremental "
         "D²min (each frame compared with the preceding one) measures current activity. (b) Peak number of "
         "simultaneous incremental-D²min clusters for the three compositions (xZn = 0.35) and the three volume "
         "fractions of the xCa = 0.05 matrix.")),

("P", "In every configuration the number of simultaneously active incremental clusters rose from near zero to a level "
      "that it then held for the rest of the loading (Fig. 9b; Table 5). It first reached 80% of its maximum 2–4 "
      "percentage points after the peak strain, and no single cluster ever dominated; the largest anywhere held 182 "
      "atoms. Averaged over strains of 10% and above, the count was similar for the three compositions (134, 118 and "
      "133 for xCa = 0.03, 0.05 and 0.07) and for the three pure glasses (124, 126 and 129), although σmax changes by "
      "16% over this range. It decreased with volume fraction (93, 78 and 56 at 5%, 10% and 15%), to 26%, 38% and 56% "
      "below the matched glass. These counts, σmax and softening come from the same configurations and are not "
      "independent. At the peak of Composition 1 the rearranging atoms form small, scattered groups (Fig. 2b)."),

("TBL", ("t14", 4,
         "Table 5.|Incremental D²min cluster populations at the default settings (cutoff 4.5 Å, threshold mean + 4 SD, "
         "DBSCAN eps 6 Å, minimum cluster size 15). “Mean” is over all frames at strain ≥ 10%, and “strain at 80% of "
         "peak” is the first frame at which the count reaches 80% of its maximum. Peak strains are 6.6–7.1% (Table 3). "
         "Each row is one trajectory.")),

("P", "The absolute count depends strongly on the analysis settings and is an index, not a number of physical zones. "
      "For Composition 1 it ranged from 47 to 241 over the threshold multipliers tested and from 7 to 289 over the "
      "minimum cluster sizes (Supplementary Table S7). Three results held in all ten non-degenerate settings, since a "
      "4 Å clustering radius finds almost no clusters. The counts fell with volume fraction, the number of flagged "
      "atoms, which does not depend on the clustering, fell in the same order (in all 11 settings), and no cluster "
      "dominated. The counts of the three compositions differed by 11–23% in every moderate setting. The analysis "
      "counts clusters interval by interval and does not track individual zones, so zone lifetimes are not known. A "
      "bounded population of many small rearrangements is consistent with shear-transformation-zone descriptions of "
      "amorphous plasticity [@falk; @argon]."),

("P", "Hypothesis (ii) was tested directly on every seed of B5 (59 seeds) and B15 (177 seeds) at 2%, 4%, 7% and 10% "
      "strain; their peak strains are 6.9% and 6.7%. With the shell at the 12 Å build radius, the mean "
      "interphase-to-bulk ratio of incremental D²min was 1.36, 1.40, 1.25 and 1.22 at the four strains in B5, and "
      "1.32, 1.33, 1.27 and 1.20 in B15 (Fig. 10a; Supplementary Table S8). The 95% interval lay above 1 at every "
      "strain (p ≤ 6×10⁻⁹), and most seeds were individually above 1 (51–59 of 59 and 144–173 of 177). Between 4% "
      "and 7% strain the ratio fell by 0.15 in B5 (paired t-test p = 0.002; 41 of 59 seeds decreased; Fig. 10b) and "
      "by 0.06 in B15 (p = 0.007), with no change between 2% and 4%. Wilcoxon tests gave the same conclusions. With "
      "the shell at the recovered radius (9.7 Å) the ratios were higher (1.54, 1.58, 1.37 and 1.30 in B5). With a "
      "bulk reference that includes the interphases of neighboring seeds they were lower (1.20 and 1.14 at 7% strain). "
      "The magnitude therefore depends on these choices, but the sign and the decline with strain do not."),

("P", "Before the peak, D²min includes elastic and thermal non-affine motion. The pre-peak excess therefore shows that "
      "the interphase responds more non-affinely than the bulk, not that it rearranges irreversibly earlier. A 20–40% "
      "excess is far below the 25-fold ratio recovered for an engineered interphase-only rearrangement (Section 2.9). "
      "The interphase is thus a site of modestly enhanced non-affine response, not the dominant site of plastic "
      "activity. Whether the excess is intrinsic or comes from the way the interphase was built (Section 3.4) cannot "
      "be decided from these data."),

("FIG", (F + "/fig15/fig10_interphase.png", 6.2,
         "Fig. 10.|(a) Mean interphase-to-bulk incremental D²min ratio (95% confidence interval) against strain over all "
         "seeds of B5 (59 seeds) and B15 (177 seeds), with the shell at the 12 Å build radius; the dashed line repeats "
         "B5 with each seed’s recovered radius. The gray band marks the peak strain. (b) The 59 seeds of B5 at 4% and "
         "7% strain; gray lines are individual seeds and the green line is the mean.")),

("H2", "3.6. Local environment and stress of Ca and Zn atoms"),

("P", "Model A, built on identity and 3.25 Å neighbor counts, predicted each Zn and Ca atom’s raw von Mises stress "
      "with a held-out R² of 0.18–0.23, and the atom’s own type dominated the SHAP ranking (Supplementary Table S9). "
      "In all thirteen configurations analyzed at 15% strain, the ten single-seed compositions and the three "
      "volume-fraction composites of the xCa = 0.05 matrix, Ca atoms had a higher mean log stress than Zn atoms. The "
      "difference was 0.622–0.701, a ratio of about 1.9–2.0 (Supplementary Table S10). The pure glasses show the same "
      "difference, within 0.002 of the composites. It is already present in the unloaded glass and at the peak, and "
      "it changes by only 0.4–1.5% over the post-peak history (Table 6). The asymmetry is therefore a property of the "
      "glass matrix that exists before loading."),

("P", "Part of this factor simply reflects the larger size of Ca atoms. Because the raw per-atom value is stress × "
      "volume, a larger atom has a larger value at equal local stress. The Voronoi volumes of the unloaded pure glasses "
      "give V_Ca/V_Zn = 1.23–1.24. Dividing each atom’s stress tensor by its own Voronoi volume lowers the Ca–Zn "
      "log-stress difference to 0.451, 0.431 and 0.417 before loading, and to 0.478, 0.460 and 0.447 at the peak "
      "(Table 6; Fig. 11b). The reduction equals ln(V_Ca/V_Zn) = 0.21, so atomic volume accounts for 31–34% of the "
      "raw difference. The remainder is a real excess. Ca atoms carry 1.5–1.6 times the volume-normalized von Mises "
      "stress of Zn atoms, most of it already present in the as-quenched glass, and the single-seed composites lie "
      "within 0.007 of their glasses. The load increment is also shared unequally. The volume-weighted mean change of "
      "axial stress per species, divided by the change of the cell average, is 1.02–1.05 for Mg, 0.77–0.82 for Zn and "
      "1.67–1.89 for Ca. This shows where the load goes, not why."),

("TBL", ("t14", 5,
         "Table 6.|Ca–Zn difference of the group means of ln σvM (Δln) in the three pure glasses, for the raw per-atom "
         "stress (stress × volume, as reported by LAMMPS) and the volume-normalized stress (stress divided by the atom’s "
         "Voronoi volume). The volume part is the difference of the two and equals ln(V_Ca/V_Zn). The last column gives "
         "the volume-normalized Δln of the single-seed composite built on the same matrix, at the nearest analyzed "
         "frame.")),

("FIG", (F + "/fig14/fig10_stress_asymmetry.png", 6.5,
         "Fig. 11.|(a) Mean absolute SHAP values of the descriptors of the volume-normalized stress model (model G) for "
         "the pure glass of xCa = 0.05 at 15% strain; the Voronoi volume and the own-type indicator lead. (b) Ca–Zn "
         "log-stress difference of the raw (stress × volume) and volume-normalized per-atom stress in the three pure "
         "glasses (solid) and the single-seed composites (hatched) at the peak; the dotted line is "
         "ln(V_Ca/V_Zn) = 0.21, the part of the raw difference due to atomic volume.")),

("P", "Ca and Zn atoms also occupy different local environments. At each composition’s peak strain the mean q6 was "
      "0.492 for Zn atoms at all three compositions and 0.312–0.317 for Ca atoms, about 35% lower (Fig. 12a); the "
      "distributions differ significantly (Mann–Whitney p < 0.001). Non-averaged q6 depends on the number of neighbors "
      "(Section 2.10), and the single 4.0 Å cutoff captures different parts of the first shell of the large Ca and "
      "the small Zn atoms. The q6 difference therefore cannot be attributed to packing order alone, and a comparison "
      "at matched neighbor number or with averaged parameters [@lechner] would be needed. A mean q6 close to the "
      "ideal HCP value (0.484) also does not by itself indicate HCP-like order around Zn, because non-crystalline "
      "clusters give values in a similar or higher range (0.663 for an ideal icosahedron)."),

("P", "Both species have fewer same-type neighbors than random mixing would give (αAA > 0). The deficit is 2.2–2.3 "
      "times larger for Ca (αCaCa = 0.511, 0.498 and 0.463 at xCa = 0.03, 0.05 and 0.07) than for Zn "
      "(αZnZn = 0.222, 0.220 and 0.211) (Fig. 12b), so Ca atoms are more strongly separated from one another than Zn "
      "atoms are."),

("FIG", (F + "/fig14/fig11_structure.png", 6.0,
         "Fig. 12.|(a) Mean non-averaged q6 of Zn and Ca atoms at the peak strain, with the ideal HCP and FCC values for "
         "reference; non-averaged q6 depends on the number of neighbors within the 4.0 Å cutoff (Section 2.10). "
         "(b) Same-type Warren–Cowley parameters of Zn and Ca; the deficit of Ca is 2.2–2.3 times that of Zn at every "
         "composition.")),

("P", "That the environments differ does not show that they account for the stress difference. Local descriptors "
      "were therefore added to the regression on the volume-normalized target, on the same atoms and split "
      "(Table 7). From a baseline of identity and coordination (R² = 0.092–0.133), bond-orientational order added at "
      "most 0.012 and chemical order at most 0.004. Mean bond length and coordination added 0.044–0.063, and the "
      "Voronoi volume alone added 0.097–0.136, more than all other structural descriptors together. With every "
      "descriptor R² reached 0.243–0.267, and the Voronoi volume carried the largest mean absolute SHAP value, above "
      "the own-type indicator (Fig. 11a). A structure-only model, with no identity or chemical-order information, "
      "reached R² = 0.150–0.183 and predicted 81–88% of the normalized Ca–Zn difference. The same ordering of "
      "descriptors holds for the raw target, where the structure-only model predicts 88–91% of a difference that still "
      "contains the volume factor (Supplementary Table S9)."),

("P", "Local spacing (volume and bond length), not bond-orientational or chemical order, therefore carries the "
      "structural information about per-atom stress that these descriptors provide. Static local structure still "
      "explains at most 27% of the per-atom variance, and less within a species (R² of 0.067–0.091 among Ca atoms and "
      "0.107–0.122 among Zn atoms). Because the split is random over atoms, these values are, if anything, "
      "optimistic."),

("TBL", ("t14", 6,
         "Table 7.|Volume-normalized regression in the three pure glasses (per-atom stress divided by the atom’s "
         "Voronoi volume) at zero strain, the stress peak and 15% strain. Held-out R² of models A–G (A: identity and "
         "3.25 Å counts; B: A + q4, q6; C: A + α; D: A + mean bond length, coordination; E: all descriptors; F: A + "
         "Voronoi volume; G: E + Voronoi volume), R² of the structure-only model, the share of the normalized Ca–Zn "
         "log-stress difference it predicts without and with the Voronoi volume, and the mean absolute SHAP value of "
         "the Voronoi volume in model G. The equivalent regression on the raw target is in Supplementary Table S9.")),

("H2", "3.7. Strain-rate and temperature sensitivity"),

("P", "Peak strength fell linearly with temperature, from 1.315 GPa at 250 K to 0.885 GPa at 400 K, a slope of "
      "−2.86 MPa K⁻¹ or −0.245% per K relative to 300 K (Fig. 13b; Table 8). The 310 K case lies on this line, and "
      "the repeat of the 300 K case reproduces the baseline run to 0.03%. Over the same range Young’s modulus fell "
      "from 31.9 to 26.5 GPa and the peak strain from 7.4% to 6.4%, while the softening ratio rose from 20.9% to "
      "24.4%. The thermal pre-stress (Section 2.4) leaves the peak stress unchanged but shifts the peak strain by up "
      "to 0.26%, so the strengths are unaffected, whereas the peak strains and moduli carry this uncertainty. The "
      "whole range lies well below Tg, and the result is for one structure."),

("P", "Peak strength rose monotonically with strain rate over a 20-fold range, from 0.864 to 1.286 GPa between "
      "1×10⁸ and 2×10⁹ s⁻¹ (+48.8%; Fig. 13a). The end points give a logarithmic rate sensitivity "
      "d ln σmax/d ln ε̇ of about 0.13. The increase per decade rose from 0.28 to 0.35 to 0.40 GPa between successive "
      "pairs, so the dependence is not log-linear over this range. With four points we report this as an observation "
      "and fit no functional form. The peak strain and Young’s modulus also rose with rate. The modulus change may "
      "partly reflect finite-rate loading over the 0.1–1% fitting window and is not interpreted."),

("P", "To test whether this is specific to the baseline, the step from 10⁹ to 2×10⁹ s⁻¹ was repeated on the relaxed "
      "structures of Compositions 6, 4 and 3. Together with Composition 1 these span xZn = 0.25–0.40 and "
      "xCa = 0.05–0.07 and differ by 13% in strength at 10⁹ s⁻¹. Doubling the rate raised the peak strength by a "
      "factor of 1.104, 1.102, 1.109 and 1.099 at Compositions 1, 6, 4 and 3, a mean of 1.103 with a standard "
      "deviation of 0.004 (Fig. 13c). Strain-rate hardening is therefore not specific to the baseline composition. "
      "Each pair reuses one relaxed structure and one velocity seed, so the spread reflects differences between "
      "compositions and not repeat-run scatter. No discrete shear band was detected at any rate or temperature "
      "(Supplementary Table S6)."),

("FIG", (F + "/fig15/fig12_rate_temperature.png", 6.5,
         "Fig. 13.|(a) Peak strength versus strain rate at 300 K for the four compositions loaded at 2×10⁹ s⁻¹ (log "
         "axis); points are offset slightly along the rate axis so that coinciding values remain visible. (b) Peak "
         "strength (filled) and flow stress at 15% strain (hollow) versus temperature at 10⁹ s⁻¹ for Composition 1; the "
         "line is the linear fit of the peak strength (−2.86 MPa K⁻¹), and the cross marks the 310 K case. "
         "(c) Increase in peak strength for a doubling of the strain rate (10⁹ to 2×10⁹ s⁻¹) at four compositions.")),

("TBL", ("t14", 7,
         "Table 8.|Temperature and strain-rate cases on the Composition 1 structure, and the 2×10⁹ s⁻¹ cases on "
         "Compositions 6, 4 and 3 with their 10⁹ s⁻¹ values for comparison.")),

("H2", "3.8. Gaussian process surrogate"),

("P", "In leave-one-out cross-validation of the 16-point model the standardized residuals had a standard deviation of "
      "1.08 (ideal 1), and 14 of the 16 held-out points fell inside their 95% intervals. The 310 K case has no "
      "informative interval, since with it held out no other training point differs in temperature. Of the remaining "
      "15, 13 lie inside their intervals (Fig. 14a); the exceptions are Composition 10 (z = −2.45) and the 18.49% case "
      "(z = −2.14)."),

("P", "The 13-point model, whose training set stopped at 15% volume fraction, predicted the two saturated cases "
      "before they were run. It gave 0.926 GPa (σ = 0.011 GPa) at 18.91%, against a simulated 0.928 GPa (error "
      "+0.2%), and 0.931 GPa (σ = 0.010 GPa) at 18.49%, against 0.917 GPa (−1.5%). Linear extrapolation of the "
      "5–15% trend in the same matrix gives 0.923 and 0.928 GPa (errors +0.6% and −1.2%), so on these two tests the "
      "surrogate did no better than the simplest trend model."),

("P", "The 23 second-campaign simulations then served as held-out data for frozen models fitted to earlier data "
      "without tuning, which makes this a retrospective test. Seventeen of the 23 lie within the range of inputs seen "
      "in training and have informative intervals (Fig. 14b; Table 9; Supplementary Table S11). Sixteen of these 17 "
      "simulated values lie inside the 95% interval, with a mean absolute error of 0.99%. The exception is the pure "
      "glass of xCa = 0.05, underpredicted by 3.9%. All three glasses are underpredicted, by 1.8–3.9%, because their "
      "volume fraction of zero lies outside the 0.16–19% range of the composites in the training set."),

("P", "The temperature series is not a test of the model. In the 25-point training set temperature took only the "
      "values 300 and 310 K, so 250, 350 and 400 K lie far from the training data in scaled units and the predictions "
      "are controlled by the prior. They were 1.091, 1.101 and 1.109 GPa (σ ≈ 0.40 GPa), against simulated values of "
      "1.315, 1.024 and 0.885 GPa. The simulated values fall inside these intervals because the intervals are wide, "
      "not because the model has predictive skill."),

("FIG", (F + "/fig15/fig13_surrogate.png", 6.5,
         "Fig. 14.|(a) Leave-one-out predictions (±1.96σ) versus simulated peak strength for the fifteen configurations "
         "with informative intervals; open points fall outside their interval. (b) Predictions (±1.96σ) of the frozen "
         "surrogates versus simulated peak strength for the 17 second-campaign simulations with informative intervals "
         "and the two high-volume-fraction cases predicted before they were run. Prior-dominated temperature "
         "predictions are listed in Table 9 and Supplementary Table S11 and are not plotted.")),

("TBL", ("t14", 8,
         "Table 9.|Held-out tests of the surrogate. Predictions come from models fitted only to data that existed before "
         "each simulation; the 23 second-campaign simulations are retrospective tests (Section 2.12). Cases whose inputs "
         "were never varied in training have intervals that reflect the prior only.")),

]
