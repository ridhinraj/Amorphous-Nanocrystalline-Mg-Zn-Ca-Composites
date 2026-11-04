# -*- coding: utf-8 -*-
# Draft 15, part A: introduction and methods
F = "/tmp/w"

BLOCKS = [

("H1", "1. Introduction"),

("P", "Magnesium alloys are studied as biodegradable materials for temporary orthopedic implants. Magnesium is "
      "naturally present in the human body, corrodes gradually in vivo, which may avoid a second operation to remove "
      "the implant, and has an elastic modulus closer to that of bone than permanent implant metals such as titanium "
      "or stainless steel [@bairagi]. Within Mg alloys, Mg-Zn-Ca has attracted particular interest. Zinc and calcium "
      "are essential trace elements, and some compositions in this ternary combine good corrosion resistance with "
      "glass-forming ability under rapid solidification [@vinogradov]. Mg-Zn-Ca glasses with more than about 28 at.% "
      "Zn, including Mg₆₀Zn₃₅Ca₅, showed no clinically observable hydrogen evolution in vivo [@zberg09], and "
      "low-alloyed Mg-Zn-Ca has degraded on timescales compatible with bone healing in osteoporotic, old and juvenile "
      "rats [@sommer]."),

("P", "Mg-based alloys, crystalline and amorphous alike, are brittle compared with what a load-bearing implant "
      "requires. Metallic glasses deform at room temperature by localized shear banding, which limits their tensile "
      "ductility [@schuh; @cheng]. Amorphous-nanocrystalline (AN) microstructures, in which nanoscale crystals are "
      "dispersed in a glassy matrix, have been proposed as a way to improve damage tolerance, drawing on shear-band "
      "control mechanisms established in other metallic glass composites [@li26]. Amorphous and partially crystalline "
      "Mg-Zn-Ca structures of this kind have been produced by laser powder bed fusion, mainly for the corrosion "
      "advantage of the amorphous phase [@wang22; @chen23]. In the present work Mg-Zn-Ca serves as a model system for "
      "the mechanics of such composites. Corrosion and biocompatibility are outside its scope."),

("P", "The atomistic mechanisms behind the mechanical response of these composites are largely unexplored. Molecular "
      "dynamics (MD) studies of Mg-Zn-Ca glasses have addressed their atomic structure [@gulenko], the effect of "
      "alloying additions on the monolithic glass [@sun20] and, from first principles, their structural and "
      "electronic properties [@erkartal]. In other systems, MD has shown how crystalline nanoprecipitates affect "
      "shear-band propagation in Cu-Zr-based glasses [@brink16], how interfaces control creep in metallic glass "
      "composites [@kalcher], and how deformation changes from glass-like to grain-boundary-mediated as the amorphous "
      "fraction vanishes [@brink18; @albe13]. How the response of a Mg-Zn-Ca AN composite depends on alloy "
      "composition, nanocrystal volume fraction and the amorphous/crystalline interphase has not been reported. "
      "Separately, segregation of Ca and Zn to grain and twin boundaries affects boundary strength in crystalline Mg "
      "alloys [@basu; @jang22]. Whether Ca and Zn atoms play distinguishable roles inside an amorphous matrix has not "
      "been examined."),

("P", "Each simulated configuration requires equilibrating and loading a cell of several hundred thousand atoms, which "
      "makes an exhaustive sweep over composition, volume fraction, temperature and strain rate expensive. Surrogate "
      "models trained on a modest set of simulations, with their uncertainty used to choose further runs [@xie23], and "
      "explainable machine learning that extracts atomic-scale correlations from the same simulations "
      "[@lu23; @zhong22] offer ways to get more out of such a campaign."),

("P", "We tested two hypotheses: (i) that composition and nanocrystal volume fraction control the strength of Mg-Zn-Ca "
      "AN composites mainly by changing the nucleation and arrest of a discrete shear band, and (ii) that the "
      "amorphous/crystalline interphase is a preferred site of plastic rearrangement. Using a second-nearest-neighbor "
      "modified embedded-atom method (2NN MEAM) potential for Mg-Zn-Ca [@jang19; @jang18], we built glass and AN cells "
      "of about 380,000 atoms (20 nm) and loaded them in tension in 48 simulations. They cover ten compositions "
      "centered on Mg₆₀Zn₃₅Ca₅ (xZn = 0.25–0.40, xCa = 0.03–0.07), pure-glass controls at three Ca contents, "
      "nanocrystal volume fractions of 5–19% in three matrices with packing replicates, and series in seed radius, "
      "strain rate and temperature. We report peak strength, modulus and post-peak softening. We test for localized "
      "deformation with per-atom stress and non-affine displacement [@falk], test directly whether the interphase "
      "concentrates rearrangement, compare the local stress and environment of Ca and Zn atoms, and check a Gaussian "
      "process surrogate against held-out runs and simple extrapolation."),

("P", "The main new results are pure-glass controls that place the calcium effect in the glass matrix, a replicated "
      "volume-fraction series in three matrices, a direct test of the interphase hypothesis over all 236 seeds of two "
      "composites, and a volume-normalized comparison of the stress carried by Ca and Zn atoms."),

("H1", "2. Computational methods"),

("P", "The workflow is summarized in Fig. 1. The simulations were run in two campaigns. The first (25 simulations) "
      "covered the ten compositions, the initial volume-fraction and seed-geometry cases and three strain rates. The "
      "second (23 simulations) added the pure-glass controls, the xCa = 0.03 volume-fraction series, the packing "
      "replicates, the temperature series and two further compositions at 2×10⁹ s⁻¹, and part of it was chosen with "
      "the surrogate of Section 2.12. All 48 used the protocols below. Validation and implementation details that are "
      "not needed to follow the results are given in the Supplementary Material."),

("FIG", (F + "/fig15/f1/fig01.png", 6.5,
         "Fig. 1.|Workflow of the study. (a) Model construction: interatomic potential, melt–quench matrices, composite "
         "construction, relaxation and tension. (b) The 48 tensile simulations, grouped by the variable changed. "
         "(c) Analyses applied to the simulations. The dashed arrow marks follow-up runs ranked by the Gaussian process "
         "surrogate. Section numbers refer to this paper.")),

("H2", "2.1. Interatomic potential"),

("P", "All simulations used LAMMPS [@lammps] (version 22 July 2025, update 4) with the 2NN MEAM potential for Mg-Zn-Ca "
      "of Jang et al. [@jang19], which builds on the Mg-Zn potential of the same formalism [@jang18]. The parameter "
      "files were taken from the NIST Interatomic Potentials Repository (potential 2019--Jang-H-S-Seol-D-Lee-B-J--"
      "Mg-Zn-Ca, implementation 2019--Jang-H-S--Mg-Zn-Ca--LAMMPS--ipr1) [@becker; @hale]. A short-range "
      "Ziegler–Biersack–Littmark (ZBL) correction [@ziegler] was added through a hybrid/overlay pair style. It acts "
      "only below 2 Å, where it prevents unphysically large forces at close contact during composite construction, "
      "and does not affect equilibrium behavior."),

("P", "The potential was checked (Section 3.1) against the lattice geometry, cohesive energy, single-crystal elastic "
      "constants and basal generalized stacking-fault energy (GSFE) of pure Mg, and against the liquid mixing enthalpy "
      "of the ternary. The GSFE was computed on a 4,000-atom slab of relaxed HCP Mg. The upper half was displaced "
      "rigidly relative to the lower half, allowing only out-of-plane relaxation, along the straight ⟨112̅0⟩ Burgers "
      "vector and, separately, along the ⟨101̅0⟩ partial-dislocation path on which first-principles values are quoted. "
      "Elastic constants were obtained from finite differences of stress with respect to strain, with positions "
      "re-minimized in each strained state. The mixing enthalpy was ΔHmix = E_{ternary} − (xMg E_{Mg} + xZn E_{Zn} + "
      "xCa E_{Ca}), from the mean potential energy per atom of separately equilibrated 4,000-atom liquids (NPT, zero "
      "pressure, 1000 K). Each liquid was prepared with one melt-and-cool protocol, confirmed to be liquid from its "
      "mean-squared displacement and run with three velocity seeds. At zero pressure and equal temperature the kinetic "
      "and pressure–volume terms cancel, so the potential-energy difference gives the mixing enthalpy."),

("H2", "2.2. Melt–quench protocol and composition matrix"),

("P", "Each composition was built as a 389,344-atom cubic cell, 4 × 46³ sites of a face-centered cubic lattice with "
      "a lattice constant of 4.5 Å (207 Å per side), populated at random with the target element ratios. All "
      "melt–quench stages used a 0.5 fs timestep. The system was heated from 300 to 1500 K (100,000 steps), held at "
      "1500 K (200,000 steps) and quenched to 300 K over 2,400,000 steps (1.2 ns, a nominal 10¹² K/s) under isotropic "
      "NPT coupling. After a 300 K settling stage (100,000 steps), each structure was reheated to 600 K (200,000 "
      "steps), held for 500 ps, cooled again to 300 K at 10¹² K/s (600,000 steps) and given a final 300 K stage "
      "(200,000 steps). All NPT stages used Nosé–Hoover coupling [@nose; @hoover] with damping constants of 0.1 ps "
      "for temperature and 5.0 ps for pressure. After quenching, the cell edge was 202.7–207.0 Å, depending on "
      "composition."),

("P", "Ten compositions were generated in this way, one quench each, centered on Mg₆₀Zn₃₅Ca₅ (Table 1). "
      "Composition 1 (Mg₆₀Zn₃₅Ca₅) is the baseline for the seed-radius, strain-rate and temperature series. "
      "Compositions 1 and 9 were resumed from a checkpoint during the first quench, so their effective quench rates "
      "were 1.14×10¹² and 1.20×10¹² K/s."),

("TBL", ("t14", 0,
         "Table 1.|Composition matrix and the single seed used for each composition in the composition sweep. Seed "
         "volume fractions are computed on the quenched cell of each composition.")),

("H2", "2.3. Composite construction and pure-glass controls"),

("P", "Composites were built by removing a spherical region from the quenched matrix and inserting a hexagonal "
      "close-packed (HCP) α-Mg crystal of the same radius. The seed was cut from an HCP lattice with a = 3.20 Å and "
      "c = 5.21 Å (c/a = 1.628), within 0.3% of the values of this potential (Section 3.1); the small mismatch relaxes "
      "during minimization and NPT equilibration. No rotation was applied, and every seed in every composite therefore has its "
      "a-axis ⟨2̅110⟩ along the tensile axis (x) and its c-axis [0001] along z. Matrix atoms within 2.0 Å of any "
      "inserted crystal atom were removed, and composite cells therefore differ slightly in atom count from their "
      "matrices (377,736–389,344 atoms; Supplementary Table S1)."),

("P", "For the ten-composition sweep, one seed was inserted at the cell center, with a radius of 15–35 Å that rises "
      "with Ca content (Table 1; seed volume fraction 0.16–2.02%). This sweep therefore varies composition and seed "
      "geometry together, a confound that the pure-glass controls below resolve. For the volume-fraction series, 12 Å "
      "seeds were placed by rejection sampling at target volume fractions of 5%, 10% and 15% in each of three "
      "matrices, reusing each matrix’s quench. These are Compositions 5, 1 and 6, with xCa = 0.03, 0.05 and 0.07 at "
      "xZn = 0.35. Seed centers were kept at least 2R + 4 Å apart, leaving a gap of at least 4 Å between seed "
      "surfaces, and at least R + 4 Å from every cell face, so that no seed crosses a periodic boundary and the "
      "crystals never form a connected path. Targets of 20% and 25% in the xCa = 0.05 matrix placed only 223 of 236 "
      "and 218 of 295 seeds, even with the attempt budget raised from 2×10⁵ to 2×10⁶, and are reported at their "
      "achieved values of 18.91% and 18.49%. Seeds of 15 and 20 Å at a 10% target in the same matrix complete a "
      "seed-radius series with the 12 Å case. All volume fractions are computed on the quenched cell volume. The 5% "
      "and 15% composites of the xCa = 0.05 matrix, with 59 and 177 seeds, are called B5 and B15 below. Structures "
      "were inspected with OVITO [@stukowski] (Fig. 2)."),

("P", "The three pure-glass controls are the matrices of Compositions 5, 1 and 6 without seeds, copied, relaxed and "
      "loaded exactly as the composites. A control therefore differs from its composites only by the absence of "
      "seeds. To measure the scatter caused by seed placement alone, three further packings with different random seed "
      "positions were built at 10% volume fraction in each matrix, giving four per matrix. All packings in a matrix "
      "share that matrix’s quench, so the replicates sample seed placement and not independent quenches. "
      "Perfect-crystal references are defect-free HCP Mg cells of 339,264 atoms, relaxed by energy minimization with "
      "box relaxation and 20,000 steps of NPT at 300 K, then loaded along ⟨2̅110⟩ and [0001]."),

("FIG", (F + "/b15/fig02_placeholder.png", 6.3,
         "Fig. 2.|Atomistic views of the simulation cells, rendered with OVITO 3.17 [@stukowski]. (a) B15 (xZn = 0.35, "
         "xCa = 0.05; 177 seeds of 12 Å radius, 15% volume fraction) at zero strain, showing only the atoms identified "
         "as HCP by polyhedral template matching [@larsen] (RMSD cutoff 0.2). These atoms form 177 separate clusters, "
         "one per seed. (b) A thin slab through Composition 1 at 7.0% strain, colored by incremental D²min over the "
         "preceding 0.5% strain interval; atoms above the calibrated threshold (Section 2.8) are shown in magenta. The "
         "tensile axis is horizontal. (c) A 4 Å slab through one seed of B5 at zero strain, showing the seed interior, "
         "the 3 Å interphase shell, the 4.5 Å buffer and the bulk reference used in the interphase analysis "
         "(Section 2.9).")),

("H2", "2.4. Relaxation and mechanical loading"),

("P", "Each structure was relaxed by energy minimization and 50,000 steps of NPT equilibration at 300 K. All stages "
      "from here on used a 1 fs timestep. Tensile loading consisted of a further 20,000-step equilibration with "
      "pressure coupling in the two lateral directions, followed by 150,000 steps of deformation along x with the "
      "LAMMPS fix deform command at an engineering strain rate of 10⁹ s⁻¹, to 15% strain. The barostat held both "
      "lateral stresses at zero, so the stress state approximates free uniaxial tension, and all boundaries remained "
      "periodic. The macroscopic stress was recorded every 100 steps and per-atom quantities every 5,000 steps (0.5% "
      "strain). All runs used the same velocity seed for loading."),

("P", "In the temperature series the relaxed 300 K structure of Composition 1 was loaded at 250, 300, 350 and 400 K, "
      "with the thermostat and initial velocities set to the test temperature during the 20,000-step equilibration. "
      "A 310 K case had been run earlier with the same protocol. In the strain-rate series the same structure was "
      "loaded at 1×10⁸, 5×10⁸ and 2×10⁹ s⁻¹, with the number of loading steps scaled inversely with rate so that each "
      "run reached 15% strain. The 2×10⁹ s⁻¹ rate was repeated on the relaxed structures of Compositions 6, 4 and 3. "
      "These series reuse an already relaxed structure and do not sample thermal history independently. The 300 K "
      "case repeats the baseline run with identical inputs and serves as a reproducibility check. The axial box "
      "length is fixed during equilibration, and a structure equilibrated away from 300 K therefore develops a thermal "
      "axial pre-stress. This shifts the stress–strain curve along itself. The peak stress is unchanged, but the "
      "strain at the peak moves by up to 0.26%."),

("H2", "2.5. Mechanical properties"),

("P", "Stress–strain curves were smoothed with a cubic Savitzky–Golay filter [@savgol] over a fixed strain width of "
      "0.01, so that runs at different rates are smoothed identically. The peak strength σmax is the maximum of the "
      "smoothed curve, and the peak strain εmax is the strain at that maximum. In MD of glass–crystal nanocomposites "
      "the peak stress is the usual strength measure because it is sensitive to the state of the glass [@brink18]. We "
      "call it peak strength rather than yield strength to avoid confusion with the 0.2% offset yield strength. "
      "Young’s modulus E is the slope of a least-squares line between 0.1% and 1.0% strain. The flow stress σ15 is "
      "the value at 15% strain of a linear fit to the final 0.01 of strain, and the post-peak softening ratio is "
      "(σmax − σ15)/σmax. For the perfect crystals the same definition gives an ideal strength. All properties were "
      "extracted by one script applied identically to all 48 curves (Supplementary Table S1). The volume-fraction "
      "slope is the least-squares slope, constrained through the origin, of the relative change from the matched pure "
      "glass against achieved volume fraction for the nine composites at 5–15%."),

("H2", "2.6. Structure and thermal history of the matrices"),

("P", "The radial distribution function g(r) of each quenched matrix was computed with a periodic implementation (bin "
      "width 0.05 Å, up to r = 10 Å), and the density from the cell volume, atom count and composition-weighted mean "
      "atomic mass. The thermodynamic output of each melt–quench, recorded every 0.5 ps, was analyzed with one script "
      "for all ten compositions. For each cooling ramp, the potential energy per atom and the specific volume were "
      "fitted against temperature with a continuous two-line model, and the breakpoint, searched on a 1 K grid, was "
      "taken as Tg. Because the liquid branch is curved, the first quench was fitted with upper window bounds of 700, "
      "800 and 1000 K (lower bound 300 K), and the mean and range are reported. The reheating and second cooling, "
      "which span 300–600 K, were fitted over 300–560, 300–580 and 300–600 K."),

("P", "The fictive temperature of the 300 K glass is Tf = (P(T0) − b_{g}T0 − a_{l})/(b_{l} − b_{g}), where P(T0) is "
      "the property at T0 = 300 K (mean over the last 50 ps of the 300 K stage), b_{g} is the slope of the glass line "
      "and a_{l} + b_{l}T is the liquid line, evaluated for each fitting window and averaged. Effective quench rates "
      "were taken from the slope of temperature against time between 400 and 1200 K. Tg and Tf obtained in this way "
      "are operational quantities of the simulated cooling rate and fitting window, not calorimetric values."),

("H2", "2.7. Per-atom stress and detection of localized deformation"),

("P", "Per-atom stress tensors were computed with the LAMMPS compute stress/atom command, which returns the product of "
      "the stress tensor and a per-atom volume [@thompson09]. LAMMPS does not divide by any volume, so these raw values "
      "have units of stress × volume. The localization test and the first regression model (Section 2.11) use the von "
      "Mises stress formed from these raw tensors. Dividing every atom by the same mean volume V/N would only shift "
      "all logarithmic values equally and would not change a comparison between species. Normalization by each atom’s "
      "own Voronoi volume is described in Section 2.11."),

("P", "To test for a discrete high-stress band, an absolute von Mises threshold was calibrated for each trajectory "
      "from its undeformed frame (mean plus four standard deviations), and atoms above it in each recorded frame were "
      "clustered with DBSCAN [@ester; @sklearn]. Clustering used wrapped coordinates without periodic images, so a "
      "region that crosses a cell face is split in two. A percentile threshold was not used because it flags a fixed "
      "fraction of atoms in every frame and reported spurious clusters even in undeformed cells. The test was applied "
      "to all 46 glass and composite simulations. Strain localization is more often identified from the atomic shear "
      "strain [@shimizu]. That measure was not computed, and the stress-based test was not checked against a cell "
      "known to form a shear band (Section 4.3)."),

("H2", "2.8. Non-affine displacement"),

("P", "Per-atom stress cannot distinguish elastic stress concentration from irreversible rearrangement, so the "
      "non-affine squared displacement D²min [@falk] was also computed for every atom at every recorded frame. For "
      "atom i with n_{i} neighbors j within 4.5 Å in the reference frame, the best-fit local affine transformation "
      "F_{i} was found by least squares, and"),

("EQ", "(1)"),

("P", "where Δr_{j} are the vectors from atom i to its neighbors in the current (cur) and reference (ref) frames, and "
      "F = XY⁻¹ with X = Σ_{j} Δr_{j}^{cur} ⊗ Δr_{j}^{ref} and Y = Σ_{j} Δr_{j}^{ref} ⊗ Δr_{j}^{ref}. Unlike the "
      "original definition [@falk], Eq. (1) is divided by n_{i}, and at least four neighbors were required. The 4.5 Å "
      "cutoff spans the first-neighbor distance of every species pair, including Ca–Ca (3.83 Å), and gives about "
      "17–18 neighbors per atom at the number densities of these glasses (0.044–0.047 Å⁻³). The neighbor cutoffs of "
      "all analyses are listed in Supplementary Table S2."),

("P", "Two reference conventions were used. Cumulative D²min compares each frame with the undeformed frame, and "
      "incremental D²min compares each frame with the preceding recorded frame, an interval of 5,000 steps or 0.5% "
      "strain. For each trajectory an absolute threshold (mean plus four standard deviations) was calibrated from the "
      "first interval (0–0.5% strain), which is elastic since the peak strains are 5.8–8.0%. Atoms above it were "
      "clustered with DBSCAN (eps = 6.0 Å, min_samples = 15) to count spatially distinct, simultaneously active "
      "rearrangement regions. The sensitivity of the counts to the cutoff, threshold multiplier and clustering "
      "parameters was examined by varying one setting at a time (11 settings; Supplementary Table S7)."),

("H2", "2.9. Interphase analysis"),

("P", "To test whether rearrangement concentrates at the amorphous/crystalline interphase, incremental D²min was "
      "compared between a shell |r| ≤ 3 Å around each seed and a bulk region farther than 7.5 Å from the surface of "
      "every seed. Here r = |x − x_{c}| − R is the signed distance from the surface of a seed with center x_{c} and "
      "radius R. The 4.5 Å gap between shell and bulk equals the D²min cutoff, so an atom’s score is not influenced by "
      "neighbors in the adjacent region. In a synthetic test without this gap, an engineered interphase-only "
      "rearrangement produced a spurious signal in the crystal interior. With the gap, the engineered rearrangement "
      "was recovered as a 25-fold interphase-to-bulk ratio."),

("P", "Seed centers and radii were found from the undeformed frame. Mg atoms with q6 > 0.40 and a local Mg fraction "
      "above 0.95 (4.0 Å cutoff) were flagged and clustered with DBSCAN (eps = 6.0 Å, min_samples = 15). On a "
      "synthetic composite this recovered the seed center within 1.1 Å and the radius within 1.4 Å, and in B5 and "
      "B15 it recovered exactly the 59 and 177 seeds placed. The recovered radius (9.7 Å) underestimates the 12 Å "
      "build radius because the surface atoms of a seed fail the selection criteria, so shells were placed at both "
      "radii. Every seed of B5 and B15 was analyzed at 2%, 4%, 7% and 10% strain. Per-seed ratios were tested against "
      "1 (t-test and Wilcoxon signed-rank test) and compared between strains with paired tests. The seeds of one "
      "composite come from a single trajectory, so all p-values are conditional on that realization."),

("H2", "2.10. Bond-orientational and chemical short-range order"),

("P", "Steinhardt bond-orientational order parameters q4 and q6 [@steinhardt] were computed for every Zn and Ca atom "
      "in their standard, non-averaged form, using all neighbors within 4.0 Å regardless of species:"),

("EQ", "(2)"),
("EQ", "(3)"),

("P", "where N_{b} is the number of neighbors and (θ_{j}, φ_{j}) are the bond angles. A 4.0 Å cutoff was used because "
      "at 3.25 Å, 37% of Ca atoms in Composition 1 had fewer than four neighbors (the Ca–Ca first-neighbor distance "
      "is 3.83 Å). The implementation reproduces the ideal-lattice values (FCC q6 = 0.575, HCP q6 = 0.484). For "
      "uncorrelated points ql scales roughly as N_{b}^{−1/2}, so non-averaged values depend on the coordination "
      "number. Averaged parameters [@lechner], which are less sensitive to this, were not computed. The "
      "Warren–Cowley parameter [@cowley; @warren] for species A and B,"),

("EQ", "(4)"),

("P", "was computed with the same 4.0 Å cutoff. Here PAB is the mean fraction of the neighbors of A atoms that are of "
      "type B, and xB is the atomic fraction of B. A value of αAB = 0 indicates random mixing, negative values favor "
      "A–B pairing, and positive values mean fewer A–B pairs than random mixing would give. Same-type parameters "
      "(αCaCa, αZnZn) are reported, and a per-atom version of αAA was used to compare the distributions of Ca and Zn "
      "atoms (Mann–Whitney U test). The finite-sample bias of αCaCa was checked directly (+0.23 at 245 Ca atoms, "
      "+0.01 at 3,971) and is negligible at the 12,000–27,000 Ca atoms per configuration analyzed here. One cutoff "
      "was used for all pairs rather than pair-specific cutoffs at the first minimum of each partial g(r)."),

("H2", "2.11. Volume normalization and per-atom stress regression"),

("P", "For every Zn and Ca atom, the numbers of Mg, Zn and Ca neighbors within 3.25 Å and an own-type indicator "
      "(Ca = 1, Zn = 0) were used to predict the natural logarithm of the atom’s von Mises stress by gradient-boosted "
      "regression [@friedman; @sklearn] with default hyperparameters and a 75/25 random train/test split (model A). "
      "Feature attributions were computed on the test set with SHAP [@shap], as in explainable machine learning for "
      "materials [@lu23; @zhong22]. The 3.25 Å cutoff lies near the Mg–Mg first peak (3.175 Å) and below the first "
      "minimum of the total g(r) (about 4.0 Å) and counts only the nearest part of the first shell. This model was "
      "evaluated at the final frame (15% strain) of thirteen configurations, and at the peak strain, 10% and 15% "
      "strain for Compositions 5, 1 and 6."),

("P", "Because the raw per-atom value is stress × volume, its von Mises value mixes stress with atomic volume, which "
      "differs between species. The Voronoi cell volume of every atom was therefore computed on the periodic cell "
      "with a Qhull [@qhull] tessellation through SciPy [@scipy], and each atom’s stress tensor was divided by its own "
      "Voronoi volume before the von Mises stress was formed. The volume-weighted mean of this normalized stress "
      "reproduced the thermodynamic stress of the run to the printed 0.001 GPa at every frame. The regression was "
      "repeated on this target in the three pure glasses, on the same atoms and split, with descriptors added in "
      "groups. Model B adds q4 and q6, C the per-atom Warren–Cowley α, and D the mean neighbor distance and "
      "coordination number within 4.0 Å. Model E has all of these plus the neighbor counts by species within 4.0 Å, "
      "and models F and G add the Voronoi volume to A and E. A structure-only model, with no own-type or "
      "chemical-order information, was fitted to ask how much of the Ca–Zn log-stress difference local structure "
      "alone predicts, using three-fold out-of-fold predictions. The split is random across atoms within one "
      "configuration, so spatially correlated neighbors fall in both sets and the held-out R² values are optimistic."),

("H2", "2.12. Gaussian process surrogate"),

("P", "A Gaussian process regression model [@rasmussen; @sklearn] was trained on peak strength with standardized "
      "composition, volume-fraction, seed-radius, temperature and strain-rate inputs. The kernel was a constant times "
      "an isotropic Matérn kernel (ν = 2.5) plus a white-noise term, with normalized targets, five optimizer restarts "
      "and hyperparameters set by maximizing the log marginal likelihood. Calibration was assessed by leave-one-out "
      "cross-validation of a 16-point model, and untested configurations were ranked by predictive standard "
      "deviation. No refit-and-reselect cycle was performed [@xie23]."),

("P", "Each model was trained only on data that existed before the simulations it was asked to predict. A 13-point "
      "model was used for the two high-volume-fraction cases and the 310 K case. Its predictions for the two "
      "high-volume-fraction cases were made before those runs, and were also compared with linear extrapolation of the "
      "5–15% trend in the same matrix. A 16-point model was used for the pure glasses, the xCa = 0.03 series and the "
      "packing replicates, and a 25-point model for the temperature series and the 2×10⁹ s⁻¹ cases at Compositions 3 "
      "and 4. The second-campaign predictions were computed from these frozen models after the runs had finished, so "
      "they are retrospective tests. An input that takes a single value in training has zero variance and is left "
      "unscaled, and predictions for new values of it are controlled by the prior. Such cases are reported but not "
      "counted as tests. An anisotropic (automatic relevance determination) kernel and a logarithmic rate input would "
      "avoid this but were not used."),

("H2", "2.13. Robustness checks and software"),

("P", "Four analyses test how far the conclusions depend on analysis choices: the one-at-a-time scan of the D²min "
      "cluster settings (Section 2.8), the regression with added structural descriptors (Section 2.11), the interphase "
      "test over every seed of two composites at four strains (Section 2.9) and the volume normalization of per-atom "
      "stress. The scans used a vectorized implementation, checked against the per-atom loop implementation used for "
      "the baseline results (Supplementary Section S4). Analyses were scripted in Python with NumPy [@numpy], SciPy "
      "[@scipy] and scikit-learn [@sklearn], and figures were drawn with Matplotlib [@mpl]."),

]
