# 2x2 G4RW

Hadronic secondary-interaction (SI) systematics for the DUNE ND-LAr
Demonstrator (2x2), computed with Geant4Reweight (G4RW). This repo holds
the patches, configs, and analysis code behind the technote in `docs/`,
organized so each result can be reproduced from the files here plus the
NERSC paths and upstream sources listed per file.

## Layout

Pipeline order, upstream to downstream:

1. `edepsim_patch/` : EDepSim trajectory patch
2. `fracs_files/` : G4Cascade exclusive-channel fraction files
3. `geant4reweight_patch/` : Geant4Reweight source patches
4. `weight_production/` : per-track G4RW weight executables
5. `validation/` : cross-checks on the above (matching, trajectory counts, proc/subproc tags)
6. `fitting/` : weight-vs-length parameterization
7. `analysis_plots/` : data/MC comparison plots with systematic bands
8. `common/` : shared ROOT style file
9. `install_notes/` : build/run instructions per component
10. `docs/` : the technote (PDF and LaTeX source)

Most ROOT macros load `common/protoDUNEStyle.C` via
`gROOT->LoadMacro("protoDUNEStyle.C")`, which expects the file in the
same working directory. Copy or symlink it in before running a macro
from elsewhere in the repo.

---

## edepsim_patch/

Patch to EDepSim's trajectory pruning. The stock `SelectTrajectoryPoints`
keeps only geometric boundary crossings and endpoints (about 2 points
per track), too coarse for G4RW's step-by-step weight calculation.

- **EDepSimPersistencyManager.cc** : bypasses pruning for charged pions (±211),
  protons (2212), neutrons (2112), and charged kaons (±321), keeping every Geant4
  transport step. Other species are unaffected. Kaons carry the patch as
  groundwork only; no reweighting is implemented for them. Replaces
  `<edep-sim source>/src/EDepSimPersistencyManager.cc`.
  Source: `/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW/patch/EDepSimPersistencyManager.cc`
- **run_g4rw_edepsim_batched.sh** : runs the patched edep-sim build over
  MiniRun5 GENIE files in groups of 10, hadds each group, cleans up
  intermediates. Uses the withRock geometry and genie.nu input only
  (`patched_withRock_genieNu`, Table III Sample 5 in the technote).
  Usage: `bash run_g4rw_edepsim_batched.sh <start_group> <end_group> [genie_offset]`.
  Source: `/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW/batch_scripts/run_g4rw_edepsim_batched.sh`

## fracs_files/

Exclusive-channel fraction files: the fraction of the total inelastic
cross section carried by each channel vs. momentum (technote Section
III B, Eq. 7). Proton and neutron files were generated in 3 momentum
chunks (G4Cascade's memory use scales with NDivisions, and the full
0 to 9000 MeV/c range exceeds NERSC's memory budget in one run) and
stitched together; π− was generated in a single run.

- **LAr_P_cascade_0_9000.root / LAr_N_cascade_0_9000.root** : proton /
  neutron fracs, 5 nucleon-counting channels (0n0p, 0n1p, 1n0p, 1n1p,
  Other), 11 to 9002 MeV/c, 881 points per channel.
  Source: `/global/cfs/cdirs/dune/users/mazam/G4RW/G4Cascade/`
- **LAr_PiMinus_cascade_0_9000.root** : π− fracs, 5 channels (abs, cex,
  dcex, inel, prod), 11 to 9001 MeV/c.
  Source: `/global/cfs/cdirs/dune/users/mazam/G4RW/G4Cascade/`
- π+ is not included; it comes from the standard dune_pardata
  distribution: `$DUNE_PARDATA_DIR/Geant4ReweightInput/LAr_PiPlus_cascade_0_9000.root`
- **g4rw_cascade_proton_part1/2/3.fcl, g4rw_cascade_neutron_part1/2/3.fcl** :
  the 3-way momentum-chunked G4Cascade configs (NDivisions 293/293/292,
  boundaries offset by one momentum step to avoid duplicate/skipped
  points).
- **g4rw_cascade_piMinus.fcl** : single-shot config, full 11 to 9001
  MeV/c range.
  Each FCL includes `g4rw_material.fcl` and `g4rw_physlist.fcl`, which
  ship with the Geant4Reweight framework and are not duplicated here.
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/fcl/geant4reweight/`
- **stitch_results.py** : combines the 3 part files per species into the
  final fracs ROOT file. Patched from the upstream Geant4Reweight script
  (NuSoftHEP/Geant4Reweight) to add a PDG 2112 (neutron) branch; upstream
  only supports pion (211) and proton (2212).
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/geant4reweight_custom/build/install/bin/stitch_results.py`
- **proton_fracs_file_extended_range.md / neutron_fracs_file_extended_range.md** :
  notes on the chunking boundaries and stitch command per species,
  folded into `install_notes/README_G4RW_fracs_files.md`.
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/fcl/geant4reweight/`

### validation/fracs_files/

- **fracs_comparison.ipynb** : plots all four species' exclusive channel
  fractions vs. momentum, reproducing technote Figure 3. Requires
  `uproot`, `numpy`, `matplotlib`.
  Source: `/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/fracs_comparison.ipynb`

## geant4reweight_patch/

Source patches to Geant4Reweight (upstream: NuSoftHEP/Geant4Reweight),
combining nucleon exclusive-channel support with the Eq. 3.10
interaction-step weight fix. This is the state used for all production
weight generation; `geant4reweight_custom` was an earlier intermediate
build (nucleon patches only) and is not shipped here.

- **G4ProtonReweighter.cc / G4NeutronReweighter.cc** : give protons and
  neutrons the same 5-channel exclusive decomposition pions already have
  (technote Section IV A). Constructor channel list changed from
  `{"total"}` to `{"0n0p", "0n1p", "1n0p", "1n1p", "Other"}`;
  `GetInteractionSubtype()` counts outgoing protons/neutrons via
  `HasChild()`. π− needs no equivalent patch (`G4PiMinusReweighter.cc`
  already implements all five channels).
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/geant4reweight/`
- **G4ReweightParameterMaker.cc** : `case 2212` and `case 2112` extended
  from `{"total", "reac"}` to the 5 exclusive channels plus `"reac"`.
  K± remain `{"total", "reac"}`.
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/geant4reweight/`
- **G4Reweighter.cc / G4Reweighter.hh** : implements the full Eq. 3.10
  interaction-step weight (technote Section II B), replacing the bare
  scale factor `W_I = f_I` (Eq. 5) with the rate-ratio form
  `W_I = f_I * (sigma/sigma') * (1-e^(-dL*sigma'))/(1-e^(-dL*sigma))`
  (Eq. 3).
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/geant4reweight_eq310/`
- `G4Cascade.cc` is not included: it is identical to the upstream
  NuSoftHEP file. `G4PiPlusReweighter.cc` and `G4MultiReweighter.cc` are
  also not included (unmodified).

## weight_production/

Executables and configs that compute per-track G4RW weights from
patched EDepSim trajectories (technote Section V), producing the
`g4rwEvt` TTree (Table V). Links against the eq3.10-patched
Geant4Reweight build and the fracs files above.

- **do_weighting_v2_PiPlus.c++ / do_weighting_v2_PiMinus.c++ /
  do_weighting_v2_protons.c++ / do_weighting_v2_neutron.c++** : one
  executable per species (two species cannot share a G4ReweightManager
  in one process). Each computes weights at all 5 scale factors from
  Eq. 8 (`f = 0.34, 0.67, 1.33, 1.66, 2.00`), hardcoded in the source via
  repeated `SetAllParameterValues()` calls, not read from the FCL.
  Writes one `g4rwEvt` entry per particle: `file_idx`, `spill_id`,
  `event_id` (= RunId*1e6 + EventId, the CAF matching key), `run_id`,
  `track_id`, `pdg`, one `weight_*` branch per scale factor, `len`.
  Rock spills (`RunId >= 1e9`) are skipped.
- **reweight_PiPlus.fcl / reweight_PiMinus.fcl / reweight_proton.fcl /
  reweight_neutron.fcl** : per-species material, FracsFile path, and
  ParameterSet (the `"reac"` cut). The top-level `Weight`/`Sigma` fields
  are parsed but not applied; the 5 scale factors are hardcoded in the
  `.c++` files above.
- **material.fcl** : shared liquid argon definition, unmodified
  Geant4Reweight default.
- **Makefile** : 4 independent `g++` targets (`make do_weighting_v2_PiPlus`,
  etc, no combined `make all`), linking against `$GEANT4REWEIGHT_LIB`/
  `$GEANT4REWEIGHT_INC` (must point at the eq3.10 build), `cetlib`,
  `fhiclcpp`, `boost::program_options`, `edepsim_io`, ROOT.
- **run_weighting_PiPlus.sh / run_weighting_PiMinus.sh /
  run_weighting_proton.sh / run_weighting_neutron.sh** : run the species
  executable over all 1000 patched EDepSim files, skip missing
  inputs/existing outputs, hadd into
  `weights_1000_all_LAr_*_cascade_0_9000_5weights.root`.
  Source: `/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/multi_eq310/`

See `install_notes/README_G4RW_weight_production.md` for the full build
and environment setup.

---

## validation/

Cross-checks confirming the pipeline above behaves as documented.
Pions and protons are validated in parallel subfolders with the same
structure.

### validation/edepsim_patch/

Confirms the EDepSim patch takes effect for all six species;
reproduces technote Table I.

- **check_patched_traj_nkaon_all1000.C** : loops all 1000 files of
  `patched_withRock_genieNu_v2`, records trajectory-point-count
  statistics per species. Run: `root -l -b -q check_patched_traj_nkaon_all1000.C`
  (processes about 941k events).
- **check_patched_traj_nkaon.C** : single-file version. Run:
  `root -l -b -q 'check_patched_traj_nkaon.C("/path/to/file.root")'`
- **check_all1000_output.txt** : saved output of the all-1000-files run
  (the numbers behind Table I).
  Source: `/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW/batch_scripts/`

### validation/caf_matching/pions/

Confirms the three-way CAF/EDepSim/G4RW-weights lookup key (technote
Section VII A, Table XII) resolves consistently.

- **caf_edepsim_weights_matcher.C** : picks 15 pions (5 each from 3
  weight-value bands) and traces each through the CAF, 3 EDepSim
  samples (official unpatched, old patched, current patched), and 3
  weight output files, printing one table per track. Confirms the
  matching convention: file index direct match, interaction key
  `RunId*1e6 + EventId`, particle key `G4ID`/`TrackId`/`track_id`, no
  offset. For `RunId` above about 2147, `RunId*1e6` overflows a 32-bit
  int; the weight-writer and this matcher use the same arithmetic, so
  the wrapped value matches on both sides and the lookup still works.
  Requires the pion CAF output (`test/on_CAFS/31Aug2026/output.root` or any `caf_output.root` file)
  and the EDepSim/weight files it points to; not shipped here (size).
- **find_example_pions_from_caf_list.C** : scans `caf_pion_list.root`
  and picks candidates in three weight-value bands (A: 0.05 to 0.30,
  B: 0.45 to 0.55, C: 0.95 to 1.05 at f=1.33), sorted by distance to
  band center. Feeds the `validate_pion_example_*` scripts and
  `caf_edepsim_weights_matcher.C`'s target list.
- **validate_pion_example_1.c++ / _2.c++ / _2_exact.c++ / _2_eq3p10.c++ /
  _3.c++ / _4.c++ / _5.c++ / _6.c++** : single-track weight validations.
  `_1` is an independent example (RunId 0, EventId 3, TrackId 6).
  `_2`/`_2_exact`/`_2_eq3p10` all target the same track (RunId 0,
  EventId 7, TrackId 3), the technote's own Table VI example, at
  different formula eras (`_2_eq3p10` is current Eq. 3.10); kept
  together for comparison. `_3` to `_6` each target a distinct band
  example (`_3` = Band A #1, `_4` = Band B #3, `_5` = Band B #1
  anomaly, `_6` = Band C #1), matching `caf_edepsim_weights_matcher.C`'s
  target list.
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation/`
- **dump_caf_pion_list_v2.C** : dumps CAF-matched truth pions with
  `file_idx`/`event_id`/`track_id` stored undecoded, the convention
  `purity_dlp_multiplicity.cc` itself uses for its weight lookup.
  Produces `caf_pion_list_v2.root`, consumed by
  `plot_pion_weights_2way_eq3p10_v2.C`.
- **dump_caf_pion_list_eq3p10.C** : dumps CAF-matched truth pions with
  `RunId`/`EventId` decoded from `nu_id` (integer division/modulo), the
  schema `scan_all_proc_subproc.c++`, `plot_caf_weight_vs_length_2d.C`,
  and `replot_proc_subproc.C` expect. Produces
  `caf_pion_list_eq3p10.root`.
  Note: `scan_all_proc_subproc.c++` and `plots/replot_proc_subproc.C`
  hardcode the input filename `caf_pion_list.root`; rename
  `caf_pion_list_eq3p10.root` to that, or edit the filename in the
  script, before running them.
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation/codes_for_slides/`

#### validation/caf_matching/pions/plots/

- **plot_running_weight_vs_length.C** : running (cumulative) and
  per-step G4RW weight vs. track length for Examples 1 to 6, patched and
  unpatched, all 5 scale factors overlaid. Labels Examples 3 and 5 as
  "no inelastic tag found" and 4/6 as "clean," illustrating the
  two-population structure behind Figure 8. Reads the TGraphs saved by
  `validate_pion_example_*.c++`.
- **plot_caf_weight_vs_length_2d.C** : 2D weight vs. track length, one
  plot per scale factor, truth and reco length, reproducing Figure 7
  (left). Also splits truth length by containment (+/-55 cm cube) to
  test whether the low-weight, no-tag population corresponds to pions
  leaving the fiducial volume. Input: `caf_pion_list_eq3p10.root`.
- **plot_pion_weights_2way_eq3p10_v2.C** : reproduces Figure 2, truth π+
  G4RW weight distribution at f=1.33, old formula (`W_I = f_I`) vs. the
  full Eq. 3.10 rate-ratio formula, looked up fresh by
  `(file_idx, event_id, track_id)`. Input: `caf_pion_list_v2.root`
  (used only as an identity list, not for its embedded weights).
- **plot_truth_pion_weights_patched.C** : reproduces Figure 6 (left),
  truth π+ weight distribution from matched interactions, all 5 scale
  factors overlaid. Reads `w_patched[5]` directly from the CAF's
  `syst_g4rw_pion_truth` tree.
- **replot_proc_subproc.C** : two outputs. A plain re-draw of the
  proc/subproc composition histogram from `scan_all_proc_subproc.root`.
  And, via `plot_weight_vs_length_by_tag()`, Figure 8: weight vs. track
  length for truth π+ at f=1.33, split by the last point's Geant4
  process/subprocess (4/121 hadron inelastic, 2/22 EM/scintillation,
  6/201 decay). Inputs: `scan_all_proc_subproc.root` (from
  `validation/proc_subproc/`) and `caf_pion_list.root`.
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_pion_Eq3p10/`
  (`plot_pion_weights_2way_eq3p10_v2.C`, `plot_truth_pion_weights_patched.C`)
  and `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation/`
  (`replot_proc_subproc.C`, `plot_caf_weight_vs_length_2d.C`,
  `plot_running_weight_vs_length.C`)

### validation/proc_subproc/

Classifies the Geant4 process/subprocess code at each track's last
recorded trajectory point (technote Table VII, Figure 8).

- **scan_all_proc_subproc.c++** : for every pion in a CAF-matched
  candidate list, opens the patched EDepSim trajectory and records the
  last point's (proc, subproc) pair, the same quantity
  `do_weighting_v2_*`'s last-step check reads. Produces the Table VII
  percentages (for example 70.3%/26.4%/3.3% for 4/121, 2/22, 6/201 in
  π+) and a per-track tree for `replot_proc_subproc.C`. Input:
  `caf_pion_list.root` (see note under `dump_caf_pion_list_eq3p10.C`
  above).
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation/scan_all_proc_subproc.c++`

### validation/caf_matching/protons/

Proton analogs of the pion suite above.

- **caf_edepsim_weights_matcher_proton.C** : traces one proton
  (RunId 10, EventId 16, TrackId 2, the same target as
  `validate_proton_example_7.c++`) through CAF, EDepSim, and the weight
  files. Uses the pre-eq3.10 (`_0_4000`) weight files; the matching
  logic itself does not depend on which weight era is used.
- **dump_caf_proton_list.C / dump_caf_proton_list_eq3p10.C** : dump
  CAF-matched truth protons with `RunId`/`EventId` decoded from
  `nu_id`, same convention in both. Produce `caf_proton_list.root` and
  `caf_proton_list_eq3p10.root`.
- **find_example_protons_from_caf_list.C** : classifies each proton as
  UPPER (`w = f*exp((1-f)S)`, correctly tagged), LOWER
  (`w = exp((1-f)S)`, missing tag), or FLAT (never reweighted), by
  inverting the weight formula and checking self-consistency against
  the saved w_2.00.
- **validate_proton_example_1.c++ through _7.c++** : single-track
  validations, each targeting a distinct example track.
- **scan_all_proc_subproc_proton.c++** : proton analog of the Table
  VII-style composition scan.
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/`

#### validation/caf_matching/protons/plots/

- **plot_caf_weight_vs_length_2d_proton.C** : proton analog of
  `plot_caf_weight_vs_length_2d.C` (pdg to 2212, proton file/tree).
  Reproduces Figure 7 (right). Makes no assumption that the pion
  explanation for any two-trend structure transfers; protons cannot
  decay, so the last-point classification differs.
- **replot_proc_subproc.C** : proton-specific, not a copy of the pion
  version. The stored `h_proc_subproc` histogram had accumulated stale
  cycles from repeated scan runs, so this version derives the tag-count
  histogram fresh from `proc_subproc_by_track` every time. Legend omits
  "decay" (protons cannot decay).
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/`

## fitting/

Continuous SI weight parameterization vs. true track length (technote
Section VI): polynomial fits of orders 3 to 7, separated by last-point
tag (hadron inelastic 4/121 vs. other), fit range [5, 290] cm.

- **pions/plot_caf_weight_vs_length_2d_piPlus_rawL_equivalent_allOrders.C** :
  produces Table XI and Figure 11 (π+). pol3 is adopted as final (Eq. 9)
  on the technote's stability argument across orders.
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation/fitting/`
- **protons/plot_caf_weight_vs_length_2d_proton_rawL_equivalent_allOrders.C** :
  produces Table X and Figure 10 (protons), same method.
  Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/fitting/`

---

## analysis_plots/proton_dataMC/

Data vs. MC proton distributions with systematic uncertainty bands
(technote Section VIII A, Figures 12 to 13).

- **purity_dlp_multiplicity.cc** : CAF-level proton selection,
  truth-backtracking, and the CAF/G4RW weight matching described in
  Section VII A (Ref. 16, commit ab109f6: reconstructed fill lines
  3300 to 3324, truth fill lines 3328 to 3347). Produces `output.root`,
  consumed by every plotting script below (`reco_matched_proton_length`,
  `reco_matched_proton_cosL`, `syst_g4rw_proton`,
  `syst_g4rw_proton_reco`, `syst_g4rw_nonproton`).
- **data.cc** : produces the data-side histograms
  (`recoProtonLengthAll`, `recoProtonCosLAll`) read from `FILE_DATA`.
  Source: https://github.com/mbilalazam/proton_plots/blob/main/data.cc
- **plot_dataMC_recoProton_truthSplit_all.C** : Figure 12, SI panel.
  Builds the +/-1 sigma SI band from species-split G4RW weights
  (proton/π+/π− by truePDG, at f=1.33/f=0.67). Also produces
  lines-only, band+lines, and mean-weight-per-bin plots verifying that
  the band's fractional width equals the per-bin mean G4RW weight
  deviation from 1.
- **plot_dataMC_proton_GENIE.C** : Figure 12, GENIE panel. Converts the
  pre-computed GENIE sum-of-dials band (from
  `nusyst_step2_reco_output.root`, produced by the nd_nusyst framework,
  external to this repo) to a fractional uncertainty, re-centered on
  the data-normalized MC.
- **plot_dataMC_proton_detector_band_MR65sys.C** : Figure 12, detector
  panel. Per-bin RMS band from a MiniRun6.5 baseline plus 10
  detector/physics variations (Birks_A_b, Birks_k_b, Long_Diff,
  Temperature, W_ion, each up/down). Box_Model, e_lifetime, and
  Tran_Diff variations exist but are excluded from this band.
- **plot_dataMC_proton_combined_band.C** : Figure 12, total panel.
  GENIE + SI + detector in quadrature, rebuilding all three bands
  internally. Also prints a per-bin table of each band's absolute and
  fractional width.
- **plot_truthProton_band.C** : Figure 13, all four panels: reco and
  truth-level proton candidates (2257, correctly backtracked, +/-1
  sigma SI band) and reco/truth non-proton candidates (weighted by true
  π+/π−, other species at w=1), raw MC counts, no data overlay or
  normalization.

All scripts read data from:
`/global/u1/m/mazam/sandbox/v11/sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root`
and the eq3.10 weight files from `weight_production/`'s output:
`patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_{P,PiPlus,PiMinus}_cascade_0_9000_5weights.root`.
Source: `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/total_systematics/`
and `/global/cfs/cdirs/dune/users/mazam/softwares/patches/2x2_g4rw/validation_proton/codes_for_slides/proton_plots/`

## analysis_plots/multiplicity_dataMC/

Data vs. MC charged-track multiplicity with systematic uncertainty
bands (technote Section VIII B, Figure 14).

- **nominal_mc.cc** : CAF-level track-multiplicity counting
  (`minTrkLength = 3`), truth-backtracking, CAF/G4RW weight matching.
  Produces `nominal_mc.root` (`track_multGood` histogram,
  `syst_g4rw_mult_reco` tree with `nu_id`/`cafFileIdx`/`G4ID`/`truePDG`).
- **nominal_data.cc** : produces `nominal_data.root` (`track_mult`
  histogram) from the sandbox_v11 data sample.
  Build and run (see `instructions.txt`):
  ```
  g++ -std=c++17 -O2 -g -Wall nominal_mc.cc -o plotter ...
  ./plotter MiniRun6p5.txt nominal_mc.root 1
  ./plotter sandbox_v11.txt nominal_data.root 1
  ```
  `MiniRun6p5.txt` / `sandbox_v11.txt` list the CAF files to process.
- **g4rw_dataMC_trackMult.C** : Figure 14, SI panel (top). Each hadron
  track (proton, π+, π−) contributes its own weight independently, no
  event-level combination. 
- **nusyst_dataMC_trackMult.C** : Figure 14, GENIE panel (bottom left).
  Sum-of-dials band from the nd_nusyst framework's per-file systematic
  weights.
- **detector_dataMC_trackMult.C** : Figure 14, detector panel (bottom
  right). Per-bin RMS band from a single MiniRun6.5 nominal
  (`nominal_mc.root` serves as both the RMS reference and the display
  curve, since all variations come from the same `nominal_mc.cc` code).
  Uses all 23 detector/physics variations.
- **run_purity_dlp_systematics.sh** : compiles `nominal_mc.cc` once,
  runs it against every `MR6p5_*.txt` variation file-list under the
  detector-systematics directory, writing one
  `nominal_mc_<label>_output.root` per variation to a single output
  folder.
- **compare_proton_SI_weights.C** : cross-check, not a technote figure.
  Compares SI weight distributions for protons as selected by this
  multiplicity study (`syst_g4rw_mult_reco`, L>2cm, no containment cut)
  against the dedicated proton study's selection
  (`syst_g4rw_proton_reco`, L>3cm, containment on). Also compares π+ vs.
  π− weight distributions within the multiplicity sample.
  Source: `/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/`
  and `/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/multiplicity_systematics/`

---

## common/

- **protoDUNEStyle.C** : shared ROOT plotting style (white background,
  no colored frames, protoDUNE font/label conventions, kBird palette,
  `TGaxis::SetMaxDigits(4)`), loaded by most plotting scripts in this
  repo. Copy or symlink into a script's working directory before
  running it.

## install_notes/

Build and run instructions, one file per component.

- **README_G4RW_edepsim_patch.md** : clone and patch edep-sim, build
  inside the `mjkramer/sim2x2:ndlar011` container (required: the 2x2
  GDML geometry uses G4Torus shapes that crash under the system
  Geant4), symlink the MiniRun5 GENIE input, run single/batch jobs,
  hadd outputs, verify the patch.
- **README_G4RW_fracs_files.md** : reproduce the proton, neutron, and
  π− fracs files: environment setup, chunk boundaries, running
  G4Cascade per chunk, stitching.
- **README_G4RW_geant4reweight_patches.md** : build Geant4Reweight with
  the nucleon exclusive-channel patches and the Eq. 3.10 fix combined.
- **README_G4RW_weight_production.md** : build and run the 4-species
  weight-production executables against the eq3.10 Geant4Reweight
  build.

## docs/

- **G4RW_Technote.pdf** : "Hadronic Secondary-Interaction Systematics
  for the DUNE ND-LAr Demonstrator (2x2) with Geant4Reweight" (M. B.
  Azam, Z. Djurcic). 
- **technote_source/** : LaTeX source. `main.tex` is the top-level
  document.
