## Building and Running the G4RW Weight Production (eq3.10, 4 species)

> Based on `instructions.txt` from
> `/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/multi_eq310/`.
> This is the actual production workflow — 4 species, linked against
> the eq3.10-patched Geant4Reweight build. (An earlier, superseded
> 2-species/pre-eq3.10 setup exists in
> `README_2x2_g4rw_install.md`/`geant4reweight_custom` and is not used
> here.)

## Phase 1 — Environment

```bash
shifter --image=fermilab/fnal-wn-sl7:latest --module=cvmfs -- /bin/bash

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup geant4reweight v01_21_02  -q e26:prof:s131
setup dune_pardata   v01_84_00
setup boost          v1_82_0    -q e26:prof
setup edepsim        v3_2_0e    -q e26:prof
setup root           v6_28_12   -q e26:p3915:prof
setup geant4         v4_11_2_p02 -q e26:prof

export SOFTWAREDIR=/global/cfs/cdirs/dune/users/mazam/softwares
export G4RW_CFS_EQ310=$SOFTWAREDIR/geant4reweight_eq310
```

## Phase 2 — Point at the eq310 library (critical override)

```bash
export GEANT4REWEIGHT_LIB=$G4RW_CFS_EQ310/build/install/lib
export GEANT4REWEIGHT_INC=$G4RW_CFS_EQ310/build/install/include
export LD_LIBRARY_PATH=$GEANT4REWEIGHT_LIB:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=/cvmfs/larsoft.opensciencegrid.org/products/boost/v1_80_0/Linux64bit+3.10-2.17-e26-prof/lib:$LD_LIBRARY_PATH

# Needed so the FracsFile relative path in reweight_PiPlus.fcl resolves
# (dune_pardata's Geant4ReweightInput/)
export FHICL_FILE_PATH=/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/multi_eq310:$FHICL_FILE_PATH
```

## Phase 3 — Build

```bash
cd /global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/multi_eq310
make do_weighting_v2_PiPlus
make do_weighting_v2_protons
make do_weighting_v2_PiMinus
make do_weighting_v2_neutron
```

Confirm each is linked against the eq310 build, not the production/UPS
one:
```bash
ldd do_weighting_v2_PiPlus   | grep ReweightBaseLib
ldd do_weighting_v2_protons  | grep ReweightBaseLib
ldd do_weighting_v2_PiMinus  | grep ReweightBaseLib
ldd do_weighting_v2_neutron  | grep ReweightBaseLib
```

## Phase 4 — Run (1000 files, merged + hadd'd per species)

```bash
chmod +x run_weighting_PiPlus.sh   && ./run_weighting_PiPlus.sh
chmod +x run_weighting_proton.sh   && ./run_weighting_proton.sh
chmod +x run_weighting_PiMinus.sh  && ./run_weighting_PiMinus.sh
chmod +x run_weighting_neutron.sh  && ./run_weighting_neutron.sh
```

Each batch script writes its merged output to
`$CFS/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/
weights_1000_all_LAr_<species>_cascade_0_9000_5weights.root`.

## Verify

```bash
root -l -e 'TFile f("<merged_output>.root"); \
  TTree* t=(TTree*)f.Get("g4rwEvt"); t->Print(); \
  t->Scan("file_idx:event_id:track_id:pdg:weight_1p33","","",10);'
```

## Notes

- The 5 scale factors (Eq. 8: `f ∈ {0.34, 0.67, 1.33, 1.66, 2.00}`) are
  hardcoded in each `.c++` file via repeated
  `reweighter->SetAllParameterValues({f})` calls — the FCL's top-level
  `Weight`/`Sigma` fields are parsed but unused.
- Input EDepSim: the patched, withRock, genie.nu sample
  (`patched_withRock_genieNu`) — see `edepsim_patch/`.
- See `weight_production/` in this repo for the source/FCL/Makefile/
  batch-script files themselves.