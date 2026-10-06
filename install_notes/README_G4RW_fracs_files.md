## Generating G4RW Exclusive-Channel Fracs Files (π−, proton, neutron)

> This note is based on `proton_fracs_file_extended_range.md`,
> `neutron_fracs_file_extended_range.md`, and `instructions.txt`,
> originally at:
> `/global/cfs/cdirs/dune/users/mazam/softwares/patches/fcl/geant4reweight/`

## Phase 0 — Paths and environment

```bash
export SOFTWAREDIR=/global/cfs/cdirs/dune/users/mazam/softwares
export PATCHDIR=$SOFTWAREDIR/patches/geant4reweight
export G4RW_CFS=$SOFTWAREDIR/geant4reweight_custom

shifter --image=fermilab/fnal-wn-sl7:latest --module=cvmfs -- /bin/bash

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup cetmodules v3_24_01 -q e26:prof
setup geant4reweight v01_21_02 -q e26:prof:s131
setup dune_pardata v01_84_00
setup boost v1_80_0 -q e26:prof

# re-export (doesn't persist into container)
export SOFTWAREDIR=/global/cfs/cdirs/dune/users/mazam/softwares
export PATCHDIR=$SOFTWAREDIR/patches/geant4reweight
export G4RW_CFS=$SOFTWAREDIR/geant4reweight_custom
```

## Phase 1 — Why chunking is needed

`G4Cascade` memory scales with `(NDivisions+1) × NCascades`. A single
0–9000 MeV/c run (NDivisions ≈ 880, NCascades 10000) exceeds NERSC's
memory budget and gets killed. The working fix: split the 11–9001 MeV/c
range into 3 chunks, each within the memory budget already proven safe
by the earlier 0–4000 MeV/c file.

π− did **not** need chunking — `g4rw_cascade_piMinus.fcl` runs the full
range in a single pass successfully.

## Phase 2 — Chunk boundaries (proton and neutron)

Same boundaries for both species (`NCascades: 10000` throughout):

| Part | Momentum range (MeV) | NDivisions | Output file |
|---|---|---|---|
| 1 | 11.0 – 3004.261364 | 293 | `LAr_{P,N}_cascade_0_9000_part1.root` |
| 2 | 3014.477273 – 6007.738636 | 293 | `LAr_{P,N}_cascade_0_9000_part2.root` |
| 3 | 6017.954545 – 9001.0 | 292 | `LAr_{P,N}_cascade_0_9000_part3.root` |

Boundaries offset by one momentum step (~10.216 MeV) between consecutive
parts, so no momentum point is duplicated or skipped.

## Phase 3 — Run G4Cascade for each part

```bash
# Make sure the FCL is where G4Cascade looks for it
cp $SOFTWAREDIR/patches/fcl/geant4reweight/g4rw_cascade_proton_part1.fcl \
   $G4RW_CFS/build/fcl/

export LD_LIBRARY_PATH=$G4RW_CFS/build/install/lib:$LD_LIBRARY_PATH
export FHICL_FILE_PATH=$G4RW_CFS/build/fcl:$FHICL_FILE_PATH
export LD_LIBRARY_PATH=/cvmfs/larsoft.opensciencegrid.org/products/boost/v1_82_0/Linux64bit+3.10-2.17-e26-prof/lib:$LD_LIBRARY_PATH

cd $G4RW_CFS/build
./install/bin/G4Cascade -c fcl/g4rw_cascade_proton_part1.fcl
```

Repeat for `_part2.fcl`, `_part3.fcl`, and the neutron equivalents.
For π−, run `g4rw_cascade_piMinus.fcl` the same way — one pass, no
chunking.

## Phase 4 — Stitch the 3 parts together (proton, neutron only)

`stitch_results.py` (patched from the upstream Geant4Reweight version to
add a PDG 2112/neutron branch — see `fracs_files/stitch_results.py`)
concatenates each channel's TGraph across the input files in momentum
order:

```bash
python $G4RW_CFS/build/install/bin/stitch_results.py \
  /global/cfs/cdirs/dune/users/mazam/G4RW/G4Cascade/LAr_P_cascade_0_9000.root \
  2212 \
  /global/cfs/cdirs/dune/users/mazam/G4RW/G4Cascade/LAr_P_cascade_0_9000_part1.root \
  /global/cfs/cdirs/dune/users/mazam/G4RW/G4Cascade/LAr_P_cascade_0_9000_part2.root \
  /global/cfs/cdirs/dune/users/mazam/G4RW/G4Cascade/LAr_P_cascade_0_9000_part3.root
```

Argument order: output file, PDG code (2212 = proton, 2112 = neutron),
then input files in ascending momentum order.

## Result

- `LAr_P_cascade_0_9000.root`, `LAr_N_cascade_0_9000.root`: 881 points
  per channel, spanning 11–9002 MeV/c (the 9002 vs. 9001 endpoint is
  normal G4Cascade behavior, also seen in the π− file — not specific to
  chunking). Sum-to-1 check passes across the full range with no
  discontinuity at either chunk boundary.
- `LAr_PiMinus_cascade_0_9000.root`: produced directly, no stitching.

## Verify

```bash
ls -lh /global/cfs/cdirs/dune/users/mazam/G4RW/G4Cascade/LAr_P_cascade_0_9000.root
```

See also `validation/fracs_files/fracs_comparison.ipynb`, which plots
all four species' fracs files (reproduces technote Figure 3).