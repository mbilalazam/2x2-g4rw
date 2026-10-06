## Building Geant4Reweight with Nucleon Exclusive Channels + Eq. 3.10 Fix

> Synthesized from `README_geant4reweight_custom_install.md` and
> `README_geant4reweight_eq3p10_install.md`.
>
> Two builds exist on NERSC: `geant4reweight_custom` (an earlier,
> intermediate build — nucleon patches only, no eq3.10 fix) and
> `geant4reweight_eq310` (the actual production build — nucleon patches
> + eq3.10 fix combined). This note documents the production path only.

## Phase 0 — Paths and environment

```bash
ssh mazam@perlmutter.nersc.gov

export SOFTWAREDIR=/global/cfs/cdirs/dune/users/mazam/softwares
export PATCHDIR=$SOFTWAREDIR/patches/geant4reweight
export G4RW_CFS_EQ310=$SOFTWAREDIR/geant4reweight_eq310

mkdir -p $G4RW_CFS_EQ310
cd $G4RW_CFS_EQ310

shifter --image=fermilab/fnal-wn-sl7:latest --module=cvmfs -- /bin/bash

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup cetmodules v3_24_01 -q e26:prof
setup geant4reweight v01_21_02 -q e26:prof:s131
setup dune_pardata v01_84_00
setup boost v1_80_0 -q e26:prof

# re-export (doesn't persist into container)
export SOFTWAREDIR=/global/cfs/cdirs/dune/users/mazam/softwares
export PATCHDIR=$SOFTWAREDIR/patches/geant4reweight
export G4RW_CFS_EQ310=$SOFTWAREDIR/geant4reweight_eq310
```

## Phase 1 — Clone and pin to the production commit

```bash
cd $G4RW_CFS_EQ310
git clone https://github.com/NuSoftHEP/Geant4Reweight.git src
cd src
git checkout fdde49eec00d1dcbb8409d2cdd1825a15bb88556   # tag v01_21_02
git rev-parse HEAD   # confirm it reads fdde49ee...
```

## Phase 2 — Apply the nucleon exclusive-channel patches

```bash
cp $PATCHDIR/G4ProtonReweighter.cc \
   $G4RW_CFS_EQ310/src/geant4reweight/ReweightBase/G4ProtonReweighter.cc

cp $PATCHDIR/G4NeutronReweighter.cc \
   $G4RW_CFS_EQ310/src/geant4reweight/ReweightBase/G4NeutronReweighter.cc

cp $PATCHDIR/G4ReweightParameterMaker.cc \
   $G4RW_CFS_EQ310/src/geant4reweight/PropBase/G4ReweightParameterMaker.cc

# Verify
grep -n "0n0p" $G4RW_CFS_EQ310/src/geant4reweight/ReweightBase/G4ProtonReweighter.cc
grep -n "0n0p" $G4RW_CFS_EQ310/src/geant4reweight/ReweightBase/G4NeutronReweighter.cc
grep -n "0n0p" $G4RW_CFS_EQ310/src/geant4reweight/PropBase/G4ReweightParameterMaker.cc
```

G4Cascade.cc is **not** patched — the upstream file already has native
nucleon-counting classification; no change needed there.

## Phase 3 — Apply the Eq. 3.10 interaction-weight fix

```bash
cp $SOFTWAREDIR/patches/geant4reweight_eq310/G4Reweighter.cc \
   $G4RW_CFS_EQ310/src/geant4reweight/ReweightBase/G4Reweighter.cc

cp $SOFTWAREDIR/patches/geant4reweight_eq310/G4Reweighter.hh \
   $G4RW_CFS_EQ310/src/geant4reweight/ReweightBase/G4Reweighter.hh

# Verify
grep -n "GetInteractionRateFactor" $G4RW_CFS_EQ310/src/geant4reweight/ReweightBase/G4Reweighter.cc
grep -n "IsOtherInteraction" $G4RW_CFS_EQ310/src/geant4reweight/ReweightBase/G4Reweighter.hh
```

## Phase 4 — Build

> **Known NERSC issue:** building directly on `$G4RW_CFS_EQ310` (CFS)
> fails. `setup_for_development` calls cetmodules, which uses `flock()`
> to lock `CMakePresets.json` before writing it. CFS is served to
> Perlmutter via DVS, which does not support POSIX/BSD file locking, so
> this fails with `errno 524` ("operation not supported") — not
> transient, not fixable by retrying. Perlmutter `$SCRATCH` is native
> Lustre and supports `flock()` normally, so the build must happen
> there; only the final installed library/headers are placed back onto
> CFS via `CMAKE_INSTALL_PREFIX`.

```bash
export G4RW_SCRATCH_EQ310=$SCRATCH/geant4reweight_eq310
mkdir -p $G4RW_SCRATCH_EQ310
cd $G4RW_SCRATCH_EQ310

# Clone fresh on scratch (faster than copying, avoids carrying .git
# history / build artifacts across filesystems)
git clone https://github.com/NuSoftHEP/Geant4Reweight.git src
cd src
git checkout fdde49eec00d1dcbb8409d2cdd1825a15bb88556
git rev-parse HEAD

# Re-apply all five patches on the scratch clone (Phase 2 + Phase 3,
# same files, same destinations, just under $G4RW_SCRATCH_EQ310 now)
cp $PATCHDIR/G4ProtonReweighter.cc \
   $G4RW_SCRATCH_EQ310/src/geant4reweight/ReweightBase/G4ProtonReweighter.cc
cp $PATCHDIR/G4NeutronReweighter.cc \
   $G4RW_SCRATCH_EQ310/src/geant4reweight/ReweightBase/G4NeutronReweighter.cc
cp $PATCHDIR/G4ReweightParameterMaker.cc \
   $G4RW_SCRATCH_EQ310/src/geant4reweight/PropBase/G4ReweightParameterMaker.cc
cp $SOFTWAREDIR/patches/geant4reweight_eq310/G4Reweighter.cc \
   $G4RW_SCRATCH_EQ310/src/geant4reweight/ReweightBase/G4Reweighter.cc
cp $SOFTWAREDIR/patches/geant4reweight_eq310/G4Reweighter.hh \
   $G4RW_SCRATCH_EQ310/src/geant4reweight/ReweightBase/G4Reweighter.hh

cd $G4RW_SCRATCH_EQ310
mkdir -p build && cd build
source $G4RW_SCRATCH_EQ310/src/ups/setup_for_development -p

# Build dir on scratch; install prefix points back to CFS (the
# permanent, quota-tracked location)
cmake -DCMAKE_INSTALL_PREFIX=$G4RW_CFS_EQ310/build/install \
      -DCMAKE_BUILD_TYPE=$CETPKG_TYPE \
      $G4RW_SCRATCH_EQ310/src

make -j4 2>&1 | tee build.log
make install
```

Verify install (lands on CFS per `CMAKE_INSTALL_PREFIX` above):
```bash
ls $G4RW_CFS_EQ310/build/install/bin/
ls $G4RW_CFS_EQ310/build/install/lib/
```

> `$SCRATCH` is subject to NERSC's purge policy (~8 weeks unaccessed).
> The scratch src/build tree is disposable and fully reproducible from
> the steps above; only `$G4RW_CFS_EQ310/build/install` (on CFS) is the
> permanent artifact.

## Runtime — load the patched library

```bash
shifter --image=fermilab/fnal-wn-sl7:latest --module=cvmfs -- /bin/bash

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup geant4reweight v01_21_02 -q e26:prof:s131
setup dune_pardata v01_84_00
setup boost v1_82_0 -q e26:prof
setup edepsim v3_2_0e -q e26:prof
setup root v6_28_12 -q e26:p3915:prof
setup geant4 v4_11_2_p02 -q e26:prof

export SOFTWAREDIR=/global/cfs/cdirs/dune/users/mazam/softwares
export G4RW_CFS_EQ310=$SOFTWAREDIR/geant4reweight_eq310

export GEANT4REWEIGHT_LIB=$G4RW_CFS_EQ310/build/install/lib
export GEANT4REWEIGHT_INC=$G4RW_CFS_EQ310/build/install/include
export LD_LIBRARY_PATH=$GEANT4REWEIGHT_LIB:$LD_LIBRARY_PATH
```

## For collaborators (no reinstallation needed)

The patched libraries are pre-built on CFS at
`$G4RW_CFS_EQ310/build/install/`. Just run the Runtime block above to
load the environment.

## Patch summary

Five files modified from stock Geant4Reweight (NuSoftHEP/Geant4Reweight,
commit `fdde49ee...`, tag `v01_21_02`):

| File | Change |
|---|---|
| `G4ProtonReweighter.cc` | Constructor: `{"total"}` → 5 nucleon channels; `GetInteractionSubtype()`: nucleon-counting logic |
| `G4NeutronReweighter.cc` | Same pattern, written fresh from the proton version |
| `G4ReweightParameterMaker.cc` | `case 2212`/`case 2112`: `{"total","reac"}` → 5 channels + `"reac"` |
| `G4Reweighter.cc` | Full Eq. 3.10 rate-ratio weight (`GetInteractionRateFactor`), replacing the bare scale-factor approximation |
| `G4Reweighter.hh` | Declarations/helpers for the Eq. 3.10 fix |

`G4Cascade.cc` is **not** patched — upstream already supports nucleon
counting natively.

See also `geant4reweight_patch/` in this repo for the five patch files
themselves, and `fracs_files/install_notes/README_G4RW_fracs_files.md`
for how G4Cascade generates the fracs files these patches consume.