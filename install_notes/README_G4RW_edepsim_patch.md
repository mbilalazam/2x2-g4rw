## Phase 0 — One-time setup: Clone 2x2 Sim and patch edep-sim
```bash
git clone https://github.com/DUNE/2x2_sim.git
```
## Phase 1 — One-time setup: Clone and patch edep-sim

> **Do this once on a login node. Skip if already done.**

```bash
# NOTE: These UPS products are only needed to build outside the container.
# Do NOT load them when running inside the container (Phase 2).
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup cmake v3_22_2
setup gcc v9_3_0
setup geant4 v4_10_6_p01f -q e20:prof
setup root v6_26_06b -q e20:p3913:prof

cd /global/homes/m/mazam
git clone https://github.com/lbl-neutrino/edep-sim.git
cd edep-sim
git checkout 3.2.0.2x2.1
```

Apply the patch to `src/EDepSimPersistencyManager.cc` as described above or copy the patched file from the following directory:

```bash
/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW/patch/EDepSimPersistencyManager.cc
```

---

## Phase 2 — One-time setup: Build patched edep-sim inside the container

> **Do this once from a login or compute node. Skip if `/global/homes/m/mazam/edepsim-container/` already exists.**
>
> Building inside the container ensures the binary links against the correct
> Geant4 and ROOT versions. The container is required because the DUNE 2x2
> GDML geometry contains G4Torus shapes that crash with the system Geant4.

```bash
# Enter the container
shifter --image=mjkramer/sim2x2:ndlar011 --module=cvmfs /bin/bash

# Source the container's Geant4 and ROOT
source /opt/generators/geant4/bin/geant4.sh
source /opt/generators/root/bin/thisroot.sh

# Build
cd /global/homes/m/mazam/edep-sim
mkdir -p build_container && cd build_container

cmake .. \
    -DCMAKE_INSTALL_PREFIX=/global/homes/m/mazam/edepsim-container \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo

make -j8
make install

# Exit the container
exit
```

Verify the install:
```bash
ls /global/homes/m/mazam/edepsim-container/bin/   # should show edep-sim
ls /global/homes/m/mazam/edepsim-container/lib/   # should show libedepsim.so
```

---

## Phase 3 — One-time setup: Create symlink for GENIE input

> **Do this once. Skip if symlink already exists.**

```bash
mkdir -p /global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test/run-genie

ln -s /global/cfs/cdirs/dune/www/data/2x2/simulation/productions/MiniRun5_1E19_RHC/MiniRun5_1E19_RHC.genie.nu \
      /global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test/run-genie/MiniRun5_1E19_RHC.genie.nu
```

---

## Phase 4 — Run patched edep-sim

> **Run this from a compute node or login node.**
> The `SHIFTER_VOLUME` mount replaces the container's edep-sim binary and
> library with our patched versions at runtime.

```bash
cd /global/u1/m/mazam/2x2_sim/run-edep-sim

export ARCUBE_RUNTIME=SHIFTER
export ARCUBE_CONTAINER=mjkramer/sim2x2:ndlar011
export ARCUBE_EDEP_MAC=macros/dune-nd.mac
export ARCUBE_GEOM_EDEP=geometry/Merged2x2MINERvA_v4/Merged2x2MINERvA_v4_withRock.gdml
export ARCUBE_GENIE_NAME=MiniRun5_1E19_RHC.genie.nu
export ARCUBE_OUTDIR_BASE=/global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test
export ARCUBE_LOGDIR_BASE=/global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test/logs
export ARCUBE_OUT_NAME=G4RW_test.edep.nu
export ARCUBE_RUN_OFFSET=0
export SHIFTER_VOLUME="/global/homes/m/mazam/edepsim-container/bin:/opt/generators/edep-sim/install/bin;/global/homes/m/mazam/edepsim-container/lib:/opt/generators/edep-sim/install/lib"
```

**Run a single file (index 0):**
```bash
ARCUBE_INDEX=0 ./run_edep_sim.sh 2>&1 | tee bilal_test/run_edep_sim_output_0.txt
```

**Run 10 files in parallel:**
```bash
for i in $(seq 0 9); do
    ARCUBE_INDEX=$i ./run_edep_sim.sh 2>&1 | \
    tee /global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test/run_edep_sim_output_${i}.txt &
done
wait
echo "All 10 done"
```

Output files will appear at:
```bash
bilal_test/run-edep-sim/G4RW_test.edep.nu/EDEPSIM/0000000/G4RW_test.edep.nu.XXXXXXX.EDEPSIM.root
```

---

## Phase 5 — Hadd output files

```bash
shifter --image=mjkramer/sim2x2:ndlar011 --module=cvmfs /bin/bash
source /opt/generators/root/bin/thisroot.sh

hadd /global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test/G4RW_test.edep.nu.hadd.root \
    /global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test/run-edep-sim/G4RW_test.edep.nu/EDEPSIM/0000000/G4RW_test.edep.nu.000000*.EDEPSIM.root

exit
```

---

## Phase 6 — Verify the patch worked

```bash
shifter --image=mjkramer/sim2x2:ndlar011 --module=cvmfs /bin/bash
source /opt/generators/root/bin/thisroot.sh
export LD_LIBRARY_PATH=/opt/generators/edep-sim/install/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=/opt/generators/edep-sim/install/include/EDepSim:$ROOT_INCLUDE_PATH

root -l -b -q check_patched_traj.C
exit
```

**Expected output if patch is working (averaged over 10 EdepSim files correspodning to 100 GENIE NU files and 10 CAF file, 9264 events):**

| Particle | Count | Avg pts | Min | Max |
|---|---|---|---|---|
| Pions | 16,722 | ~45 | 2 | 853 |
| Protons | 268,132 | ~9 | 2 | 3,638 |

Unpatched MiniRun5 baseline: pion avg ~2.5, proton avg ~2.1.

---

## File size comparison

| Sample | Size | Notes |
|---|---|---|
| Patched (hadded) | ~460 MB | Full π/p trajectory steps |
| MiniRun5 unpatched (hadded) | ~489 MB | Geometric pruning active |

The patch does **not** inflate file sizes.

---

## Cleanup (if re-running)

```bash
rm -f /global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test/run-edep-sim/G4RW_test.edep.nu/EDEPSIM/0000000/*.root
rm -f /global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test/tmp/run-edep-sim/G4RW_test.edep.nu/*.root
rm -f /global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test/run_edep_sim_output_*.txt
```

---

## Phase 7 — Submit as a Batch Job
Copy the following script and run it using `bash run_g4rw_edepsim_batched.sh <start_group> <end_group>` from the compute node.

```bash
/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW/batch_scripts/run_g4rw_edepsim_batched.sh
```
A group means a set of 10 files. If one wants to process 20 GENIE files (from 0-19) to produce 2 hadded Edepsim files, one would use `bash run_g4rw_edepsim_batched.sh 0 1`.

Additionlly, following script can be used to compare patched and unpatched Edepsim files:

```bash
/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW/batch_scripts/check_patched_traj.C
```

## Key file locations

| Item | Path |
|---|---|
| Patched edep-sim source | `/global/homes/m/mazam/edep-sim/src/EDepSimPersistencyManager.cc` |
| Container build install | `/global/homes/m/mazam/edepsim-container/` |
| UPS build install | `/global/homes/m/mazam/edepsim-local/` |
| 2x2_sim repo | `/global/u1/m/mazam/2x2_sim/` |
| Test output directory | `/global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test/` |
| Diagnostic script | `check_patched_traj.C` (in run-edep-sim dir) |
| MiniRun5 GENIE GTRAC | `/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/MiniRun5_1E19_RHC/MiniRun5_1E19_RHC.genie.nu/GTRAC/` |
| MiniRun5 EDEPSIM spills | `/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/MiniRun5_1E19_RHC/MiniRun5_1E19_RHC.spill/EDEPSIM_SPILLS/` |
| Production container | `mjkramer/sim2x2:ndlar011` |
| edep-sim inside container | `/opt/generators/edep-sim/install/bin/edep-sim` |

---
