#!/usr/bin/env bash
# run_g4rw_edepsim_batched.sh
#
# Runs patched edep-sim over GENIE files in groups of 10.
# Each group produces 10 individual EDEPSIM ROOT files, hadds them
# into one file, then deletes the 10 intermediates.
#
# Output: /global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/
# Naming: MiniRun5_1E19_RHC.nu.XXXXXXX.G4RW.EDEPSIM.root
#
# Usage:
#   bash run_g4rw_edepsim_batched.sh <start_group> <end_group> [genie_offset]
#
# Examples:
#   bash run_g4rw_edepsim_batched.sh 0 99          # groups 0-99,    GENIE files 0-999     (1st batch)
#   bash run_g4rw_edepsim_batched.sh 100 199 1000  # groups 100-199, GENIE files 1000-1999 (2nd batch)
#   bash run_g4rw_edepsim_batched.sh 200 299 2000  # groups 200-299, GENIE files 2000-2999 (3rd batch)
#   bash run_g4rw_edepsim_batched.sh 300 399 3000  # groups 300-399, GENIE files 3000-3999 (4th batch)
#   bash run_g4rw_edepsim_batched.sh 400 499 4000  # groups 400-499, GENIE files 4000-4999 (5th batch)
#   bash run_g4rw_edepsim_batched.sh 500 599 5000  # groups 500-599, GENIE files 5000-5999 (6th batch)
#   bash run_g4rw_edepsim_batched.sh 600 699 6000  # groups 600-699, GENIE files 6000-6999 (7th batch)
#   bash run_g4rw_edepsim_batched.sh 700 799 7000  # groups 700-799, GENIE files 7000-7999 (8th batch)
#   bash run_g4rw_edepsim_batched.sh 800 899 8000  # groups 800-899, GENIE files 8000-8999 (9th batch)
#   bash run_g4rw_edepsim_batched.sh 900 999 9000  # groups 900-999, GENIE files 9000-9999 (10th batch)

set -e

# ── User settings ──────────────────────────────────────────────────────────
EDEPSIM_DIR=/global/u1/m/mazam/2x2_sim/run-edep-sim
# OUTDIR=/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu
# TMPDIR=/global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test
LOGDIR=${TMPDIR}/logs_batched

# Extended patch (pi/p/n/K) output — was: EdepSim_G4RW_nu
OUTDIR=/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/patched_withRock_genieNu_v2
# Separate staging area for the extended-patch run — was: bilal_test
TMPDIR=/global/u1/m/mazam/2x2_sim/run-edep-sim/bilal_test_v2

# ──────────────────────────────────────────────────────────────────────────

START_GROUP=${1:-0}   # first group index  (default 0)
END_GROUP=${2:-1}     # last group index   (default 1, i.e. 2 groups = 20 files)
GENIE_OFFSET=${3:-0}  # offset added to GENIE file indices (e.g. 1000 for second batch)

mkdir -p "$OUTDIR" "$LOGDIR"

# Set Shifter environment (same as before)
export ARCUBE_RUNTIME=SHIFTER
export ARCUBE_CONTAINER=mjkramer/sim2x2:ndlar011
export ARCUBE_EDEP_MAC=macros/dune-nd.mac
# export ARCUBE_GEOM_EDEP=geometry/Merged2x2MINERvA_v4/Merged2x2MINERvA_v4_noRock.gdml
# withRock confirmed as the physically correct geometry — was: noRock
export ARCUBE_GEOM_EDEP=geometry/Merged2x2MINERvA_v4/Merged2x2MINERvA_v4_withRock.gdml

export ARCUBE_GENIE_NAME=MiniRun5_1E19_RHC.genie.nu
export ARCUBE_OUTDIR_BASE=${TMPDIR}
export ARCUBE_LOGDIR_BASE=${LOGDIR}
export ARCUBE_OUT_NAME=G4RW_test.edep.nu
export ARCUBE_RUN_OFFSET=0
export SHIFTER_VOLUME="/global/homes/m/mazam/edepsim-container/bin:/opt/generators/edep-sim/install/bin;/global/homes/m/mazam/edepsim-container/lib:/opt/generators/edep-sim/install/lib"

cd "$EDEPSIM_DIR"

for group in $(seq "$START_GROUP" "$END_GROUP"); do

    FIRST_IDX=$(( (group - START_GROUP) * 10 + GENIE_OFFSET ))
    LAST_IDX=$(( FIRST_IDX + 9 ))

    # Output hadd file name — zero-padded group index
    GROUP_STR=$(printf "%07d" "$group")
    HADD_OUT=${OUTDIR}/MiniRun5_1E19_RHC.nu.${GROUP_STR}.G4RW.EDEPSIM.root

    echo "========================================"
    echo "Group ${group}: GENIE files ${FIRST_IDX}-${LAST_IDX}"
    echo "Output: ${HADD_OUT}"
    echo "========================================"

    # Skip if output already exists
    if [[ -f "$HADD_OUT" ]]; then
        echo "Output already exists, skipping group ${group}."
        continue
    fi

    # ── Step 1: Run 10 edep-sim jobs in parallel ──────────────────────────
    echo "Running edep-sim for indices ${FIRST_IDX}-${LAST_IDX}..."
    for i in $(seq "$FIRST_IDX" "$LAST_IDX"); do
        ARCUBE_INDEX=$i ./run_edep_sim.sh \
            > "${LOGDIR}/edepsim_${i}.log" 2>&1 &
    done
    wait
    echo "All 10 edep-sim jobs done for group ${group}."

    # ── Step 2: Collect the 10 output files ──────────────────────────────
    INPUT_FILES=""
    for i in $(seq "$FIRST_IDX" "$LAST_IDX"); do
        IDX_STR=$(printf "%07d" "$i")
        SUB_DIR=$(printf "%07d" $(( i / 1000 * 1000 )))
        F="${TMPDIR}/run-edep-sim/G4RW_test.edep.nu/EDEPSIM/${SUB_DIR}/G4RW_test.edep.nu.${IDX_STR}.EDEPSIM.root"
        if [[ ! -f "$F" ]]; then
            echo "ERROR: Expected output file not found: $F"
            exit 1
        fi
        INPUT_FILES="${INPUT_FILES} ${F}"
    done

    # ── Step 3: hadd into one file ────────────────────────────────────────
    echo "Hadding into ${HADD_OUT}..."
    shifter --image=mjkramer/sim2x2:ndlar011 --module=cvmfs \
        /bin/bash -c "
            cd /opt/generators/root && source bin/thisroot.sh && cd -
            hadd ${HADD_OUT} ${INPUT_FILES}
        "
    echo "hadd done."

    # ── Step 4: Delete the 10 intermediate files ─────────────────────────
    echo "Cleaning up intermediate files..."
    for i in $(seq "$FIRST_IDX" "$LAST_IDX"); do
        IDX_STR=$(printf "%07d" "$i")
        SUB_DIR=$(printf "%07d" $(( i / 1000 * 1000 )))
        rm -f "${TMPDIR}/run-edep-sim/G4RW_test.edep.nu/EDEPSIM/${SUB_DIR}/G4RW_test.edep.nu.${IDX_STR}.EDEPSIM.root"
        rm -f "${TMPDIR}/tmp/run-edep-sim/G4RW_test.edep.nu/G4RW_test.edep.nu.${IDX_STR}.EDEPSIM.root"
    done
    echo "Cleanup done for group ${group}."

    echo "Group ${group} complete → ${HADD_OUT}"
    ls -lh "$HADD_OUT"
    echo ""

done

echo "All groups done."
echo "Output files:"
ls -lh "$OUTDIR"/