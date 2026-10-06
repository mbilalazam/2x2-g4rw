#!/bin/bash

# Run via
# shifter --image=fermilab/fnal-wn-sl7:latest --module=cvmfs -- /bin/bash run_purity_dlp_systematics.sh

set -o pipefail   # dropped -u (unbound var) and -e (added back later) — vendor UPS scripts aren't nounset-safe

echo "==> Inside container"
echo "Hostname: $(hostname)"

# -----------------------------
# Load DUNE stack
# -----------------------------
echo "==> Sourcing setup_dune.sh ..."
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
if [ $? -ne 0 ] || ! command -v setup &> /dev/null; then
    echo "FATAL: setup_dune.sh failed or 'setup' command not available after sourcing." >&2
    exit 1
fi
echo "==> setup_dune.sh OK, 'setup' command is available."

run_setup () {
    echo "==> setup $*"
    setup "$@"
    if [ $? -ne 0 ]; then
        echo "FATAL: 'setup $*' failed." >&2
        exit 1
    fi
}

run_setup cmake v3_22_2
run_setup gcc v9_3_0
run_setup pycurl
run_setup ifdhc
run_setup geant4 v4_11_0_p01c -q e20:debug
run_setup dk2nugenie   v01_10_01k -q debug:e20
run_setup genie_xsec   v3_04_00 -q AR2320i00000:e1000:k250
run_setup genie_phyopt v3_04_00 -q dkcharmtau
run_setup jobsub_client
run_setup eigen v3_3_5
run_setup duneanaobj v03_06_01b -q e20:prof
run_setup hdf5 v1_10_5a -q e20
run_setup fhiclcpp v4_15_03 -q debug:e20

echo "==> All UPS products set up successfully."

set -e   # now go strict for the rest of the script (our own commands, not vendor scripts)

G4_cmake_file=$(find "${GEANT4_FQ_DIR}/lib" -name 'Geant4Config.cmake')
export Geant4_DIR=$(dirname "$G4_cmake_file")
export PATH=$PATH:$GEANT4_FQ_DIR/bin
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$GENIE_INC/GENIE

echo "==> Geant4_DIR = $Geant4_DIR"
echo "==> DUNEANAOBJ_INC = ${DUNEANAOBJ_INC:-UNSET}"
echo "==> DUNEANAOBJ_LIB = ${DUNEANAOBJ_LIB:-UNSET}"

# -----------------------------
# Paths
# -----------------------------
# Per-variation MR6p5_*.txt file lists still live here, one per subfolder --
# INPUT side only, unchanged from the proton-study production.
DST_ROOT="/global/u1/m/mazam/MiniRuns/detector_systematics/anl_2026_samples"

# New source: nominal_mc.cc (track multiplicity code), not the
# proton-specific purity_dlp_multiplicity.cc used previously.
SRC_CC="/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/nominal_mc.cc"

# ALL outputs (compiled binary, logs, ROOT files) now go into ONE shared
# folder, instead of being scattered back into each variation's own
# subfolder as before.
OUT_DIR="/global/cfs/cdirs/dune/users/mazam/talks/workshop_Sep2026/multiplicity_systematics/detector_systematic_output_root"
mkdir -p "$OUT_DIR"

# -----------------------------
# Compile ONCE, directly into OUT_DIR
# -----------------------------
echo "==> Compiling plotter ..."
g++ -std=c++17 -O2 -g -Wall "$SRC_CC" -o "$OUT_DIR/plotter" \
    $(root-config --cflags --glibs) \
    -I"$DUNEANAOBJ_INC" \
    -L"$DUNEANAOBJ_LIB" \
    -lduneanaobj_StandardRecord \
    -lduneanaobj_StandardRecord_dict \
    -lduneanaobj_StandardRecordProxy

echo "==> Compile succeeded: $OUT_DIR/plotter"

MC_ONLY=1

# -----------------------------
# For each variation's MR6p5_*.txt (read from its own subfolder under
# DST_ROOT), run the SAME plotter binary, but write the ROOT output and
# log into the single central OUT_DIR, named nominal_mc_<label>_output.root
# -----------------------------
find "$DST_ROOT" -type f -name "MR6p5_*.txt" -print0 | while IFS= read -r -d '' txt; do
    label=$(basename "$txt" .txt)

    out="$OUT_DIR/nominal_mc_${label}_output.root"
    log="$OUT_DIR/nominal_mc_${label}.log"

    echo "==> Running for $label (input: $txt, log: $log)"
    ( cd "$OUT_DIR" && ./plotter "$txt" "$out" "$MC_ONLY" > "$log" 2>&1 )
    status=$?
    if [ $status -ne 0 ]; then
        echo "    FAILED (exit $status) — see $log"
    else
        echo "    OK -> $out"
    fi
done

echo "==> All done."