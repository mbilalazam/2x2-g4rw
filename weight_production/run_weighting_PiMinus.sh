#!/bin/bash
# run_weighting.sh
# Runs do_weighting_v2 over the first N EDepSim files and merges output.
#
# Usage:
#   chmod +x run_weighting.sh
#   ./run_weighting.sh
#
# Output: weights_all.root  (merged weights tree with file_idx branch)

# ── Configure these paths ─────────────────────────────────────────────────
# EDEPSIM_BASE=/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/MiniRun5_1E19_RHC/MiniRun5_1E19_RHC.spill/EDEPSIM_SPILLS/0000000 # Unpatched (original) Edepsim files

# EDEPSIM_BASE=/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/unpatched_withRock_genieNu # Unatched Edepsim files (SELF) 

EDEPSIM_BASE=/global/cfs/cdirs/dune/users/mazam/EdepSim_G4RW_nu/patched_withRock_genieNu # Patched Edepsim files with intermediate trajectories 

BINARY=./do_weighting_v2_PiMinus
FCL=reweight_PiMinus.fcl
N_FILES=1000          # number of files to process (0000000 .. 0000019)

# OUTPUT_MERGED=/global/cfs/cdirs/dune/users/mazam/G4RW/output_files_unpatched_official/weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root # Unpatched (original) Edepsim files

# OUTPUT_MERGED=/global/cfs/cdirs/dune/users/mazam/G4RW/output_files_unpatched_self/weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root # Unpatched (SELF) Edepsim files

# OUTPUT_MERGED=/global/cfs/cdirs/dune/users/mazam/G4RW/output_files_patched/weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root # Patched Edepsim files with intermediate trajectories

# withRock paths
# OUTPUT_MERGED=/global/cfs/cdirs/dune/users/mazam/G4RW/output_files_unpatched_official/weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root # Unpatched (original) Edepsim files

# OUTPUT_MERGED=/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/unpatched_withRock_genieNu/weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights_addChild_range_0_9000.root # Unpatched (SELF) Edepsim files

# OUTPUT_MERGED=/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu/weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root # Patched Edepsim files with intermediate trajectories

OUTPUT_MERGED=/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu_eq3p10/weights_1000_all_LAr_PiMinus_cascade_0_9000_5weights.root # Patched Edepsim files with intermediate trajectories
# ─────────────────────────────────────────────────────────────────────────

echo "=========================================="
echo "  G4RW batch weighting: $N_FILES files"
echo "=========================================="

OUTFILES=""
N_OK=0
N_FAIL=0

for IDX in $(seq 0 $((N_FILES - 1))); do

    # Zero-pad to 7 digits: 0 -> 0000000, 3 -> 0000003, etc.
    FILEIDX=$(printf "%07d" $IDX)

    # EDEPSIMFILE=$EDEPSIM_BASE/MiniRun5_1E19_RHC.spill.${FILEIDX}.EDEPSIM_SPILLS.root # Unpatched (original) Edepsim files
    EDEPSIMFILE=$EDEPSIM_BASE/MiniRun5_1E19_RHC.nu.${FILEIDX}.G4RW.EDEPSIM.root # Patched Edepsim files with intermediate trajectories 
    OUTFILE=weighted_test_v2.${FILEIDX}.root

    # Check input exists
    if [ ! -f "$EDEPSIMFILE" ]; then
        echo "[$FILEIDX] SKIP -- file not found: $EDEPSIMFILE"
        N_FAIL=$((N_FAIL + 1))
        continue
    fi

    # Skip if output already exists
    if [ -f "$OUTFILE" ]; then
        echo "[$FILEIDX] SKIP -- output already exists: $OUTFILE"
        OUTFILES="$OUTFILES $OUTFILE"
        N_OK=$((N_OK + 1))
        continue
    fi

    echo "[$FILEIDX] Processing $EDEPSIMFILE ..."
    $BINARY --c $FCL --i $EDEPSIMFILE --o $OUTFILE

    if [ $? -eq 0 ]; then
        echo "[$FILEIDX] OK --> $OUTFILE"
        OUTFILES="$OUTFILES $OUTFILE"
        N_OK=$((N_OK + 1))
    else
        echo "[$FILEIDX] FAILED"
        N_FAIL=$((N_FAIL + 1))
    fi

done

echo ""
echo "=========================================="
echo "  Processed: $N_OK OK,  $N_FAIL failed"
echo "=========================================="

# Merge all individual outputs into one file
if [ $N_OK -gt 0 ]; then
    echo "Merging into $OUTPUT_MERGED ..."
    hadd -f $OUTPUT_MERGED $OUTFILES
    echo "Done. $OUTPUT_MERGED ready."
    # Clean up individual files
    echo "Cleaning up individual files..."
    rm -f $OUTFILES
    echo "Removed individual weight files."
    echo ""
    echo "Verify with:"
    echo "  root -l -e 'TFile f(\"$OUTPUT_MERGED\"); TTree* t=(TTree*)f.Get(\"g4rwEvt\"); t->Print(); t->Scan(\"file_idx:event_id:track_id:pdg:weight\",\"\",\"\",10);'"
else
    echo "No output files to merge."
fi
