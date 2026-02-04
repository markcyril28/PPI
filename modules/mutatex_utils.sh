#!/usr/bin/env bash
#######################################################################
# MutateX Utility Functions
#
# Common functions for MutateX batch processing.
# Source this file in your scripts:
#   source "$SCRIPT_DIR/modules/mutatex_utils.sh"
#######################################################################

# Background cleanup function - runs during mutatex to keep disk usage low
# Usage: cleanup_worker <mutations_dir> [interval_seconds]
# Returns: PID of background process (use to kill later)
cleanup_worker() {
    local mutations_dir="$1"
    local interval="${2:-60}"  # Check every 60 seconds
    
    while true; do
        sleep "$interval"
        if [[ -d "$mutations_dir" ]]; then
            # Remove PDB files from completed mutation directories (those with Average_*.fxout)
            find "$mutations_dir" -type d -mindepth 2 -maxdepth 2 | while read -r mutdir; do
                if ls "$mutdir"/Average_*.fxout 1>/dev/null 2>&1; then
                    # Mutation completed, safe to delete PDBs
                    find "$mutdir" -type f -name "*.pdb" -delete 2>/dev/null
                    rm -rf "$mutdir/molecules" 2>/dev/null
                fi
            done
        fi
    done
}

# Post-run cleanup function - aggressive cleanup after mutatex completes
# Usage: cleanup_mutations_dir <mutations_dir>
cleanup_mutations_dir() {
    local mutations_dir="$1"
    
    if [[ ! -d "$mutations_dir" ]]; then
        return 0
    fi
    
    # Delete all PDB files in mutation subdirectories (both mutant and WT structures)
    # These are regenerated per mutation and consume ~280MB × number_of_mutations
    find "$mutations_dir" -type f -name "*.pdb" -delete
    
    # Delete raw FoldX output files (keep only Average/Dif summaries)
    find "$mutations_dir" -type f \( \
        -name "Raw_*.fxout" \
        -o -name "PdbList_*.fxout" \
    \) -delete
    
    # Delete molecules subdirectories (temporary FoldX files)
    find "$mutations_dir" -type d -name "molecules" -exec rm -rf {} + 2>/dev/null || true
}

# Start cleanup worker in background and return PID
# Usage: start_cleanup_worker <mutations_dir> [interval]
#        cleanup_pid=$(start_cleanup_worker "$outdir/mutations" 30)
start_cleanup_worker() {
    local mutations_dir="$1"
    local interval="${2:-60}"
    
    cleanup_worker "$mutations_dir" "$interval" &
    echo $!
}

# Stop cleanup worker by PID
# Usage: stop_cleanup_worker $cleanup_pid
stop_cleanup_worker() {
    local pid="$1"
    
    if [[ -n "$pid" ]]; then
        kill "$pid" 2>/dev/null || true
        wait "$pid" 2>/dev/null || true
    fi
}
