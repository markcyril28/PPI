#!/usr/bin/env bash
#set -euo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source logging utilities
source "$SCRIPT_DIR/modules/logging/logging_utils.sh"
source "$SCRIPT_DIR/modules/logging/gpu_utils.sh"
source "$SCRIPT_DIR/modules/mutatex_utils.sh"

# Set up logging with absolute paths (needed because script changes directories)
RUN_ID="mutatex_$(date +%Y%m%d_%H%M%S)"
LOG_DIR="$SCRIPT_DIR/logs/log_files"
TIME_DIR="$SCRIPT_DIR/logs/time_logs"
SPACE_DIR="$SCRIPT_DIR/logs/space_logs"
SPACE_TIME_DIR="$SCRIPT_DIR/logs/space_time_logs"
ERROR_WARN_DIR="$SCRIPT_DIR/logs/error_warn_logs"
SOFTWARE_CATALOG_DIR="$SCRIPT_DIR/logs/software_catalogs"
GPU_LOG_DIR="$SCRIPT_DIR/logs/gpu_log"
setup_logging

log_step "Starting MutateX batch processing"

# Activate the conda environment
CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"
ENV_NAME="${ENV_NAME:-PPI}"
conda activate "$ENV_NAME"

# Source mutatex environment hints if available
ENV_HINT_FILE="${ENV_HINT_FILE:-$SCRIPT_DIR/mutatex.env}"
if [[ -f "$ENV_HINT_FILE" ]]; then
  log_info "Loading environment hints from $ENV_HINT_FILE"
  source "$ENV_HINT_FILE"
fi

# Configuration variables
THREADS=32
NRUNS=5
FORCE_RERUN="${FORCE_RERUN:-true}"  # Set to true to force re-run even if completed
CLEANUP_DURING_RUN="${CLEANUP_DURING_RUN:-true}"  # Set to true for incremental cleanup (saves ~195GB)

PDB_FILES=(
  #"$SCRIPT_DIR/inputs/SmelGRF-GIF/fold_1_x_1_model_0.pdb"
  #"$SCRIPT_DIR/inputs/SmelGRF-GIF/fold_2_x_1_model_0.pdb"
  #"$SCRIPT_DIR/inputs/SmelGRF-GIF/fold_recruitment_of_smelgif_to_swi_model_0.pdb"
  
  # SmelDMP structures (new shorter filenames)
  #"$SCRIPT_DIR/inputs/SmelDMP/SmelDMP01.730_SmelHAP2.pdb"
  #"$SCRIPT_DIR/inputs/SmelDMP/SmelDMP01.990_SmelHAP2.pdb"
  #"$SCRIPT_DIR/inputs/SmelDMP/SmelDMP02_SmelHAP2.pdb"
  #"$SCRIPT_DIR/inputs/SmelDMP/SmelDMP04_SmelHAP2.pdb"
  "$SCRIPT_DIR/inputs/SmelDMP/SmelDMP10.200_SmelHAP2.pdb"
  #"$SCRIPT_DIR/inputs/SmelDMP/SmelDMP10_550_560_SmelHAP2.pdb"
  #"$SCRIPT_DIR/inputs/SmelDMP/SmelDMP12_SmelHAP2.pdb"
)

# FoldX configuration (use local modules)
FOLDX_BINARY="$SCRIPT_DIR/modules/foldx_20270131_5.1"
ROTABASE="$SCRIPT_DIR/modules/rotabase.txt"

# Templates (use local modules)
REPAIR_TEMPLATE="$SCRIPT_DIR/modules/repair_runfile_template.txt"
MUTATE_TEMPLATE="$SCRIPT_DIR/modules/mutate_runfile_template.txt"
INTERFACE_TEMPLATE="$SCRIPT_DIR/modules/interface_runfile_template.txt"

# Absolute path for log file
LOG_FILE_ABS="$SCRIPT_DIR/$LOG_FILE"

log_info "Using FoldX: $FOLDX_BINARY"
log_info "Using templates from: $SCRIPT_DIR/modules/"
log_info "Threads: $THREADS, Runs: $NRUNS"

# Loop through each PDB file
for pdb in "${PDB_FILES[@]}"; do
  # Get the base name without path and extension for output directory
  base_name=$(basename "$pdb" .pdb)
  outdir="$SCRIPT_DIR/results_mutatex/${base_name}"
  completion_marker="$outdir/.completed"
  
  # Check if run was already completed successfully
  if [[ -f "$completion_marker" && "$FORCE_RERUN" != "true" ]]; then
    log_warn "Skipping: $pdb (already completed - marker found at $completion_marker)"
    echo ""
    continue
  fi
  
  # If directory exists but no completion marker, it's a partial/cancelled run - continue from where we left off
  if [[ -d "$outdir" && ! -f "$completion_marker" ]]; then
    log_info "Resuming interrupted run for: $pdb"
  fi
  
  log_step "Processing: $pdb"
  
  # Change to output directory (mutatex runs in current directory)
  mkdir -p "$outdir"
  pushd "$outdir" > /dev/null
  
  # Start background cleanup worker if enabled (prevents 196GB disk usage)
  cleanup_pid=""
  if [[ "$CLEANUP_DURING_RUN" == "true" ]]; then
    log_info "Starting incremental cleanup worker (saves disk space during run)"
    cleanup_worker "$outdir/mutations" 30 &
    cleanup_pid=$!
  fi
  
  mutatex "$pdb" \
    -p "$THREADS" \
    -n "$NRUNS" \
    -x "$FOLDX_BINARY" \
    -f suite5 \
    -R "$REPAIR_TEMPLATE" \
    -M "$MUTATE_TEMPLATE" \
    -I "$INTERFACE_TEMPLATE" \
    -B -v \
    -l 2>&1 | tee -a "$LOG_FILE_ABS"
  
  # Check exit status
  exit_status=${PIPESTATUS[0]}
  
  # Stop the cleanup worker
  if [[ -n "$cleanup_pid" ]]; then
    kill "$cleanup_pid" 2>/dev/null || true
    wait "$cleanup_pid" 2>/dev/null || true
    log_info "Stopped incremental cleanup worker"
  fi
  
  if [[ $exit_status -eq 0 ]]; then
    log_info "MutateX completed successfully for $pdb"
    # Generate PyMOL visualization
    mutatex plot --pymol
    # Create completion marker
    touch "$completion_marker"
    log_info "Created completion marker: $completion_marker"
  else
    log_error "MutateX failed for $pdb"
    # Remove completion marker if it exists (in case of force re-run failure)
    rm -f "$completion_marker"
  fi
  
  popd > /dev/null
  
  # Aggressive cleanup to save disk space (196GB+ per run without this)
  # Delete intermediate PDB files but keep essential results
  if [[ -d "$outdir/mutations" ]]; then
    log_info "Cleaning up intermediate files to save disk space..."
    cleanup_mutations_dir "$outdir/mutations"
    log_info "Cleanup complete. Space saved."
  fi
    
  log_step "Completed: $pdb"
  echo ""
done

log_step "Batch processing complete"

python "$SCRIPT_DIR/mutatex_extract_all_critical_residues.py"

python "$SCRIPT_DIR/mutatex_generate_mutatex_visualizations.py"

log_info "Full log saved to: $LOG_FILE_ABS"
