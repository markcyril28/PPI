#!/usr/bin/env bash
#######################################################################
# GROMACS Common Functions and Configuration
#
# Source this file in other GROMACS scripts:
#   source "$(dirname "$0")/modules/gromacs_common.sh"
#
# Provides:
# - Common environment variables and paths
# - GPU configuration for AMD HIP
# - Logging functions
# - Structure preparation functions
# - GROMACS wrapper functions
#######################################################################

#------------------------------------------------------------------------------
# PATH CONFIGURATION
#------------------------------------------------------------------------------

# Detect script directory (where this file lives)
if [[ -z "${GROMACS_COMMON_DIR:-}" ]]; then
    GROMACS_COMMON_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Project root (parent of modules)
if [[ -z "${PROJECT_ROOT:-}" ]]; then
    PROJECT_ROOT="$(dirname "$GROMACS_COMMON_DIR")"
fi

# Modules directory
MODULES_DIR="${GROMACS_COMMON_DIR}/gromacs_utils"

# Add modules to Python path for CLI access
export PYTHONPATH="${GROMACS_COMMON_DIR}:${PYTHONPATH:-}"

#------------------------------------------------------------------------------
# DEFAULT GROMACS CONFIGURATION
#------------------------------------------------------------------------------

# Fix MPI library conflicts (use conda's OpenMPI, not system)
# Always prefer the gromacs_HIP environment for GROMACS binaries
GROMACS_HIP_ENV="/mnt/local3.5tb/home/mcmercado/miniconda3/envs/gromacs_HIP"

# Determine the best GROMACS environment to use
if [[ -z "${GMX_BIN:-}" ]]; then
    if [[ -x "${CONDA_PREFIX:-}/bin/gmx" ]]; then
        # Current conda env has gmx
        GMX_CONDA_PATH="${CONDA_PREFIX}"
    elif [[ -x "${GROMACS_HIP_ENV}/bin/gmx" ]]; then
        # Use gromacs_HIP environment (preferred for AMD GPU)
        GMX_CONDA_PATH="${GROMACS_HIP_ENV}"
    else
        # Fall back to current CONDA_PREFIX even if gmx isn't there (will error later)
        GMX_CONDA_PATH="${CONDA_PREFIX:-${GROMACS_HIP_ENV}}"
    fi
else
    # GMX_BIN already set, derive conda path from it
    GMX_CONDA_PATH="$(dirname "$(dirname "${GMX_BIN}")" 2>/dev/null)" || GMX_CONDA_PATH="${GROMACS_HIP_ENV}"
fi

# CRITICAL: Prepend conda lib path to avoid system MPI conflicts
# This fixes: "undefined symbol: ompi_mpi_errors_throw_exceptions"
export LD_LIBRARY_PATH="${GMX_CONDA_PATH}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# Use non-MPI gmx for single-node runs (avoids MPI library conflicts)
# For multi-node MPI runs, set GMX_BIN to gmx_mpi and ensure LD_PRELOAD is set
# Try conda path first, then fall back to system gmx
if [[ -z "${GMX_BIN:-}" ]]; then
    if [[ -x "${GMX_CONDA_PATH}/bin/gmx" ]]; then
        GMX_BIN="${GMX_CONDA_PATH}/bin/gmx"
    elif command -v gmx &>/dev/null; then
        GMX_BIN="$(command -v gmx)"
    else
        GMX_BIN="${GMX_CONDA_PATH}/bin/gmx"  # Will fail later with clear error
    fi
fi

# For MPI version, uncomment below and ensure conda MPI libs are preloaded:
# export LD_PRELOAD="${CONDA_PREFIX}/lib/libmpi.so${LD_PRELOAD:+:$LD_PRELOAD}"
# GMX_BIN="${GMX_BIN:-${CONDA_PREFIX}/bin/gmx_mpi}"

# Force field and water model
FORCEFIELD="${FORCEFIELD:-amber99sb-ildn}"
WATERMODEL="${WATERMODEL:-tip3p}"

# Box and solvation
BOX_DISTANCE="${BOX_DISTANCE:-1.0}"
BOX_TYPE="${BOX_TYPE:-dodecahedron}"
ION_CONCENTRATION="${ION_CONCENTRATION:-0.15}"

# Simulation parameters
EM_STEPS="${EM_STEPS:-5000}"
NVT_STEPS="${NVT_STEPS:-50000}"
NPT_STEPS="${NPT_STEPS:-50000}"
MD_STEPS="${MD_STEPS:-250000}"

#------------------------------------------------------------------------------
# GPU CONFIGURATION (AMD HIP)
#------------------------------------------------------------------------------

# Number of threads
NTHREADS="${NTHREADS:-8}"

# GPU flags for different simulation types
# EM: Can't use -pme gpu with steep integrator
GPU_EM_FLAGS="${GPU_EM_FLAGS:--nb gpu -pme cpu -bonded cpu}"
# MD: Maximum GPU offload (bonded cpu required for HIP backend)
GPU_MD_FLAGS="${GPU_MD_FLAGS:--nb gpu -pme gpu -bonded cpu -update gpu}"

# Environment variables for AMD GPU performance
setup_amd_gpu_env() {
    export GMX_ENABLE_DIRECT_GPU_COMM=1
    export GPU_MAX_HW_QUEUES=8
    export HIP_VISIBLE_DEVICES="${HIP_VISIBLE_DEVICES:-0}"
    # Note: Do NOT set HSA_OVERRIDE_GFX_VERSION - causes kernel arch mismatch
}

# Call setup by default
setup_amd_gpu_env

#------------------------------------------------------------------------------
# LOGGING FUNCTIONS
#------------------------------------------------------------------------------

# Simple log with timestamp
log() {
    echo "[$(date '+%H:%M:%S')] $1"
}

log_info() {
    echo "[$(date '+%H:%M:%S')] INFO: $1"
}

log_warn() {
    echo "[$(date '+%H:%M:%S')] WARN: $1" >&2
}

log_error() {
    echo "[$(date '+%H:%M:%S')] ERROR: $1" >&2
}

log_section() {
    echo ""
    echo "================================================================"
    echo " $1"
    echo "================================================================"
}

#------------------------------------------------------------------------------
# REQUIREMENT CHECKS
#------------------------------------------------------------------------------

check_gromacs() {
    if [[ ! -x "$GMX_BIN" ]]; then
        log_error "GROMACS not found at: $GMX_BIN"
        log_error "Please activate the gromacs_HIP conda environment: conda activate gromacs_HIP"
        return 1
    fi
    
    # Verify GROMACS can actually run (check for library issues)
    local version_output
    version_output=$("$GMX_BIN" --version 2>&1) || {
        log_error "GROMACS failed to run: $version_output"
        log_error "This may be a library path issue. Try: conda activate gromacs_HIP"
        return 1
    }
    
    log_info "GROMACS: $(echo "$version_output" | head -1)"
    return 0
}

check_file() {
    local file="$1"
    local desc="${2:-File}"
    if [[ ! -f "$file" ]]; then
        log_error "$desc not found: $file"
        return 1
    fi
    return 0
}

check_directory() {
    local dir="$1"
    local desc="${2:-Directory}"
    if [[ ! -d "$dir" ]]; then
        log_error "$desc not found: $dir"
        return 1
    fi
    return 0
}

check_python_modules() {
    if [[ ! -d "$MODULES_DIR" ]]; then
        log_error "Python modules not found at: $MODULES_DIR"
        return 1
    fi
    return 0
}

#------------------------------------------------------------------------------
# STRUCTURE PREPARATION
#------------------------------------------------------------------------------

# Prepare structure using Python CLI (handles CIF/PDB, cleaning)
prepare_structure() {
    local input="$1"
    local output="${2:-clean.pdb}"
    
    python3 -m gromacs_utils.cli prepare-structure "$input" -o "$output"
}

# Quick PDB cleaning (bash-only, for simple cases)
clean_pdb_quick() {
    local input="$1"
    local output="$2"
    
    grep -E "^(ATOM|TER|END)" "$input" | \
        sed 's/HSD/HIS/g; s/HSE/HIS/g; s/HSP/HIS/g' > "$output"
}

# Convert CIF to PDB using Python CLI
convert_cif_to_pdb() {
    local cif_file="$1"
    local pdb_file="$2"
    
    python3 -m gromacs_utils.cli prepare-structure "$cif_file" -o "$pdb_file"
}

#------------------------------------------------------------------------------
# MDP FILE GENERATION
#------------------------------------------------------------------------------

# Generate MDP files using Python CLI
generate_mdp() {
    local mdp_type="$1"
    local output="$2"
    shift 2
    
    python3 -m gromacs_utils.cli generate-mdp "$mdp_type" -o "$output" "$@"
}

# Generate all MDP files to a directory
generate_all_mdp() {
    local output_dir="$1"
    local em_steps="${2:-$EM_STEPS}"
    local nvt_steps="${3:-$NVT_STEPS}"
    local npt_steps="${4:-$NPT_STEPS}"
    local md_steps="${5:-$MD_STEPS}"
    
    python3 -m gromacs_utils.cli generate-mdp all -o "$output_dir" \
        --em-steps "$em_steps" \
        --nvt-steps "$nvt_steps" \
        --npt-steps "$npt_steps" \
        --md-steps "$md_steps"
}

#------------------------------------------------------------------------------
# GROMACS WRAPPER FUNCTIONS
#------------------------------------------------------------------------------

# Run mdrun with GPU, fall back to CPU if needed
run_mdrun() {
    local name="$1"
    local gpu_flags="$2"
    local log_dir="${3:-.}"
    
    mkdir -p "$log_dir"
    
    # GROMACS 2025+ requires -ntmpi when using GPU with OpenMP threads
    if $GMX_BIN mdrun -v -deffnm "$name" -ntmpi 1 -ntomp $NTHREADS $gpu_flags 2>&1 | tee "${log_dir}/mdrun_${name}.log"; then
        return 0
    else
        log_warn "GPU failed, using CPU..."
        $GMX_BIN mdrun -deffnm "$name" -ntmpi 1 -ntomp 32 2>&1 | tee "${log_dir}/mdrun_${name}_cpu.log"
    fi
}

# Run energy minimization
run_em() {
    local name="${1:-em}"
    local log_dir="${2:-.}"
    
    run_mdrun "$name" "$GPU_EM_FLAGS" "$log_dir"
}

# Run MD simulation (NVT, NPT, or production)
run_md() {
    local name="$1"
    local log_dir="${2:-.}"
    
    run_mdrun "$name" "$GPU_MD_FLAGS" "$log_dir"
}

#------------------------------------------------------------------------------
# CHAIN HANDLING
#------------------------------------------------------------------------------

# Create chain index using Python CLI
create_chain_index() {
    local pdb="$1"
    local gro="$2"
    local index="$3"
    
    python3 -m gromacs_utils.cli chain-index --pdb "$pdb" --gro "$gro" --index "$index"
}

# Get chain info
get_chain_info() {
    local pdb="$1"
    local outdir="${2:-.}"
    
    python3 -m gromacs_utils.cli chain-info "$pdb" --outdir "$outdir"
}

#------------------------------------------------------------------------------
# ANALYSIS HELPERS
#------------------------------------------------------------------------------

# Extract metrics from GROMACS outputs
extract_metrics() {
    local workdir="$1"
    local output="${2:-metrics.txt}"
    
    python3 -m gromacs_utils.cli extract-metrics --workdir "$workdir" -o "$output"
}

# Generate plots from analysis
generate_plots() {
    local analysis_dir="$1"
    local plots_dir="$2"
    
    python3 -m gromacs_utils.cli generate-plots md -i "$analysis_dir" -o "$plots_dir"
}

# Generate gnuplot scripts and run if available
run_gnuplot_scripts() {
    local plots_dir="$1"
    
    if command -v gnuplot &> /dev/null; then
        cd "$plots_dir"
        for gp in *.gp; do
            [[ -f "$gp" ]] && gnuplot "$gp" 2>/dev/null || true
        done
        cd - > /dev/null
    fi
}

#------------------------------------------------------------------------------
# OUTPUT ORGANIZATION
#------------------------------------------------------------------------------

# Create standard output structure
setup_output_dirs() {
    local base_dir="$1"
    
    mkdir -p "${base_dir}"/{logs,structures,analysis,trajectories,statistics,visualization,plots}
}

# Setup output structure using Python CLI
setup_output_structure() {
    local workdir="$1"
    
    python3 -m gromacs_utils.cli setup-output --workdir "$workdir"
}

#------------------------------------------------------------------------------
# VISUALIZATION
#------------------------------------------------------------------------------

# Generate visualization scripts
generate_visualization() {
    local workdir="$1"
    local viz_type="${2:-all}"
    
    python3 -m gromacs_utils.visualization_generator --workdir "$workdir" --type "$viz_type"
}

# Run batch analysis
run_batch_analysis() {
    local workdir="$1"
    local with_plots="${2:-true}"
    
    if [[ "$with_plots" == "true" ]]; then
        python3 -m gromacs_utils.cli batch-analysis --workdir "$workdir" --plots
    else
        python3 -m gromacs_utils.cli batch-analysis --workdir "$workdir"
    fi
}

#------------------------------------------------------------------------------
# MODULE INFO
#------------------------------------------------------------------------------

gromacs_common_info() {
    echo "GROMACS Common Configuration"
    echo "============================="
    echo "Project Root: $PROJECT_ROOT"
    echo "Modules Dir:  $MODULES_DIR"
    echo "GMX Binary:   $GMX_BIN"
    echo "Force Field:  $FORCEFIELD"
    echo "Water Model:  $WATERMODEL"
    echo "GPU EM Flags: $GPU_EM_FLAGS"
    echo "GPU MD Flags: $GPU_MD_FLAGS"
    echo "Threads:      $NTHREADS"
}

# Export for use in subshells
export GMX_BIN FORCEFIELD WATERMODEL BOX_DISTANCE BOX_TYPE ION_CONCENTRATION
export EM_STEPS NVT_STEPS NPT_STEPS MD_STEPS
export NTHREADS GPU_EM_FLAGS GPU_MD_FLAGS
export PROJECT_ROOT MODULES_DIR
