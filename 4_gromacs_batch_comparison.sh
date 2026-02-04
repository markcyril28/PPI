#!/usr/bin/env bash
set -e

#######################################################################
# GROMACS Batch PPI Comparison Script
#
# Analyzes ALL structures in input folder and ranks them.
# Uses Python modules for analysis and reporting.
#
# Usage: ./gromacs_batch_comparison.sh [DATASET] [OUTPUT_BASE_DIR]
#   DATASET: "SmelDMP", "SmelGRF-GIF", or "all" (default: all)
#   OUTPUT_BASE_DIR: Base directory for results (default: results_gromacs)
#######################################################################

#------------------------------------------------------------------------------
# CONFIGURATION
#------------------------------------------------------------------------------

# Base paths
BASE_INPUT_DIR="/mnt/local3.5tb/home/mcmercado/PPI/inputs"
BASE_OUTPUT_DIR="${2:-/mnt/local3.5tb/home/mcmercado/PPI/results_gromacs}"

# Override existing results (true) or skip if already processed (false)
OVERRIDE_EXISTING=true

# Dataset selection (SmelDMP, SmelGRF-GIF, or all)
DATASET="${1:-SmelDMP}"
#DATASET="${1:-SmelGRF-GIF}"

# Available datasets
SMELDMP_INPUT="${BASE_INPUT_DIR}/SmelDMP"
SMELGRF_GIF_INPUT="${BASE_INPUT_DIR}/SmelGRF-GIF"

# GROMACS settings (can override defaults from gromacs_common.sh)
EM_STEPS=2000
BOX_DISTANCE=0.8
BOX_TYPE=cubic

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/modules/gromacs_common.sh"

#------------------------------------------------------------------------------
# FUNCTIONS
#------------------------------------------------------------------------------

analyze_structure() {
    local input_file="$1"
    local struct_name="$2"
    local struct_dir="$3"
    
    # Skip if already processed (unless OVERRIDE_EXISTING is true)
    if [[ "$OVERRIDE_EXISTING" != "true" && -f "$struct_dir/metrics.json" && -f "$struct_dir/em.gro" ]]; then
        log "  ✓ Already processed, skipping (found metrics.json and em.gro)"
        log "    Set OVERRIDE_EXISTING=true to reprocess"
        echo "OK" > "$struct_dir/status.txt"
        return 0
    fi
    
    # If overriding, clean up old files
    if [[ "$OVERRIDE_EXISTING" == "true" && -d "$struct_dir" ]]; then
        log "  Overriding existing results..."
        rm -rf "$struct_dir"
    fi
    
    mkdir -p "$struct_dir"
    cd "$struct_dir"
    
    log "  Preparing structure..."
    
    # Use Python CLI module for structure conversion and cleaning
    python3 -m gromacs_utils.cli prepare-structure "$input_file" -o clean.pdb
    
    # Generate topology
    echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh > pdb2gmx.log 2>&1 || {
        echo "ERROR" > status.txt
        return 1
    }
    
    # Create index
    echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx > make_ndx.log 2>&1
    
    # Setup box and solvate
    $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt $BOX_TYPE > editconf.log 2>&1
    $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top > solvate.log 2>&1
    
    # Generate EM MDP using Python CLI module
    python3 -m gromacs_utils.cli generate-mdp em -o em.mdp --em-steps $EM_STEPS --em-tolerance 500.0
    
    # Run energy minimization
    $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 3 > grompp.log 2>&1
    run_em em .
    
    # Extract metrics using GROMACS
    echo -e "Potential\n\n" | $GMX_BIN energy -f em.edr -o energy.xvg 2>&1 > energy.log || true
    echo "Protein" | $GMX_BIN gyrate -s em.tpr -f em.gro -o gyrate.xvg 2>&1 > gyrate.log || true
    echo "Protein" | $GMX_BIN sasa -s em.tpr -f em.gro -o sasa.xvg 2>&1 > sasa.log || true
    
    # Parse and save metrics using Python CLI module
    python3 -m gromacs_utils.cli extract-metrics --workdir . -o metrics.json
    
    echo "OK" > status.txt
}

#------------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------------

process_dataset() {
    local input_dir="$1"
    local dataset_name="$2"
    local output_dir="${BASE_OUTPUT_DIR}/4_batch_comparison/${dataset_name}"
    
    log_section "Processing Dataset: $dataset_name"
    log "Input directory: $input_dir"
    log "Output directory: $output_dir"
    
    if [[ ! -d "$input_dir" ]]; then
        log "WARNING: Input directory does not exist: $input_dir"
        return 1
    fi
    
    mkdir -p "$output_dir"
    cd "$output_dir"
    
    # Find all PDB/CIF files in this dataset
    mapfile -t STRUCTURES < <(find "$input_dir" -maxdepth 1 -type f \( -name "*.pdb" -o -name "*.cif" \) | sort)
    
    if [[ ${#STRUCTURES[@]} -eq 0 ]]; then
        log "No structures found in $input_dir"
        return 1
    fi
    
    log "Found ${#STRUCTURES[@]} structures to analyze"
    echo ""
    
    local dataset_start=$(date +%s)
    
    # Analyze each structure
    for struct_file in "${STRUCTURES[@]}"; do
        struct_name=$(basename "$struct_file" | sed 's/\.[^.]*$//')
        struct_dir="$output_dir/$struct_name"
        
        log "Analyzing: $struct_name"
        
        if analyze_structure "$struct_file" "$struct_name" "$struct_dir"; then
            log "  ✓ Complete"
        else
            log "  ✗ Failed"
        fi
        echo ""
    done
    
    # Generate comparison report and plots using Python modules
    log "Generating comparison report for $dataset_name..."
    cd "$output_dir"
    
    python3 -m gromacs_utils.cli batch-analysis --workdir "$output_dir" --plots
    
    # Run gnuplot scripts if available
    run_gnuplot_scripts "$output_dir/plots"
    
    local dataset_end=$(date +%s)
    
    log_section "$dataset_name Analysis Complete"
    log "Time: $((dataset_end - dataset_start)) seconds"
    log "Structures analyzed: ${#STRUCTURES[@]}"
    log "Output: $output_dir"
    log ""
    log "Key files:"
    log "  - comparison_report.txt"
    log "  - comparison_data.csv"
    log "  - comparison_data.json"
    log "  - plots/"
    
    # Show top results
    echo ""
    head -15 "$output_dir/comparison_report.txt" 2>/dev/null || true
}

main() {
    log_section "GROMACS Batch PPI Comparison"
    log "Dataset selection: $DATASET"
    echo ""
    
    START=$(date +%s)
    
    case "$DATASET" in
        SmelDMP|smeldmp)
            process_dataset "$SMELDMP_INPUT" "SmelDMP"
            ;;
        SmelGRF-GIF|smelgrf-gif|SmelGRF_GIF|smelgrf_gif)
            process_dataset "$SMELGRF_GIF_INPUT" "SmelGRF-GIF"
            ;;
        all|ALL)
            log "Processing all datasets..."
            echo ""
            process_dataset "$SMELDMP_INPUT" "SmelDMP"
            echo ""
            echo "=========================================="
            echo ""
            process_dataset "$SMELGRF_GIF_INPUT" "SmelGRF-GIF"
            ;;
        *)
            echo "ERROR: Unknown dataset '$DATASET'"
            echo "Usage: $0 [DATASET] [OUTPUT_BASE_DIR]"
            echo "  DATASET: SmelDMP, SmelGRF-GIF, or all (default: all)"
            exit 1
            ;;
    esac
    
    END=$(date +%s)
    
    log_section "All Processing Complete"
    log "Total time: $((END - START)) seconds"
    log "Results saved to: $BASE_OUTPUT_DIR"
}

main "$@"
