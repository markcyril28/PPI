#!/usr/bin/env bash
set -e

#######################################################################
# GROMACS Protein-Protein Interaction Stability Comparison Script
# GPU-ACCELERATED VERSION
#
# Compares the stability of chain pairings (protein-protein interactions)
# between two PDB structures using GROMACS with GPU acceleration.
#
# GPU NOTES:
# - Uses GPU acceleration for MD simulations
# - Includes GPU detection and fallback to CPU if unavailable
# - For AMD GPUs: Use GROMACS built with native HIP (not SYCL/AdaptiveCpp)
#
# Performs:
# 1. Structure preparation and topology generation
# 2. Energy minimization (GPU-accelerated with -nb gpu)
# 3. Short equilibration (NVT and NPT) with full GPU offload
# 4. Short production MD with full GPU offload
# 5. Interaction energy calculation between chains
# 6. Comparative analysis
#
# Usage: ./gromacs_compare_chain_stability_gpu.sh [OPTIONS] [PDB1] [PDB2] [WORKDIR]
#
# Options:
#   -m, --mode MODE    Run mode: full (default), sim-only, viz-only, compare-only
#   -h, --help         Show this help message
#
# Modes:
#   full         - Run simulation + analysis + visualization (default)
#   sim-only     - Run simulation and analysis only (no plots/viz)
#   viz-only     - Generate visualization only (requires existing data)
#   compare-only - Run comparison on existing structures (skip completed)
#
# Examples:
#   # Using default PDB files and output directory (defined in script)
#   ./gromacs_compare_chain_stability_gpu.sh
#   ./gromacs_compare_chain_stability_gpu.sh --mode sim-only
#   ./gromacs_compare_chain_stability_gpu.sh --mode viz-only
#
#   # With custom PDB files and output directory
#   ./gromacs_compare_chain_stability_gpu.sh PDB1 PDB2 OUTPUT_DIR
#   ./gromacs_compare_chain_stability_gpu.sh --mode sim-only PDB1 PDB2 OUTPUT_DIR
#   ./gromacs_compare_chain_stability_gpu.sh --mode viz-only PDB1 PDB2 OUTPUT_DIR
#
# Requirements: GROMACS 2020+ built with GPU support (CUDA/HIP/SYCL)
#######################################################################

#------------------------------------------------------------------------------
# ARGUMENT PARSING
#------------------------------------------------------------------------------

show_help() {
    sed -n '3,38p' "$0" | sed 's/^# //' | sed 's/^#//'
    exit 0
}

# Default mode
MODE="full"

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--mode)
            MODE="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        -*)
            echo "Unknown option: $1"
            show_help
            ;;
        *)
            break
            ;;
    esac
done

# Validate mode
case $MODE in
    full|sim-only|viz-only|compare-only) ;;
    *)
        echo "ERROR: Invalid mode '$MODE'. Use: full, sim-only, viz-only, or compare-only"
        exit 1
        ;;
esac

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------

INPUT_BASE="/mnt/local3.5tb/home/mcmercado/PPI/inputs"

# (minimum 2 structures required)
PDB_LIST=(
    # SmelDMP files:
    "${INPUT_BASE}/SmelDMP/SmelDMP01.730_SmelHAP2.pdb"
    "${INPUT_BASE}/SmelDMP/SmelDMP01.990_SmelHAP2.pdb"
    "${INPUT_BASE}/SmelDMP/SmelDMP02_SmelHAP2.pdb"
    "${INPUT_BASE}/SmelDMP/SmelDMP04_SmelHAP2.pdb"
    "${INPUT_BASE}/SmelDMP/SmelDMP10.200_SmelHAP2.pdb"
    "${INPUT_BASE}/SmelDMP/SmelDMP10_550_560_SmelHAP2.pdb"
    "${INPUT_BASE}/SmelDMP/SmelDMP12_SmelHAP2.pdb"
)

# Output base directory
OUTPUT_BASE="${1:-/mnt/local3.5tb/home/mcmercado/PPI/results_gromacs}"

# Generate dynamic output folder name based on structures being compared
generate_workdir_name() {
    local names=()
    for pdb in "${PDB_LIST[@]}"; do
        # Extract base name without extension and path
        local basename=$(basename "$pdb" .pdb)
        # Shorten the name for readability (take first meaningful part)
        local shortname=$(echo "$basename" | sed 's/_model_0.*//; s/_2025_01_09.*//; s/fold_//')
        names+=("$shortname")
    done
    # Join names with "_vs_"
    local joined=$(IFS="_vs_"; echo "${names[*]}")
    echo "${OUTPUT_BASE}/2_compare/${joined}"
}

# Set WORKDIR dynamically
WORKDIR=$(generate_workdir_name)

# Validate we have at least 2 structures
if [ ${#PDB_LIST[@]} -lt 2 ]; then
    echo "ERROR: At least 2 structures must be uncommented in PDB_LIST"
    exit 1
fi

# Simulation parameters (override defaults before sourcing)
BOX_DISTANCE=1.2
ION_CONCENTRATION=0.15
EM_STEPS=5000
NVT_STEPS=50000     # 100 ps
NPT_STEPS=50000     # 100 ps
MD_STEPS=250000     # 500 ps

# GPU Configuration
GPU_ID=0

# Override existing results (true) or skip if already processed (false)
OVERRIDE_EXISTING=true

# Maximum threads to use for GROMACS
MAX_THREADS=32
NTHREADS=$MAX_THREADS

# Source common functions and configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/modules/gromacs_common.sh"

# GPU mode flag
USE_GPU=true
GPU_INFO=""

#------------------------------------------------------------------------------
# GPU Functions
#------------------------------------------------------------------------------

check_gpu() {
    log "Checking GPU availability..."
    
    # Check for AMD GPU with rocm-smi
    if command -v rocm-smi &> /dev/null; then
        GPU_INFO=$(rocm-smi --showproductname 2>/dev/null | grep -i "GPU" | head -1 || echo "")
        if [ -n "$GPU_INFO" ]; then
            log "  Found AMD GPU: $GPU_INFO"
            return 0
        fi
    fi
    
    # Check for NVIDIA GPU
    if command -v nvidia-smi &> /dev/null; then
        GPU_INFO=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || echo "")
        if [ -n "$GPU_INFO" ]; then
            log "  Found NVIDIA GPU: $GPU_INFO"
            return 0
        fi
    fi
    
    log_warn "No GPU detected, will use CPU mode"
    USE_GPU=false
    return 1
}

build_mdrun_cmd() {
    local stage="$1"  # em, nvt, npt, md
    local deffnm="$2"
    
    # GROMACS 2025+ requires -ntmpi when using GPU with OpenMP threads
    local cmd="$GMX_BIN mdrun -v -deffnm $deffnm -ntmpi 1 -ntomp $NTHREADS"
    
    if [ "$USE_GPU" = true ]; then
        case "$stage" in
            em)
                cmd="$cmd -gpu_id $GPU_ID $GPU_EM_FLAGS"
                ;;
            nvt|npt|md)
                cmd="$cmd -gpu_id $GPU_ID $GPU_MD_FLAGS"
                ;;
        esac
    fi
    
    echo "$cmd"
}

#------------------------------------------------------------------------------
# Stage Completion Check Functions
#------------------------------------------------------------------------------

# Check if simulation is complete (has MD output files)
check_simulation_complete() {
    local outdir="$1"
    [[ -f "$outdir/md.gro" && -f "$outdir/md.edr" ]]
}

# Check if EM is complete (minimum for analysis)
check_em_complete() {
    local outdir="$1"
    [[ -f "$outdir/em.gro" && -f "$outdir/em.edr" ]]
}

# Check if interaction energy is calculated
check_ie_complete() {
    local outdir="$1"
    [[ -f "$outdir/ie.edr" || -f "$outdir/interaction_energy.xvg" ]]
}

# Check if trajectory outputs exist
check_trajectory_complete() {
    local outdir="$1"
    [[ -f "$outdir/statistics/md_statistics.json" ]]
}

# Check if visualization is complete
check_visualization_complete() {
    local outdir="$1"
    [[ -f "$outdir/visualization/visualize_trajectory.pml" || -f "$outdir/visualization/visualize_interface.pml" ]]
}

# Get structure completion status
get_structure_status() {
    local outdir="$1"
    local status="incomplete"
    
    if check_simulation_complete "$outdir"; then
        status="simulation_complete"
    elif check_em_complete "$outdir"; then
        status="em_complete"
    fi
    
    if check_trajectory_complete "$outdir"; then
        status="${status}+stats"
    fi
    
    if check_visualization_complete "$outdir"; then
        status="${status}+viz"
    fi
    
    echo "$status"
}

#------------------------------------------------------------------------------
# Processing Functions
#------------------------------------------------------------------------------

process_structure() {
    local pdb_input="$1"
    local name="$2"
    local outdir="$3"
    
    log_section "Processing: $name"
    
    # Check current status
    local status=$(get_structure_status "$outdir")
    log "  Current status: $status"
    
    # Skip if simulation is already complete
    if check_simulation_complete "$outdir"; then
        log "  ✓ Simulation complete, skipping (found md.gro and md.edr)"
        return 0
    fi
    
    setup_output_dirs "$outdir"
    pushd "$outdir" > /dev/null
    
    # Clean PDB using Python CLI module
    log "Cleaning PDB file..."
    python3 -m gromacs_utils.cli prepare-structure "$pdb_input" -o clean.pdb
    
    # Generate topology
    log "Generating topology with pdb2gmx..."
    echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh 2>&1 | tee logs/pdb2gmx.log
    
    # Create index file with chain selections
    log "Creating index file with chain groups..."
    echo "q" | $GMX_BIN make_ndx -f protein.gro -o index.ndx 2>&1 | tee logs/make_ndx.log
    python3 -m gromacs_utils.cli chain-index --pdb clean.pdb --gro protein.gro --index index.ndx
    
    # Define box
    log "Defining simulation box..."
    $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron 2>&1 | tee logs/editconf.log
    
    # Solvate
    log "Solvating the system..."
    $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 | tee logs/solvate.log
    
    # Generate all MDP files using Python CLI module
    log "Generating MDP files..."
    python3 -m gromacs_utils.cli generate-mdp all -o . \
        --em-steps $EM_STEPS \
        --nvt-steps $NVT_STEPS \
        --npt-steps $NPT_STEPS \
        --md-steps $MD_STEPS
    
    # Add ions
    log "Adding ions..."
    $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 2>&1 | tee logs/grompp_ions.log
    echo "SOL" | $GMX_BIN genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc $ION_CONCENTRATION 2>&1 | tee logs/genion.log
    
    # Energy minimization
    log "Running energy minimization..."
    $GMX_BIN grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr 2>&1 | tee logs/grompp_em.log
    local gmx_cmd=$(build_mdrun_cmd "em" "em")
    $gmx_cmd 2>&1 | tee logs/mdrun_em.log
    echo "Potential" | $GMX_BIN energy -f em.edr -o em_potential.xvg 2>&1 | tee logs/energy_em.log
    
    # NVT equilibration
    log "Running NVT equilibration..."
    $GMX_BIN grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr 2>&1 | tee logs/grompp_nvt.log
    gmx_cmd=$(build_mdrun_cmd "nvt" "nvt")
    $gmx_cmd 2>&1 | tee logs/mdrun_nvt.log
    
    # NPT equilibration  
    log "Running NPT equilibration..."
    $GMX_BIN grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr 2>&1 | tee logs/grompp_npt.log
    gmx_cmd=$(build_mdrun_cmd "npt" "npt")
    $gmx_cmd 2>&1 | tee logs/mdrun_npt.log
    
    # Production MD
    log "Running production MD..."
    $GMX_BIN grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr 2>&1 | tee logs/grompp_md.log
    gmx_cmd=$(build_mdrun_cmd "md" "md")
    $gmx_cmd 2>&1 | tee logs/mdrun_md.log
    
    popd > /dev/null
    log "✓ MD simulation completed for $name"
}

calculate_interaction_energy() {
    local outdir="$1"
    local name="$2"
    
    # Skip if already calculated
    if check_ie_complete "$outdir"; then
        log "  ✓ Interaction energy already calculated, skipping"
        return 0
    fi
    
    # Need at least EM to calculate IE
    if ! check_em_complete "$outdir"; then
        log_warn "  ⚠ No EM data found, skipping IE calculation"
        return 1
    fi
    
    log "Calculating Interaction Energy: $name"
    pushd "$outdir" > /dev/null
    
    # Generate IE MDP if not exists
    [ -f ie.mdp ] || python3 -m gromacs_utils.cli generate-mdp ie -o ie.mdp
    
    # Calculate interaction energy using rerun
    $GMX_BIN grompp -f ie.mdp -c md.gro -p topol.top -n index.ndx -o ie.tpr 2>&1 | tee logs/grompp_ie.log || true
    $GMX_BIN mdrun -s ie.tpr -rerun md.xtc -e ie.edr 2>&1 | tee logs/mdrun_ie.log || true
    
    # Extract interaction energies
    echo -e "Coul-SR:ChainA-ChainB\nLJ-SR:ChainA-ChainB\n\n" | $GMX_BIN energy -f ie.edr -o interaction_energy.xvg 2>&1 | tee logs/energy_ie.log || true
    
    popd > /dev/null
}

calculate_binding_metrics() {
    local outdir="$1"
    local name="$2"
    
    log "Calculating binding metrics: $name"
    cd "$outdir"
    
    # Use Python modules for analysis
    python3 -m gromacs_utils.md_statistics --workdir "$outdir"
}

generate_trajectory_outputs() {
    local outdir="$1"
    local name="$2"
    
    # Skip if already generated
    if check_trajectory_complete "$outdir"; then
        log "  ✓ Trajectory outputs already generated, skipping"
        return 0
    fi
    
    # Check what data is available
    local has_md=false
    local has_em=false
    [[ -f "$outdir/md.tpr" && -f "$outdir/md.xtc" ]] && has_md=true
    check_em_complete "$outdir" && has_em=true
    
    if [[ "$has_md" != true && "$has_em" != true ]]; then
        log_warn "  ⚠ No simulation data found, skipping trajectory outputs"
        return 1
    fi
    
    log "Generating trajectory outputs: $name"
    pushd "$outdir" > /dev/null
    
    # Center and fit trajectory
    echo "Protein Protein" | $GMX_BIN trjconv -s md.tpr -f md.xtc -o trajectories/md_center.xtc -center -pbc mol 2>&1 | tee logs/trjconv.log || true
    echo "Protein" | $GMX_BIN trjconv -s md.tpr -f md.gro -o structures/final_structure.pdb 2>&1 >> logs/trjconv.log || true
    
    # Analysis
    echo "Backbone Backbone" | $GMX_BIN rms -s md.tpr -f trajectories/md_center.xtc -o analysis/rmsd.xvg -tu ps 2>&1 > logs/rmsd.log || true
    echo "Backbone" | $GMX_BIN rmsf -s md.tpr -f trajectories/md_center.xtc -o analysis/rmsf.xvg -res 2>&1 > logs/rmsf.log || true
    echo "Protein" | $GMX_BIN gyrate -s md.tpr -f trajectories/md_center.xtc -o analysis/gyrate.xvg 2>&1 > logs/gyrate.log || true
    echo "Protein" | $GMX_BIN sasa -s md.tpr -f trajectories/md_center.xtc -o analysis/sasa.xvg 2>&1 > logs/sasa.log || true
    echo -e "ChainA\nChainB" | $GMX_BIN hbond -s md.tpr -f trajectories/md_center.xtc -n index.ndx -num analysis/hbonds.xvg 2>&1 > logs/hbond.log || true
    echo -e "ChainA\nChainB" | $GMX_BIN mindist -s md.tpr -f trajectories/md_center.xtc -n index.ndx -od analysis/mindist.xvg 2>&1 > logs/mindist.log || true
    
    # Generate statistics
    python3 -m gromacs_utils.md_statistics --workdir "$outdir"
    
    popd > /dev/null
}

generate_structure_visualization() {
    local outdir="$1"
    local name="$2"
    
    # Skip if already generated
    if check_visualization_complete "$outdir"; then
        log "  ✓ Visualization already generated, skipping"
        return 0
    fi
    
    log "Generating visualization for $name..."
    pushd "$outdir" > /dev/null
    
    # Check if analysis data exists
    if [[ ! -d "analysis" ]]; then
        log_warn "  No analysis directory found, skipping visualization"
        popd > /dev/null
        return 1
    fi
    
    # Generate visualization scripts
    log "  Generating visualization scripts..."
    generate_visualization "$outdir" all
    
    # Generate plots
    log "  Generating plots..."
    run_gnuplot_scripts plots
    python3 -m gromacs_utils.cli generate-plots md -i "$outdir/analysis" -o "$outdir/plots" 2>&1 || true
    
    popd > /dev/null
    log "  ✓ $name visualization complete"
}

compare_results() {
    log_section "Comparing Stability Results"
    
    local num_structs=${#PDB_LIST[@]}
    
    # Generate summary of all structures
    log "Generating summary for $num_structs structures..."
    
    # Create a combined summary file
    local summary_file="$WORKDIR/structures_summary.txt"
    {
        echo "=============================================="
        echo "MULTI-STRUCTURE STABILITY COMPARISON SUMMARY"
        echo "=============================================="
        echo ""
        echo "Structures analyzed: $num_structs"
        echo "Generated: $(date)"
        echo ""
        
        for i in "${!PDB_LIST[@]}"; do
            local struct_num=$((i+1))
            local struct_name="structure_${struct_num}"
            local struct_dir="$WORKDIR/$struct_name"
            local pdb_name=$(basename "${PDB_LIST[$i]}" .pdb)
            
            echo "--- Structure $struct_num: $pdb_name ---"
            
            if [[ -f "$struct_dir/statistics/md_statistics.json" ]]; then
                python3 -c "
import json
with open('$struct_dir/statistics/md_statistics.json') as f:
    data = json.load(f)
for key, val in data.items():
    if isinstance(val, (int, float)):
        print(f'  {key}: {val:.4f}' if isinstance(val, float) else f'  {key}: {val}')
    elif isinstance(val, dict):
        for k, v in val.items():
            if isinstance(v, (int, float)):
                print(f'  {key}.{k}: {v:.4f}' if isinstance(v, float) else f'  {key}.{k}: {v}')
" 2>/dev/null || echo "  (statistics unavailable)"
            else
                echo "  (no statistics file found)"
            fi
            echo ""
        done
    } > "$summary_file"
    
    log "Summary saved to: $summary_file"
    
    # Run pairwise comparisons for all structure pairs
    if [ $num_structs -ge 2 ]; then
        log "Running pairwise comparisons..."
        
        for ((i=0; i<num_structs-1; i++)); do
            for ((j=i+1; j<num_structs; j++)); do
                local struct1_num=$((i+1))
                local struct2_num=$((j+1))
                local struct1_name="structure_${struct1_num}"
                local struct2_name="structure_${struct2_num}"
                local struct1_dir="$WORKDIR/$struct1_name"
                local struct2_dir="$WORKDIR/$struct2_name"
                
                # Check if both structures have data
                if [[ -d "$struct1_dir" && -d "$struct2_dir" ]]; then
                    log "  Comparing $struct1_name vs $struct2_name..."
                    
                    python3 -m gromacs_utils.md_results_comparator \
                        --workdir "$WORKDIR" \
                        --pdb1 "${PDB_LIST[$i]}" \
                        --pdb2 "${PDB_LIST[$j]}" \
                        --struct1-dir "$struct1_name" \
                        --struct2-dir "$struct2_name" \
                        --output "comparison_${struct1_num}_vs_${struct2_num}.txt" 2>&1 || {
                            log_warn "  Comparison $struct1_name vs $struct2_name failed (may be missing data)"
                        }
                else
                    log_warn "  Skipping $struct1_name vs $struct2_name (missing data)"
                fi
            done
        done
        
        # Also generate the default comparison report (structure_1 vs structure_2)
        if [[ -d "$WORKDIR/structure_1" && -d "$WORKDIR/structure_2" ]]; then
            python3 -m gromacs_utils.md_results_comparator \
                --workdir "$WORKDIR" \
                --pdb1 "${PDB_LIST[0]}" \
                --pdb2 "${PDB_LIST[1]}" \
                --struct1-dir "structure_1" \
                --struct2-dir "structure_2" 2>&1 || true
        fi
        
        # Generate combined multi-structure comparison (radar plot, ranking) for 3+ structures
        if [ $num_structs -gt 2 ]; then
            log "Generating combined multi-structure comparison plots..."
            
            # Build PDB list and structure directory arguments
            local pdb_args=""
            local struct_args=""
            for i in "${!PDB_LIST[@]}"; do
                pdb_args="$pdb_args ${PDB_LIST[$i]}"
                struct_args="$struct_args structure_$((i+1))"
            done
            
            python3 -m gromacs_utils.md_results_comparator \
                --workdir "$WORKDIR" \
                --multi \
                --pdb-list $pdb_args \
                --struct-dirs $struct_args 2>&1 || true
        fi
    fi
    
    log "Comparison complete. Reports saved to $WORKDIR/"
}

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

main() {
    log_section "GROMACS Chain Stability Comparison (GPU-Accelerated)"
    
    log "Mode: $MODE"
    log "Structures to compare: ${#PDB_LIST[@]}"
    for i in "${!PDB_LIST[@]}"; do
        log "  Structure $((i+1)): $(basename ${PDB_LIST[$i]})"
    done
    log "Output: $WORKDIR"
    
    # Create working directory
    mkdir -p "$WORKDIR/logs"
    cd "$WORKDIR"
    
    START_TIME=$(date +%s)
    
    case $MODE in
        full)
            # Full run: simulation + analysis + visualization
            check_gromacs || exit 1
            for pdb in "${PDB_LIST[@]}"; do
                check_file "$pdb" "$(basename $pdb)" || exit 1
            done
            check_python_modules || exit 1
            check_gpu
            
            if [ "$USE_GPU" = true ]; then
                log "Running in GPU-accelerated mode (GPU $GPU_ID)"
                log "EM flags: $GPU_EM_FLAGS"
                log "MD flags: $GPU_MD_FLAGS"
            else
                log_warn "Running in CPU-only mode ($NTHREADS threads)"
            fi
            echo ""
            
            # Process all structures
            for i in "${!PDB_LIST[@]}"; do
                local struct_num=$((i+1))
                local struct_name="structure_${struct_num}"
                process_structure "${PDB_LIST[$i]}" "$struct_name" "$WORKDIR/$struct_name"
                calculate_interaction_energy "$WORKDIR/$struct_name" "$struct_name"
                generate_trajectory_outputs "$WORKDIR/$struct_name" "$struct_name"
                generate_structure_visualization "$WORKDIR/$struct_name" "$struct_name"
            done
            
            # Compare results
            compare_results
            ;;
            
        sim-only)
            # Simulation and analysis only (no plots/visualization)
            check_gromacs || exit 1
            for pdb in "${PDB_LIST[@]}"; do
                check_file "$pdb" "$(basename $pdb)" || exit 1
            done
            check_python_modules || exit 1
            check_gpu
            
            if [ "$USE_GPU" = true ]; then
                log "Running in GPU-accelerated mode (GPU $GPU_ID)"
            else
                log_warn "Running in CPU-only mode ($NTHREADS threads)"
            fi
            log "NOTE: Visualization skipped (use --mode viz-only to generate later)"
            echo ""
            
            # Process all structures
            for i in "${!PDB_LIST[@]}"; do
                local struct_num=$((i+1))
                local struct_name="structure_${struct_num}"
                process_structure "${PDB_LIST[$i]}" "$struct_name" "$WORKDIR/$struct_name"
                calculate_interaction_energy "$WORKDIR/$struct_name" "$struct_name"
                generate_trajectory_outputs "$WORKDIR/$struct_name" "$struct_name"
            done
            
            # Compare results
            compare_results
            ;;
            
        viz-only)
            # Visualization only (requires existing simulation data)
            check_python_modules || exit 1
            log "Generating visualization from existing data..."
            echo ""
            
            # Check if data exists for all structures
            for i in "${!PDB_LIST[@]}"; do
                local struct_name="structure_$((i+1))"
                if [[ ! -d "$WORKDIR/$struct_name" ]]; then
                    log_error "No simulation data found for $struct_name in $WORKDIR"
                    log_error "Run with --mode full or --mode sim-only first"
                    exit 1
                fi
            done
            
            # Generate visualization for all structures
            for i in "${!PDB_LIST[@]}"; do
                local struct_name="structure_$((i+1))"
                generate_structure_visualization "$WORKDIR/$struct_name" "$struct_name"
            done
            echo ""
            
            # Regenerate comparison report
            compare_results
            ;;
            
        compare-only)
            # Compare existing structures, skip completed steps
            check_python_modules || exit 1
            log "Running comparison on existing structures (skipping completed steps)..."
            echo ""
            
            # Show status of each structure
            log "Structure status:"
            local structures_found=0
            for i in "${!PDB_LIST[@]}"; do
                local struct_num=$((i+1))
                local struct_name="structure_${struct_num}"
                local struct_dir="$WORKDIR/$struct_name"
                if [[ -d "$struct_dir" ]]; then
                    local status=$(get_structure_status "$struct_dir")
                    log "  $struct_name: $status"
                    structures_found=$((structures_found + 1))
                else
                    log "  $struct_name: not found"
                fi
            done
            echo ""
            
            if [[ $structures_found -lt 2 ]]; then
                log_error "Need at least 2 structures for comparison. Found: $structures_found"
                exit 1
            fi
            
            # Try to complete any incomplete steps for existing structures
            for i in "${!PDB_LIST[@]}"; do
                local struct_num=$((i+1))
                local struct_name="structure_${struct_num}"
                local struct_dir="$WORKDIR/$struct_name"
                
                if [[ -d "$struct_dir" ]]; then
                    # Try IE calculation if not done
                    calculate_interaction_energy "$struct_dir" "$struct_name" || true
                    
                    # Try trajectory outputs if not done
                    generate_trajectory_outputs "$struct_dir" "$struct_name" || true
                    
                    # Try visualization if not done
                    generate_structure_visualization "$struct_dir" "$struct_name" || true
                fi
            done
            echo ""
            
            # Run comparison
            compare_results
            ;;
    esac
    
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    
    log_section "Comparison Complete"
    log "Mode: $MODE"
    log "Structures compared: ${#PDB_LIST[@]}"
    log "Total time: $((ELAPSED / 3600))h $((ELAPSED % 3600 / 60))m $((ELAPSED % 60))s"
    log "Results saved to: $WORKDIR"
    if [[ "$MODE" != "viz-only" ]]; then
        log "GPU Mode: $([ "$USE_GPU" = true ] && echo 'GPU-accelerated' || echo 'CPU-only')"
    fi
    log ""
    log "Key files per structure:"
    log "  - statistics/md_statistics.json"
    log "  - trajectories/md_center.xtc"
    log "  - structures/final_structure.pdb"
    if [[ "$MODE" != "sim-only" ]]; then
        log "  - visualization/visualize_trajectory.pml"
        log "  - plots/*.png"
    fi
    if [[ ${#PDB_LIST[@]} -gt 2 ]]; then
        log ""
        log "Multi-structure comparison plots:"
        log "  - plots/md_combined_radar.png"
        log "  - plots/md_combined_ranking.png"
        log "  - plots/md_combined_metrics.png"
        log "  - md_comparison_summary.txt"
    fi
    
    # Show comparison summary
    echo ""
    if [[ ${#PDB_LIST[@]} -gt 2 ]]; then
        cat "$WORKDIR/md_comparison_summary.txt" 2>/dev/null | head -60 || cat "$WORKDIR/md_comparison_report.txt" 2>/dev/null | head -60 || true
    else
        cat "$WORKDIR/md_comparison_report.txt" 2>/dev/null | head -60 || true
    fi
}

main "$@"
