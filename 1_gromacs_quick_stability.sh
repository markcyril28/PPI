#!/usr/bin/env bash
set -e

#######################################################################
# GROMACS Quick Stability Comparison Script (GPU)
#
# Quick comparison of multiple PPI structures using energy minimization.
# Uses Python modules for analysis and reporting.
#
# Usage: ./gromacs_quick_stability.sh [OUTPUT_DIR]
#
# Configure structures to compare by editing the PDB_LIST array below.
#######################################################################

#------------------------------------------------------------------------------
# CONFIGURATION
#------------------------------------------------------------------------------

# Override existing results (true) or skip if already processed (false)
OVERRIDE_EXISTING=true

# Maximum threads to use for GROMACS
MAX_THREADS=32

# Input PDB files - Comment/uncomment to select which structures to compare
# All uncommented files will be processed and compared against each other
INPUT_BASE="/mnt/local3.5tb/home/mcmercado/PPI/inputs"

# Structure list - uncomment the files you want to compare
# (minimum 2 structures required)
PDB_LIST=(
    # SmelGRF-GIF files:
    #"${INPUT_BASE}/SmelGRF-GIF/fold_1_x_1_model_0.pdb"
    #"${INPUT_BASE}/SmelGRF-GIF/fold_2_x_1_model_0.pdb"
    #"${INPUT_BASE}/SmelGRF-GIF/fold_recruitment_of_smelgif_to_swi_model_0.pdb"
    
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
generate_output_dir_name() {
    local names=()
    for pdb in "${PDB_LIST[@]}"; do
        # Extract base name without extension and path
        local basename=$(basename "$pdb" .pdb)
        # Shorten the name for readability (take first meaningful part)
        local shortname=$(echo "$basename" | sed 's/_model_0.*//; s/_2025_01_09.*//; s/fold_//')
        names+=("$shortname")
    done
    # Join names with "_"
    local joined=$(IFS="_"; echo "${names[*]}")
    echo "${OUTPUT_BASE}/1_quick_stability/${joined}"
}

# Set OUTPUT_DIR dynamically
OUTPUT_DIR=$(generate_output_dir_name)

# Validate we have at least 2 structures
if [ ${#PDB_LIST[@]} -lt 2 ]; then
    echo "ERROR: At least 2 structures must be uncommented in PDB_LIST"
    exit 1
fi

# GROMACS settings (override defaults)
BOX_DISTANCE=1.0
EM_STEPS=10000
NTHREADS=$MAX_THREADS

# Source common functions and configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/modules/gromacs_common.sh"

#------------------------------------------------------------------------------
# FUNCTIONS
#------------------------------------------------------------------------------

check_requirements() {
    check_gromacs || exit 1
    for pdb in "${PDB_LIST[@]}"; do
        check_file "$pdb" "$(basename $pdb)" || exit 1
    done
    check_python_modules || exit 1
}

process_structure() {
    local pdb="$1"
    local name="$2"
    local outdir="$3"
    
    log "Processing $name: $(basename $pdb)"
    
    # Skip if already processed (unless OVERRIDE_EXISTING is true)
    if [[ "$OVERRIDE_EXISTING" != "true" && -f "$outdir/metrics.txt" && -f "$outdir/em.gro" ]]; then
        log "  ✓ Already processed, skipping (found metrics.txt and em.gro)"
        log "    Set OVERRIDE_EXISTING=true to reprocess"
        return 0
    fi
    
    # If overriding, clean up old files
    if [[ "$OVERRIDE_EXISTING" == "true" && -d "$outdir" ]]; then
        log "  Overriding existing results..."
        rm -rf "$outdir"
    fi
    
    mkdir -p "$outdir/logs"
    cd "$outdir"
    
    # Clean PDB using Python CLI module
    python3 -m gromacs_utils.cli prepare-structure "$pdb" -o clean.pdb
    
    # Get chain info
    get_chain_info "$pdb" "$outdir"
    
    # Generate topology
    log "  Generating topology..."
    echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh > logs/pdb2gmx.log 2>&1
    
    # Create index
    create_chain_index clean.pdb protein.gro index.ndx
    
    # Setup box and solvate
    log "  Setting up simulation box..."
    $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron > logs/editconf.log 2>&1
    $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top > logs/solvate.log 2>&1
    
    # Generate EM MDP using Python CLI module
    python3 -m gromacs_utils.cli generate-mdp em -o em.mdp --em-steps $EM_STEPS --em-tolerance 10.0
    
    # Add ions
    $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 2>&1 > logs/grompp_ions.log
    echo "SOL" | $GMX_BIN genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral 2>&1 > logs/genion.log
    
    # Energy minimization
    log "  Running energy minimization (GPU)..."
    $GMX_BIN grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr 2>&1 > logs/grompp_em.log
    run_em em logs
    
    # Extract energies
    log "  Extracting metrics..."
    echo -e "Potential\nCoul-SR\nLJ-SR\nPressure\n\n" | $GMX_BIN energy -f em.edr -o energies.xvg 2>&1 > logs/energy.log || true
    
    # Interface analysis
    echo -e "ChainA\nChainB" | $GMX_BIN hbond -s em.tpr -f em.gro -n index.ndx -num hbonds.xvg 2>&1 > logs/hbond.log || true
    echo -e "ChainA\nChainB" | $GMX_BIN mindist -s em.tpr -f em.gro -n index.ndx -od mindist.xvg -on numcont.xvg -d 0.6 2>&1 > logs/mindist.log || true
    echo "Protein" | $GMX_BIN sasa -s em.tpr -f em.gro -o sasa.xvg 2>&1 > logs/sasa.log || true
    echo "ChainA" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o sasa_chainA.xvg 2>&1 > logs/sasa_a.log || true
    echo "ChainB" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o sasa_chainB.xvg 2>&1 > logs/sasa_b.log || true
    
    # Extract metrics using Python module
    extract_metrics "$outdir" metrics.txt
    
    # Generate visualization scripts
    generate_visualization "$outdir" interface
}

compare_results() {
    log "Comparing results..."
    
    local num_structs=${#PDB_LIST[@]}
    
    # Create a combined summary file
    local summary_file="$OUTPUT_DIR/structures_summary.txt"
    {
        echo "=============================================="
        echo "QUICK STABILITY COMPARISON SUMMARY"
        echo "=============================================="
        echo ""
        echo "Structures analyzed: $num_structs"
        echo "Generated: $(date)"
        echo ""
        
        for i in "${!PDB_LIST[@]}"; do
            local struct_num=$((i+1))
            local struct_name="structure_${struct_num}"
            local struct_dir="$OUTPUT_DIR/$struct_name"
            local pdb_name=$(basename "${PDB_LIST[$i]}" .pdb)
            
            echo "--- Structure $struct_num: $pdb_name ---"
            
            if [[ -f "$struct_dir/metrics.txt" ]]; then
                cat "$struct_dir/metrics.txt" | sed 's/^/  /'
            else
                echo "  (metrics unavailable)"
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
                local struct1_dir="$OUTPUT_DIR/$struct1_name"
                local struct2_dir="$OUTPUT_DIR/$struct2_name"
                
                # Check if both structures have data
                if [[ -d "$struct1_dir" && -d "$struct2_dir" ]]; then
                    log "  Comparing $struct1_name vs $struct2_name..."
                    
                    python3 -m gromacs_utils.results_comparator \
                        --workdir "$OUTPUT_DIR" \
                        --pdb1 "${PDB_LIST[$i]}" \
                        --pdb2 "${PDB_LIST[$j]}" \
                        --struct1-dir "$struct1_name" \
                        --struct2-dir "$struct2_name" \
                        --output "comparison_${struct1_num}_vs_${struct2_num}.txt" 2>&1 || true
                else
                    log_warn "  Skipping $struct1_name vs $struct2_name (missing data)"
                fi
            done
        done
        
        # Generate combined multi-structure comparison (radar plot, ranking) for 3+ structures
        if [ $num_structs -gt 2 ]; then
            log "Generating combined multi-structure comparison..."
            
            # Build PDB list arguments
            local pdb_args=""
            local struct_args=""
            for i in "${!PDB_LIST[@]}"; do
                pdb_args="$pdb_args ${PDB_LIST[$i]}"
                struct_args="$struct_args structure_$((i+1))"
            done
            
            python3 -m gromacs_utils.results_comparator \
                --workdir "$OUTPUT_DIR" \
                --multi \
                --pdb-list $pdb_args \
                --struct-dirs $struct_args 2>&1 || true
        fi
        
        # Also generate the default comparison report (structure_1 vs structure_2)
        if [[ -d "$OUTPUT_DIR/structure_1" && -d "$OUTPUT_DIR/structure_2" ]]; then
            python3 -m gromacs_utils.results_comparator \
                --workdir "$OUTPUT_DIR" \
                --pdb1 "${PDB_LIST[0]}" \
                --pdb2 "${PDB_LIST[1]}" \
                --struct1-dir "structure_1" \
                --struct2-dir "structure_2" 2>&1 || true
        fi
    fi
}

#------------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------------

main() {
    log_section "GROMACS Quick Stability Comparison (GPU)"
    
    check_requirements
    
    log "Structures to compare: ${#PDB_LIST[@]}"
    for i in "${!PDB_LIST[@]}"; do
        log "  Structure $((i+1)): $(basename ${PDB_LIST[$i]})"
    done
    log "Output: $OUTPUT_DIR"
    log "GROMACS: $($GMX_BIN --version 2>&1 | head -1)"
    echo ""
    
    mkdir -p "$OUTPUT_DIR"
    
    START=$(date +%s)
    
    # Process all structures
    for i in "${!PDB_LIST[@]}"; do
        local struct_num=$((i+1))
        local struct_name="structure_${struct_num}"
        process_structure "${PDB_LIST[$i]}" "$struct_name" "$OUTPUT_DIR/$struct_name"
        echo ""
    done
    
    # Compare results
    compare_results
    
    END=$(date +%s)
    
    log_section "Complete"
    log "Structures compared: ${#PDB_LIST[@]}"
    log "Total time: $((END - START)) seconds"
    log "Results: $OUTPUT_DIR"
    log ""
    log "Key files:"
    log "  - structures_summary.txt"
    log "  - comparison_report.txt"
    for i in "${!PDB_LIST[@]}"; do
        log "  - structure_$((i+1))/metrics.txt"
    done
    
    # Show summary
    echo ""
    cat "$OUTPUT_DIR/structures_summary.txt" 2>/dev/null || true
    echo ""
    cat "$OUTPUT_DIR/comparison_report.txt" 2>/dev/null || true
}

main "$@"
