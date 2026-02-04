#!/usr/bin/env bash
set -e

#######################################################################
# GROMACS Interface Analysis Script
#
# Performs detailed PPI interface analysis using Python modules.
# This is a minimal bash wrapper that calls the gromacs_utils modules.
#
# Usage: ./gromacs_interface_analysis.sh <PDB_FILE> [OUTPUT_DIR]
#######################################################################

#------------------------------------------------------------------------------
# CONFIGURATION - Modify these as needed
#------------------------------------------------------------------------------

# Input directory containing PDB files
INPUTS_DIR="${INPUTS_DIR:-/mnt/local3.5tb/home/mcmercado/PPI/inputs}"

# Override existing results (true) or skip if already processed (false)
OVERRIDE_EXISTING=false

# Output directory
OUTPUT_DIR="${2:-/mnt/local3.5tb/home/mcmercado/PPI/results_gromacs/3_interface_analysis}"

# GROMACS settings (override defaults)
BOX_DISTANCE=1.0
EM_STEPS=5000

# Source common functions and configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/modules/gromacs_common.sh"

#------------------------------------------------------------------------------
# PDB FILE LIST - Comment in/out as needed
#------------------------------------------------------------------------------

# Input directory base path
INPUTS_DIR="/mnt/local3.5tb/home/mcmercado/PPI/inputs"

# List of PDB files to process 
PDB_FILES=(
    # SmelDMP structures
    "${INPUTS_DIR}/SmelDMP/SmelDMP01.730_SmelHAP2.pdb"
    "${INPUTS_DIR}/SmelDMP/SmelDMP01.990_SmelHAP2.pdb"
    "${INPUTS_DIR}/SmelDMP/SmelDMP02_SmelHAP2.pdb"
    #"${INPUTS_DIR}/SmelDMP/SmelDMP04_SmelHAP2.pdb"
    "${INPUTS_DIR}/SmelDMP/SmelDMP10.200_SmelHAP2.pdb"
    #"${INPUTS_DIR}/SmelDMP/SmelDMP10_550_560_SmelHAP2.pdb"
    #"${INPUTS_DIR}/SmelDMP/SmelDMP12_SmelHAP2.pdb"
)

#------------------------------------------------------------------------------
# HELPER FUNCTIONS
#------------------------------------------------------------------------------

check_requirements() {
    check_gromacs || exit 1
    check_file "$PDB_INPUT" "Input file" || exit 1
    check_python_modules || exit 1
}

#------------------------------------------------------------------------------
# GROMACS FUNCTIONS
#------------------------------------------------------------------------------

prepare_structure() {
    log "Preparing structure..."
    
    # Use Python CLI module for structure conversion and cleaning
    python3 -m gromacs_utils.cli prepare-structure "$PDB_INPUT" -o clean.pdb
    
    # Generate topology
    log "  Generating topology..."
    echo "1" | $GMX_BIN pdb2gmx -f clean.pdb -o protein.gro -water $WATERMODEL -ff $FORCEFIELD -ignh > logs/pdb2gmx.log 2>&1
    
    # Create chain index using Python module
    create_chain_index clean.pdb protein.gro index.ndx
    
    # Solvate and minimize
    log "  Solvating..."
    $GMX_BIN editconf -f protein.gro -o boxed.gro -c -d $BOX_DISTANCE -bt dodecahedron > logs/editconf.log 2>&1
    $GMX_BIN solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top > logs/solvate.log 2>&1
    
    # Generate EM MDP using Python CLI module
    python3 -m gromacs_utils.cli generate-mdp em -o em.mdp --em-steps $EM_STEPS --em-tolerance 100.0
    
    log "  Running energy minimization..."
    $GMX_BIN grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 2 > logs/grompp.log 2>&1
    run_em em logs
}

run_gromacs_analysis() {
    log "Running GROMACS analysis tools..."
    
    # Inter-chain distances and contacts
    echo -e "ChainA\nChainB" | $GMX_BIN mindist -s em.tpr -f em.gro -n index.ndx \
        -od analysis/mindist.xvg -on analysis/numcont.xvg -d 0.6 2>&1 > logs/mindist.log || true
    
    # Hydrogen bonds
    echo -e "ChainA\nChainB" | $GMX_BIN hbond -s em.tpr -f em.gro -n index.ndx \
        -num analysis/hbonds.xvg 2>&1 > logs/hbond.log || true
    
    # SASA
    echo "Protein" | $GMX_BIN sasa -s em.tpr -f em.gro -o analysis/sasa_complex.xvg 2>&1 > logs/sasa_complex.log || true
    echo "ChainA" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o analysis/sasa_chainA.xvg 2>&1 > logs/sasa_A.log || true
    echo "ChainB" | $GMX_BIN sasa -s em.tpr -f em.gro -n index.ndx -o analysis/sasa_chainB.xvg 2>&1 > logs/sasa_B.log || true
    
    # Energy
    echo -e "Potential\nCoul-SR\nLJ-SR\n\n" | $GMX_BIN energy -f em.edr -o analysis/energies.xvg 2>&1 > logs/energy.log || true
}

#------------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------------

process_single_structure() {
    local PDB_INPUT="$1"
    
    STRUCT_NAME=$(basename "$PDB_INPUT" | sed 's/\.[^.]*$//')
    WORKDIR="${OUTPUT_DIR}/${STRUCT_NAME}"
    
    log "Structure: $STRUCT_NAME"
    log "Output: $WORKDIR"
    
    # Skip if already processed (unless OVERRIDE_EXISTING is true)
    if [[ "$OVERRIDE_EXISTING" != "true" && -f "$WORKDIR/interface_metrics.json" ]]; then
        log "  ✓ Already processed, skipping..."
        log "    Set OVERRIDE_EXISTING=true to reprocess"
        return 0
    fi
    
    # If overriding, clean up old files
    if [[ "$OVERRIDE_EXISTING" == "true" && -d "$WORKDIR" ]]; then
        log "  Overriding existing results..."
        rm -rf "$WORKDIR"
    fi
    
    # Create output structure and change to it
    # Use pushd/popd to ensure we return to original directory
    mkdir -p "$WORKDIR"
    pushd "$WORKDIR" > /dev/null
    setup_output_structure "$WORKDIR"
    
    START=$(date +%s)
    
    # Prepare and minimize structure
    prepare_structure
    
    # Run GROMACS analysis
    run_gromacs_analysis
    
    # Generate contact map using Python module
    log "Generating contact map..."
    python3 -m gromacs_utils.contact_map --workdir "$WORKDIR" --gro em.gro --pdb clean.pdb
    
    # Extract interface metrics and generate report
    log "Generating interface report..."
    python3 -m gromacs_utils.interface_analyzer --workdir "$WORKDIR" --name "$STRUCT_NAME" --output "$WORKDIR"
    
    # Generate visualization scripts
    log "Generating visualization scripts..."
    generate_visualization "$WORKDIR" interface
    
    # Convert structure to PDB for visualization
    echo "Protein" | $GMX_BIN trjconv -s em.tpr -f em.gro -o structures/structure_minimized.pdb > logs/trjconv.log 2>&1 || true
    
    # Generate plots
    log "Generating plots..."
    run_gnuplot_scripts plots
    
    # Generate matplotlib plots using Python CLI module
    python3 -m gromacs_utils.cli generate-plots contact -i analysis/contact_map.txt -o plots/contact_map_heatmap.png 2>&1 || true
    
    # Generate focused contact map with only critical residues
    python3 -c "
from gromacs_utils.plotting import plot_focused_contact_map
plot_focused_contact_map(
    'analysis/interface_residues.txt',
    'plots/contact_map_focused.png',
    max_pairs=50
)
" 2>&1 || true
    
    # Return to original directory
    popd > /dev/null
    
    END=$(date +%s)
    log "  ✓ Completed in $((END - START)) seconds"
}

main() {
    log_section "GROMACS Interface Analysis (Batch Mode)"
    
    # Check requirements once
    check_gromacs || exit 1
    check_python_modules || exit 1
    
    if [[ ${#PDB_FILES[@]} -eq 0 ]]; then
        log_error "No PDB files defined in PDB_FILES array"
        exit 1
    fi
    
    log "Found ${#PDB_FILES[@]} PDB files to process:"
    for pdb in "${PDB_FILES[@]}"; do
        log "  - $(basename "$pdb")"
    done
    echo ""
    
    TOTAL_START=$(date +%s)
    PROCESSED=0
    FAILED=0
    
    # Process each PDB file
    for PDB_INPUT in "${PDB_FILES[@]}"; do
        log_section "Processing: $(basename "$PDB_INPUT")"
        
        if check_file "$PDB_INPUT" "Input file"; then
            if process_single_structure "$PDB_INPUT"; then
                PROCESSED=$((PROCESSED + 1))
            else
                log_warn "Failed to process: $(basename "$PDB_INPUT")"
                FAILED=$((FAILED + 1))
            fi
        else
            FAILED=$((FAILED + 1))
        fi
        echo ""
    done
    
    TOTAL_END=$(date +%s)
    TOTAL_TIME=$((TOTAL_END - TOTAL_START))
    
    log_section "Batch Analysis Complete"
    log "Total structures: ${#PDB_FILES[@]}"
    log "Processed: $PROCESSED"
    log "Failed: $FAILED"
    log "Total time: $((TOTAL_TIME / 3600))h $((TOTAL_TIME % 3600 / 60))m $((TOTAL_TIME % 60))s"
    log "Output directory: $OUTPUT_DIR"
    log ""
    log "Key files per structure:"
    log "  - analysis/interface_analysis.txt (report)"
    log "  - interface_metrics.json (data)"
    log "  - visualization/visualize_interface.pml (PyMOL)"
    log "  - plots/contact_map_focused.png (contact heatmap)"
}

main "$@"
