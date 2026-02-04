#!/usr/bin/env python3
"""
Automatically extract critical residues from all MutateX results
Processes all results in run_results/ directory
"""

import os
import sys
from pathlib import Path
import subprocess

def find_pdb_file(pdb_name, inputs_dir):
    """Find the corresponding PDB file in inputs directory."""
    pdb_path = Path(inputs_dir) / f"{pdb_name}.pdb"
    if pdb_path.exists():
        return pdb_path
    
    # Try with _model0_checked suffix
    pdb_path_checked = Path(inputs_dir) / f"{pdb_name}_model0_checked.pdb"
    if pdb_path_checked.exists():
        return pdb_path_checked
    
    return None


def process_all_results(run_results_dir, inputs_dir, output_subdir, thresholds):
    """Process all MutateX results in run_results directory."""
    run_results_path = Path(run_results_dir)
    inputs_path = Path(inputs_dir)
    
    if not run_results_path.exists():
        print(f"[ERROR] Results directory not found: {run_results_dir}")
        return
    
    # Find all result directories
    result_dirs = [d for d in run_results_path.iterdir() if d.is_dir()]
    
    if not result_dirs:
        print(f"[WARN] No result directories found in {run_results_dir}")
        return
    
    print(f"[INFO] Found {len(result_dirs)} result directories")
    print(f"[INFO] Thresholds: destabilizing >{thresholds['destabilizing']}, "
          f"stabilizing <{thresholds['stabilizing']}, top N={thresholds['top_n']}")
    print("="*80)
    
    success_count = 0
    failed_count = 0
    
    for result_dir in result_dirs:
        pdb_name = result_dir.name
        print(f"\n[PROCESSING] {pdb_name}")
        
        # Check if results exist
        results_subdir = result_dir / "results"
        if not results_subdir.exists():
            print(f"  [SKIP] No results directory found in {result_dir}")
            failed_count += 1
            continue
        
        # Find corresponding PDB file
        pdb_file = find_pdb_file(pdb_name, inputs_path)
        if not pdb_file:
            print(f"  [SKIP] PDB file not found for {pdb_name}")
            failed_count += 1
            continue
        
        # Create output directory within the result directory
        output_dir = result_dir / output_subdir
        
        # Run extraction script
        cmd = [
            "python3",
            str(Path(__file__).parent / "extract_critical_residues.py"),
            str(result_dir),
            str(pdb_file),
            "-o", str(output_dir),
            "--destabilizing", str(thresholds['destabilizing']),
            "--stabilizing", str(thresholds['stabilizing']),
            "--top-n", str(thresholds['top_n'])
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(result.stdout)
            success_count += 1
        except subprocess.CalledProcessError as e:
            print(f"  [ERROR] Failed to process {pdb_name}")
            print(f"  {e.stderr}")
            failed_count += 1
        except Exception as e:
            print(f"  [ERROR] Unexpected error for {pdb_name}: {e}")
            failed_count += 1
    
    print("\n" + "="*80)
    print(f"[SUMMARY] Processed {len(result_dirs)} directories")
    print(f"  Success: {success_count}")
    print(f"  Failed:  {failed_count}")
    print(f"[OUTPUT] Critical residue files saved in each run_results subdirectory under '{output_subdir}/'")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Automatically extract critical residues from all MutateX results"
    )
    parser.add_argument(
        "--run-results",
        default="run_results",
        help="Directory containing MutateX result folders (default: run_results)"
    )
    parser.add_argument(
        "--inputs",
        default="inputs",
        help="Directory containing PDB files (default: inputs)"
    )
    parser.add_argument(
        "-o", "--output",
        default="critical_residues",
        help="Output subdirectory name within each result folder (default: critical_residues)"
    )
    parser.add_argument(
        "--destabilizing",
        type=float,
        default=2.0,
        help="DDG threshold for destabilizing mutations (default: 2.0 kcal/mol)"
    )
    parser.add_argument(
        "--stabilizing",
        type=float,
        default=-2.0,
        help="DDG threshold for stabilizing mutations (default: -2.0 kcal/mol)"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top residues to extract per category (default: 20)"
    )
    
    args = parser.parse_args()
    
    thresholds = {
        'destabilizing': args.destabilizing,
        'stabilizing': args.stabilizing,
        'top_n': args.top_n
    }
    
    process_all_results(args.run_results, args.inputs, args.output, thresholds)


if __name__ == "__main__":
    main()
