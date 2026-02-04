#!/usr/bin/env python3
"""
MD Results Comparator Module for GROMACS Stability Analysis

Compares full MD simulation stability metrics between two structures
and generates comprehensive comparison reports with visualizations.
"""

import os
import json
import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional, Any, List


# Metric definitions: (display_name, unit, better_when, explanation)
METRIC_DESCRIPTIONS = {
    'rmsd_mean': ('RMSD (mean)', 'nm', 'lower',
        'Root Mean Square Deviation from initial structure. Lower values indicate '
        'the protein complex maintains its initial conformation better throughout simulation, '
        'suggesting structural stability.'),
    'rmsd_std': ('RMSD (std dev)', 'nm', 'lower',
        'Standard deviation of RMSD. Lower values indicate consistent stability with '
        'fewer structural fluctuations during simulation.'),
    'rmsd_final': ('RMSD (final)', 'nm', 'lower',
        'Final RMSD value at end of simulation. Lower values indicate the complex '
        'remained close to its initial structure.'),
    'rg_mean': ('Radius of Gyration (mean)', 'nm', 'similar',
        'Measure of overall protein compactness. Stable values indicate '
        'the complex maintains its overall shape during simulation.'),
    'rg_std': ('Radius of Gyration (std)', 'nm', 'lower',
        'Fluctuation in radius of gyration. Lower values indicate the complex '
        'does not expand/contract significantly during simulation.'),
    'hbonds_mean': ('H-bonds (mean)', 'count', 'higher',
        'Average number of hydrogen bonds between chains. More H-bonds indicate '
        'stronger, more specific protein-protein interactions at the interface.'),
    'hbonds_std': ('H-bonds (std)', 'count', 'lower',
        'Fluctuation in H-bond count. Lower values indicate stable, '
        'persistent hydrogen bonding throughout simulation.'),
    'mindist_mean': ('Min Distance (mean)', 'nm', 'lower',
        'Average minimum distance between chains. Lower values indicate '
        'tighter packing at the protein-protein interface.'),
    'mindist_std': ('Min Distance (std)', 'nm', 'lower',
        'Fluctuation in minimum distance. Lower values indicate consistent, '
        'stable interface contact throughout simulation.'),
    'sasa_mean': ('SASA (mean)', 'nm²', 'context',
        'Solvent Accessible Surface Area. Informational metric - the buried '
        'surface area (difference from unbound) is more relevant for binding.'),
    'sasa_std': ('SASA (std)', 'nm²', 'lower',
        'SASA fluctuation. Lower values indicate the exposed surface area '
        'remains consistent during simulation.'),
    'coul_mean': ('Coulomb Energy (mean)', 'kJ/mol', 'lower',
        'Electrostatic interaction energy between chains. More negative values '
        'indicate stronger favorable electrostatic interactions.'),
    'lj_mean': ('LJ Energy (mean)', 'kJ/mol', 'lower',
        'Lennard-Jones (van der Waals) interaction energy. More negative values '
        'indicate better shape complementarity at the interface.'),
    'total_ie_mean': ('Total Interaction Energy', 'kJ/mol', 'lower',
        'Sum of Coulomb and LJ energies. More negative values indicate '
        'stronger overall binding affinity between the protein chains.'),
    'potential': ('Potential Energy', 'kJ/mol', 'lower',
        'Total potential energy of the system. Lower values indicate a more '
        'thermodynamically stable complex conformation.'),
    'contacts': ('Interface Contacts', 'count', 'higher',
        'Number of atomic contacts at the interface. More contacts indicate '
        'a larger, more extensive binding interface.'),
    'min_distance': ('Minimum Distance', 'nm', 'lower',
        'Closest approach between chains. Lower values indicate tighter binding.'),
}


def read_metrics(filepath: str) -> Dict[str, float]:
    """Read stability metrics from a metrics file."""
    metrics = {}
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                if ':' in line and not line.startswith('='):
                    parts = line.strip().split(':')
                    if len(parts) == 2:
                        try:
                            metrics[parts[0].strip()] = float(parts[1].strip())
                        except ValueError:
                            pass
    return metrics


def read_json_metrics(filepath: str) -> Dict[str, Any]:
    """Read metrics from JSON file."""
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return json.load(f)
    return {}


def collect_all_metrics(workdir: str, struct_dir: str) -> Dict[str, float]:
    """Collect metrics from all available sources in a structure directory."""
    struct_path = os.path.join(workdir, struct_dir)
    metrics = {}
    
    # Try stability_metrics.txt first (MD analysis)
    stability_file = os.path.join(struct_path, "stability_metrics.txt")
    if os.path.exists(stability_file):
        metrics.update(read_metrics(stability_file))
    
    # Try metrics.json (quick stability)
    json_file = os.path.join(struct_path, "metrics.json")
    if os.path.exists(json_file):
        data = read_json_metrics(json_file)
        for k, v in data.items():
            if isinstance(v, (int, float)) and v is not None:
                metrics[k] = float(v)
    
    # Try statistics/md_statistics.json
    stats_file = os.path.join(struct_path, "statistics", "md_statistics.json")
    if os.path.exists(stats_file):
        data = read_json_metrics(stats_file)
        for k, v in data.items():
            if isinstance(v, (int, float)) and v is not None:
                metrics[k] = float(v)
            elif isinstance(v, dict):
                for sub_k, sub_v in v.items():
                    if isinstance(sub_v, (int, float)) and sub_v is not None:
                        metrics[f"{k}_{sub_k}"] = float(sub_v)
    
    return metrics


def generate_md_comparison_plots(
    workdir: str,
    m1: Dict[str, float],
    m2: Dict[str, float],
    name1: str,
    name2: str
) -> List[str]:
    """Generate comparison plots for MD simulation results."""
    generated_plots = []
    plots_dir = os.path.join(workdir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib not available, skipping plots")
        return generated_plots
    
    # Shorten names
    short1 = name1[:20] + "..." if len(name1) > 20 else name1
    short2 = name2[:20] + "..." if len(name2) > 20 else name2
    
    # Get common metrics
    all_keys = sorted(set(list(m1.keys()) + list(m2.keys())))
    common_keys = [k for k in all_keys 
                   if k in m1 and k in m2 
                   and m1[k] is not None and m2[k] is not None]
    
    if not common_keys:
        return generated_plots
    
    # 1. Multi-panel bar chart comparison
    n_metrics = len(common_keys)
    fig, axes = plt.subplots(min(n_metrics, 8), 1, figsize=(12, min(n_metrics, 8) * 2.5))
    if n_metrics == 1:
        axes = [axes]
    
    for ax, key in zip(axes, common_keys[:8]):
        info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', ''))
        display_name = info[0]
        unit = info[1]
        better = info[2]
        
        v1, v2 = m1[key], m2[key]
        
        # Determine colors
        if better == 'lower':
            colors = ['green' if v1 < v2 else 'coral', 
                      'green' if v2 < v1 else 'coral']
        elif better == 'higher':
            colors = ['green' if v1 > v2 else 'coral',
                      'green' if v2 > v1 else 'coral']
        else:
            colors = ['steelblue', 'steelblue']
        
        bars = ax.barh(['S1', 'S2'], [v1, v2], color=colors)
        ax.set_xlabel(f'{display_name} ({unit})' if unit else display_name)
        
        for bar, val in zip(bars, [v1, v2]):
            ax.text(bar.get_width() + 0.02 * max(abs(v1), abs(v2)), 
                    bar.get_y() + bar.get_height()/2,
                    f'{val:.3f}', va='center', fontsize=9)
    
    plt.tight_layout()
    metrics_plot = os.path.join(plots_dir, "md_comparison_metrics.png")
    plt.savefig(metrics_plot, dpi=150, bbox_inches='tight')
    plt.close()
    generated_plots.append(metrics_plot)
    
    # 2. Score summary plot
    scores = {'Structure 1': 0, 'Structure 2': 0}
    for key in common_keys:
        info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', ''))
        better = info[2]
        v1, v2 = m1[key], m2[key]
        
        if better == 'lower' and v1 < v2:
            scores['Structure 1'] += 1
        elif better == 'lower' and v2 < v1:
            scores['Structure 2'] += 1
        elif better == 'higher' and v1 > v2:
            scores['Structure 1'] += 1
        elif better == 'higher' and v2 > v1:
            scores['Structure 2'] += 1
    
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['green' if scores['Structure 1'] > scores['Structure 2'] else 'coral',
              'green' if scores['Structure 2'] > scores['Structure 1'] else 'coral']
    
    bars = ax.bar([f'Structure 1\n{short1}', f'Structure 2\n{short2}'],
                  [scores['Structure 1'], scores['Structure 2']],
                  color=colors, edgecolor='black', linewidth=2)
    
    ax.set_ylabel('Favorable Metrics Count', fontsize=12)
    ax.set_title('MD Simulation Stability Score Comparison', fontsize=14)
    
    for bar, score in zip(bars, [scores['Structure 1'], scores['Structure 2']]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                str(score), ha='center', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    summary_plot = os.path.join(plots_dir, "md_comparison_summary.png")
    plt.savefig(summary_plot, dpi=150)
    plt.close()
    generated_plots.append(summary_plot)
    
    return generated_plots


def compare_md_results(
    workdir: str,
    pdb1_path: str,
    pdb2_path: str,
    struct1_dir: str = "structure_1",
    struct2_dir: str = "structure_2",
    output_file: str = "md_comparison_report.txt"
) -> Tuple[Dict[str, int], Optional[str]]:
    """
    Compare MD simulation stability results between two structures.
    
    Returns:
        Tuple of (scores dict, winner string or None)
    """
    # Collect all metrics from both structures
    m1 = collect_all_metrics(workdir, struct1_dir)
    m2 = collect_all_metrics(workdir, struct2_dir)
    
    pdb1_name = os.path.basename(pdb1_path)
    pdb2_name = os.path.basename(pdb2_path)
    
    # Build detailed report
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append(" MD SIMULATION STABILITY COMPARISON")
    report_lines.append("=" * 80)
    report_lines.append("")
    report_lines.append(f"Structure 1: {pdb1_name}")
    report_lines.append(f"  Path: {pdb1_path}")
    report_lines.append(f"Structure 2: {pdb2_name}")
    report_lines.append(f"  Path: {pdb2_path}")
    report_lines.append("")
    report_lines.append("-" * 80)
    report_lines.append(" METRIC-BY-METRIC ANALYSIS")
    report_lines.append("-" * 80)
    
    scores = {"s1": 0, "s2": 0}
    detailed_comparisons = []
    
    # Get all metrics
    all_keys = sorted(set(list(m1.keys()) + list(m2.keys())))
    
    for key in all_keys:
        v1 = m1.get(key)
        v2 = m2.get(key)
        
        if v1 is None and v2 is None:
            continue
        
        info = METRIC_DESCRIPTIONS.get(key, (key, '', 'context', f'Metric: {key}'))
        display_name = info[0]
        unit = info[1]
        better = info[2]
        explanation = info[3] if len(info) > 3 else ''
        
        report_lines.append("")
        report_lines.append(f"▸ {display_name}")
        report_lines.append(f"  {'─' * 70}")
        
        v1_str = f"{v1:.4f} {unit}" if v1 is not None else "N/A"
        v2_str = f"{v2:.4f} {unit}" if v2 is not None else "N/A"
        
        report_lines.append(f"  Structure 1: {v1_str}")
        report_lines.append(f"  Structure 2: {v2_str}")
        
        winner = ""
        winner_explanation = ""
        
        if v1 is not None and v2 is not None:
            diff = abs(v1 - v2)
            pct_diff = (diff / max(abs(v1), abs(v2), 0.0001)) * 100
            
            if better == 'lower':
                if v1 < v2:
                    winner = "Structure 1"
                    scores['s1'] += 1
                    winner_explanation = f"Structure 1 is {pct_diff:.1f}% lower (more favorable)"
                elif v2 < v1:
                    winner = "Structure 2"
                    scores['s2'] += 1
                    winner_explanation = f"Structure 2 is {pct_diff:.1f}% lower (more favorable)"
            elif better == 'higher':
                if v1 > v2:
                    winner = "Structure 1"
                    scores['s1'] += 1
                    winner_explanation = f"Structure 1 is {pct_diff:.1f}% higher (more favorable)"
                elif v2 > v1:
                    winner = "Structure 2"
                    scores['s2'] += 1
                    winner_explanation = f"Structure 2 is {pct_diff:.1f}% higher (more favorable)"
            elif better == 'similar':
                winner_explanation = f"Similar values preferred (diff: {pct_diff:.1f}%)"
            else:
                winner_explanation = "Informational metric"
        
        if winner:
            report_lines.append(f"  ★ Winner: {winner} — {winner_explanation}")
            detailed_comparisons.append((display_name, winner, winner_explanation))
        else:
            report_lines.append(f"  ○ {winner_explanation}")
        
        if explanation:
            report_lines.append(f"  ")
            report_lines.append(f"  Interpretation: {explanation}")
    
    # Conclusion
    report_lines.append("")
    report_lines.append("=" * 80)
    report_lines.append(" OVERALL CONCLUSION")
    report_lines.append("=" * 80)
    report_lines.append("")
    report_lines.append(f"  Favorable metrics for Structure 1: {scores['s1']}")
    report_lines.append(f"  Favorable metrics for Structure 2: {scores['s2']}")
    report_lines.append("")
    
    winner_str = None
    if scores['s1'] > scores['s2']:
        winner_str = "Structure 1"
        report_lines.append(f"  ╔══════════════════════════════════════════════════════════════════════╗")
        report_lines.append(f"  ║  WINNER: Structure 1 ({pdb1_name[:50]})  ")
        report_lines.append(f"  ╚══════════════════════════════════════════════════════════════════════╝")
        report_lines.append("")
        report_lines.append(f"  Structure 1 demonstrates MORE STABLE protein-protein interaction")
        report_lines.append(f"  based on {scores['s1']} out of {scores['s1'] + scores['s2']} comparable metrics.")
    elif scores['s2'] > scores['s1']:
        winner_str = "Structure 2"
        report_lines.append(f"  ╔══════════════════════════════════════════════════════════════════════╗")
        report_lines.append(f"  ║  WINNER: Structure 2 ({pdb2_name[:50]})  ")
        report_lines.append(f"  ╚══════════════════════════════════════════════════════════════════════╝")
        report_lines.append("")
        report_lines.append(f"  Structure 2 demonstrates MORE STABLE protein-protein interaction")
        report_lines.append(f"  based on {scores['s2']} out of {scores['s1'] + scores['s2']} comparable metrics.")
    else:
        winner_str = "Tie"
        report_lines.append(f"  ╔══════════════════════════════════════════════════════════════════════╗")
        report_lines.append(f"  ║  RESULT: Both structures show SIMILAR stability                      ║")
        report_lines.append(f"  ╚══════════════════════════════════════════════════════════════════════╝")
    
    # Key factors
    if detailed_comparisons:
        report_lines.append("")
        report_lines.append("  Key factors:")
        for metric_name, winner, explanation in detailed_comparisons[:5]:
            report_lines.append(f"    • {metric_name}: {winner} — {explanation}")
    
    # Methodology
    report_lines.append("")
    report_lines.append("-" * 80)
    report_lines.append(" METHODOLOGY & INTERPRETATION")
    report_lines.append("-" * 80)
    report_lines.append("""
  This comparison evaluates protein-protein interaction stability based on
  molecular dynamics simulation results. Key stability indicators include:

  1. RMSD (Root Mean Square Deviation)
     Measures structural deviation from the starting structure.
     Lower RMSD indicates the complex maintains its fold during simulation.

  2. HYDROGEN BONDS
     More H-bonds at the interface indicate stronger, specific interactions.
     Stable H-bond count throughout simulation suggests persistent binding.

  3. INTERACTION ENERGY (Coulomb + LJ)
     More negative values indicate stronger binding between chains.
     This is a key indicator of binding affinity.

  4. MINIMUM DISTANCE
     Lower values indicate tighter packing at the interface.
     Stable minimum distance suggests consistent binding.

  5. RADIUS OF GYRATION
     Stable values indicate the complex maintains its shape.
     Large fluctuations may indicate structural instability.

  NOTE: For publication-quality binding affinity predictions, consider:
  - Longer MD simulations (100+ ns)
  - MM-PBSA or MM-GBSA free energy calculations
  - Multiple independent replicas
  - Experimental validation
""")
    report_lines.append("=" * 80)
    
    # Print and save
    report_text = "\n".join(report_lines)
    print(report_text)
    
    report_file = os.path.join(workdir, output_file)
    with open(report_file, 'w') as f:
        f.write(report_text)
    print(f"\n  Report saved to: {report_file}")
    
    # Generate plots
    print("\n  Generating comparison plots...")
    plots = generate_md_comparison_plots(workdir, m1, m2, pdb1_name, pdb2_name)
    if plots:
        print(f"  Generated {len(plots)} plots in {os.path.join(workdir, 'plots')}/")
    
    return scores, winner_str


def main():
    """Command-line interface for MD results comparison."""
    parser = argparse.ArgumentParser(
        description='Compare MD simulation stability results between two structures'
    )
    parser.add_argument('--workdir', required=True,
                       help='Base working directory')
    parser.add_argument('--pdb1', required=True,
                       help='Path to first PDB file')
    parser.add_argument('--pdb2', required=True,
                       help='Path to second PDB file')
    parser.add_argument('--struct1-dir', default='structure_1',
                       help='Subdirectory for structure 1')
    parser.add_argument('--struct2-dir', default='structure_2',
                       help='Subdirectory for structure 2')
    parser.add_argument('--output', '-o', default='md_comparison_report.txt',
                       help='Output report filename')
    
    args = parser.parse_args()
    
    compare_md_results(
        args.workdir,
        args.pdb1,
        args.pdb2,
        args.struct1_dir,
        args.struct2_dir,
        args.output
    )


if __name__ == '__main__':
    main()
