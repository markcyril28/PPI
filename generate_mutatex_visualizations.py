#!/usr/bin/env python3
"""
MutateX Visualization Generator
Generates comprehensive heatmaps and analysis plots from MutateX results.

Supports both:
1. Selfmutation data format (.dat files with avg/std/min/max columns)
2. Full mutation scan data (individual mutation files per residue)

Usage:
    python generate_mutatex_visualizations.py                    # Process all results in run_results/
    python generate_mutatex_visualizations.py <results_dir>      # Process single result directory
    python generate_mutatex_visualizations.py <results_dir> <pdb_name>  # Legacy mode
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Standard amino acid order (MutateX mutation_list.txt order)
# This matches the default mutation_list.txt template used by MutateX
AA_ORDER = ['G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S', 'T', 'C', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H']

# Color palettes matching MutateX examples (plasma for heatmaps, viridis for distributions)
HEATMAP_CMAP = 'plasma'
DIVERGING_CMAP = 'RdBu_r'  # For diverging data (positive/negative)
DISTRIB_CMAP = 'viridis'

# Default run_results directory (relative to script location)
DEFAULT_RUN_RESULTS_DIR = "run_results"


def parse_selfmutation_dat_file(dat_file_path, data_type='interface'):
    """Parse a selfmutation .dat file (format: residue_id avg std min max).
    
    Args:
        dat_file_path: Path to the .dat file
        data_type: 'interface' or 'folding'
    
    Returns:
        DataFrame with parsed data
    """
    data = []
    
    try:
        with open(dat_file_path, 'r') as f:
            lines = f.readlines()
        
        # Extract chain pair from filename if interface data
        chain_pair = None
        if 'selfmutation_energies_' in os.path.basename(dat_file_path):
            # Format: selfmutation_energies_A-B.dat
            chain_pair = os.path.basename(dat_file_path).replace('selfmutation_energies_', '').replace('.dat', '')
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) < 5:
                continue
            
            residue_id = parts[0]
            avg_ddg = float(parts[1])
            std_ddg = float(parts[2])
            min_ddg = float(parts[3])
            max_ddg = float(parts[4])
            
            # Parse residue_id: e.g., "MA1" -> M (aa), A (chain), 1 (position)
            # or "EA10" -> E (aa), A (chain), 10 (position)
            if len(residue_id) < 3:
                continue
            
            original_aa = residue_id[0]
            chain_id = residue_id[1]
            position_str = residue_id[2:]
            
            try:
                position = int(position_str)
            except ValueError:
                continue
            
            record = {
                "original_aa": original_aa,
                "position": position,
                "position_label": f"{original_aa}{chain_id}{position}",
                "chain": chain_id,
                "ddg": avg_ddg,
                "ddg_std": std_ddg,
                "ddg_min": min_ddg,
                "ddg_max": max_ddg,
                "data_type": data_type
            }
            
            if chain_pair:
                record["chain_pair"] = chain_pair
            
            data.append(record)
    
    except (IOError, ValueError) as e:
        print(f"[WARNING] Error parsing {dat_file_path}: {e}")
    
    return pd.DataFrame(data)


def find_selfmutation_data(results_dir):
    """Find selfmutation .dat files in the results directory.
    
    Looks for the standard MutateX output structure:
    - results_selfmutation/interface_ddgs/<model>/selfmutation_energies*.dat
    - results_selfmutation/mutation_ddgs/<model>/selfmutation_energies.dat
    
    Returns:
        tuple: (interface_dat_files, folding_dat_files)
    """
    interface_files = []
    folding_files = []
    
    # Check for selfmutation results structure
    selfmut_dir = os.path.join(results_dir, "results_selfmutation")
    
    if os.path.isdir(selfmut_dir):
        # Interface DDGs
        interface_ddg_dir = os.path.join(selfmut_dir, "interface_ddgs")
        if os.path.isdir(interface_ddg_dir):
            for model_dir in os.listdir(interface_ddg_dir):
                model_path = os.path.join(interface_ddg_dir, model_dir)
                if os.path.isdir(model_path):
                    for f in os.listdir(model_path):
                        if f.endswith('.dat'):
                            interface_files.append(os.path.join(model_path, f))
        
        # Folding/Mutation DDGs
        mutation_ddg_dir = os.path.join(selfmut_dir, "mutation_ddgs")
        if os.path.isdir(mutation_ddg_dir):
            for model_dir in os.listdir(mutation_ddg_dir):
                model_path = os.path.join(mutation_ddg_dir, model_dir)
                if os.path.isdir(model_path):
                    for f in os.listdir(model_path):
                        if f.endswith('.dat'):
                            folding_files.append(os.path.join(model_path, f))
    
    return interface_files, folding_files


def parse_selfmutation_results(results_dir):
    """Parse all selfmutation data from a results directory.
    
    Returns:
        tuple: (interface_df, folding_df)
    """
    interface_files, folding_files = find_selfmutation_data(results_dir)
    
    interface_dfs = []
    for f in interface_files:
        df = parse_selfmutation_dat_file(f, 'interface')
        if not df.empty:
            interface_dfs.append(df)
    
    folding_dfs = []
    for f in folding_files:
        df = parse_selfmutation_dat_file(f, 'folding')
        if not df.empty:
            folding_dfs.append(df)
    
    interface_df = pd.concat(interface_dfs, ignore_index=True) if interface_dfs else pd.DataFrame()
    folding_df = pd.concat(folding_dfs, ignore_index=True) if folding_dfs else pd.DataFrame()
    
    return interface_df, folding_df


def parse_interface_data(interface_dir):
    """Parse MutateX interface results (full mutation scan format)."""
    interface_data = []
    
    for chain_pair in os.listdir(interface_dir):
        chain_dir = os.path.join(interface_dir, chain_pair)
        if not os.path.isdir(chain_dir):
            continue
        
        for mutation_file in os.listdir(chain_dir):
            file_path = os.path.join(chain_dir, mutation_file)
            if not os.path.isfile(file_path):
                continue
            
            try:
                with open(file_path, 'r') as f:
                    lines = f.read().strip().split('\n')
                
                if len(lines) < 2:
                    continue
                
                position = ''.join(filter(str.isdigit, mutation_file))
                letters = mutation_file.replace(position, "")
                
                if len(letters) >= 2:
                    original_aa = letters[0]
                    chain_id = letters[1]
                else:
                    continue
                
                if not position:
                    continue
                
                position_label = f"{original_aa}{chain_id}{position}"
                
                # Parse all 20 mutations
                for i, line in enumerate(lines[1:21]):
                    values = line.split()
                    if values:
                        ddg_value = float(values[0])
                        target_aa = AA_ORDER[i]
                        
                        interface_data.append({
                            "original_aa": original_aa,
                            "target_aa": target_aa,
                            "position": int(position),
                            "position_label": position_label,
                            "chain": chain_id,
                            "ddg": ddg_value,
                            "chain_pair": chain_pair,
                            "mutation_label": f"{original_aa}{position}{target_aa}"
                        })
            except (ValueError, IOError) as e:
                continue
    
    return pd.DataFrame(interface_data)


def parse_folding_data(folding_dir):
    """Parse MutateX folding stability results."""
    folding_data = []
    
    for mutation_file in os.listdir(folding_dir):
        file_path = os.path.join(folding_dir, mutation_file)
        if not os.path.isfile(file_path):
            continue
        
        try:
            with open(file_path, 'r') as f:
                lines = f.read().strip().split('\n')
            
            if len(lines) < 2:
                continue
            
            position = ''.join(filter(str.isdigit, mutation_file))
            letters = mutation_file.replace(position, "")
            
            if len(letters) >= 2:
                original_aa = letters[0]
                chain_id = letters[1]
            else:
                continue
            
            if not position:
                continue
            
            position_label = f"{original_aa}{chain_id}{position}"
            
            for i, line in enumerate(lines[1:21]):
                values = line.split()
                if values:
                    ddg_value = float(values[0])
                    target_aa = AA_ORDER[i]
                    
                    folding_data.append({
                        "original_aa": original_aa,
                        "target_aa": target_aa,
                        "position": int(position),
                        "position_label": position_label,
                        "chain": chain_id,
                        "ddg": ddg_value,
                        "mutation_label": f"{original_aa}{position}{target_aa}"
                    })
        except (ValueError, IOError) as e:
            continue
    
    return pd.DataFrame(folding_data)


def generate_interface_heatmaps(interface_df, output_dir, pdbname):
    """Generate interface ΔΔG heatmaps per chain with square cells."""
    for chain in interface_df['chain'].unique():
        chain_data = interface_df[interface_df['chain'] == chain].copy()
        chain_data = chain_data.sort_values('position')
        chain_data['pos_label'] = chain_data['original_aa'] + chain_data['position'].astype(str)
        
        pivot = chain_data.pivot_table(index="target_aa", columns="pos_label", values="ddg", aggfunc='mean')
        pivot = pivot.reindex(AA_ORDER)
        
        # Split into chunks if too many columns
        max_cols = 100
        n_chunks = (len(pivot.columns) + max_cols - 1) // max_cols
        
        for chunk_idx in range(n_chunks):
            start_col = chunk_idx * max_cols
            end_col = min((chunk_idx + 1) * max_cols, len(pivot.columns))
            pivot_chunk = pivot.iloc[:, start_col:end_col]
            
            # Calculate figure size for square cells (cell_size in inches)
            cell_size = 0.4
            n_rows, n_cols = pivot_chunk.shape
            fig_width = n_cols * cell_size + 3  # Extra space for labels and colorbar
            fig_height = n_rows * cell_size + 2  # Extra space for title and labels
            
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            sns.heatmap(pivot_chunk, cmap=HEATMAP_CMAP, center=0, annot=False, 
                        cbar_kws={'label': r'ΔΔG (kcal/mol)', 'shrink': 0.8},
                        vmin=-3, vmax=5, square=True, ax=ax)
            
            part_info = f" (Part {chunk_idx+1}/{n_chunks})" if n_chunks > 1 else ""
            ax.set_title(f"{pdbname}\nChain {chain} Interface Mutation Landscape{part_info}", fontsize=12, fontweight='bold')
            ax.set_xlabel("Wild-type Residue Position", fontsize=10)
            ax.set_ylabel("Substituted Amino Acid", fontsize=10)
            plt.xticks(rotation=90, fontsize=7)
            plt.yticks(rotation=0, fontsize=8)
            plt.tight_layout()
            
            suffix = f"_part{chunk_idx+1}" if n_chunks > 1 else ""
            out_png = os.path.join(output_dir, f"1_interface_ddg_heatmap_chain_{chain}{suffix}.png")
            plt.savefig(out_png, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"[INFO] Saved interface heatmap for chain {chain}{suffix}")


def generate_folding_heatmaps(folding_df, output_dir, pdbname):
    """Generate folding stability ΔΔG heatmaps per chain with square cells."""
    for chain in folding_df['chain'].unique():
        chain_data = folding_df[folding_df['chain'] == chain].copy()
        chain_data = chain_data.sort_values('position')
        chain_data['pos_label'] = chain_data['original_aa'] + chain_data['position'].astype(str)
        
        pivot = chain_data.pivot_table(index="target_aa", columns="pos_label", values="ddg", aggfunc='mean')
        pivot = pivot.reindex(AA_ORDER)
        
        max_cols = 100
        n_chunks = (len(pivot.columns) + max_cols - 1) // max_cols
        
        for chunk_idx in range(n_chunks):
            start_col = chunk_idx * max_cols
            end_col = min((chunk_idx + 1) * max_cols, len(pivot.columns))
            pivot_chunk = pivot.iloc[:, start_col:end_col]
            
            # Calculate figure size for square cells
            cell_size = 0.4
            n_rows, n_cols = pivot_chunk.shape
            fig_width = n_cols * cell_size + 3
            fig_height = n_rows * cell_size + 2
            
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            sns.heatmap(pivot_chunk, cmap=HEATMAP_CMAP, center=0, annot=False, 
                        cbar_kws={'label': r'ΔΔG (kcal/mol)', 'shrink': 0.8},
                        vmin=-3, vmax=5, square=True, ax=ax)
            
            part_info = f" (Part {chunk_idx+1}/{n_chunks})" if n_chunks > 1 else ""
            ax.set_title(f"{pdbname}\nChain {chain} Folding Stability Mutation Landscape{part_info}", fontsize=12, fontweight='bold')
            ax.set_xlabel("Wild-type Residue Position", fontsize=10)
            ax.set_ylabel("Substituted Amino Acid", fontsize=10)
            plt.xticks(rotation=90, fontsize=7)
            plt.yticks(rotation=0, fontsize=8)
            plt.tight_layout()
            
            suffix = f"_part{chunk_idx+1}" if n_chunks > 1 else ""
            out_png = os.path.join(output_dir, f"2_folding_ddg_heatmap_chain_{chain}{suffix}.png")
            plt.savefig(out_png, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"[INFO] Saved folding heatmap for chain {chain}{suffix}")


def generate_chain_interface_profiles(interface_df, output_dir, pdbname):
    """Generate chain interaction analysis."""
    interface_summary = interface_df.groupby(['chain', 'position', 'original_aa']).agg({
        'ddg': ['mean', 'std', 'max', 'min']
    }).reset_index()
    interface_summary.columns = ['chain', 'position', 'original_aa', 'mean_ddg', 'std_ddg', 'max_ddg', 'min_ddg']
    
    fig, axes = plt.subplots(2, 1, figsize=(16, 12), sharex=False)
    
    for idx, chain in enumerate(sorted(interface_df['chain'].unique())):
        chain_summary = interface_summary[interface_summary['chain'] == chain].sort_values('position')
        x_labels = chain_summary['original_aa'] + chain_summary['position'].astype(str)
        
        axes[idx].bar(range(len(chain_summary)), chain_summary['mean_ddg'], 
                      yerr=chain_summary['std_ddg'], capsize=2, alpha=0.7,
                      color=['red' if x > 0.5 else 'blue' if x < -0.5 else 'gray' for x in chain_summary['mean_ddg']])
        axes[idx].axhline(y=0, color='black', linestyle='--', linewidth=1)
        axes[idx].axhline(y=0.5, color='red', linestyle=':', linewidth=1, alpha=0.5, label='Destabilizing threshold')
        axes[idx].axhline(y=-0.5, color='blue', linestyle=':', linewidth=1, alpha=0.5, label='Stabilizing threshold')
        axes[idx].set_xticks(range(len(chain_summary)))
        axes[idx].set_xticklabels(x_labels, rotation=90, fontsize=6)
        axes[idx].set_ylabel(r'FoldX $\Delta\Delta$G (kcal/mol)')
        axes[idx].set_title(f'Chain {chain} Interface Residues - Effect on Binding')
        axes[idx].legend(loc='upper right', fontsize=8)
        axes[idx].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    out_png = os.path.join(output_dir, "3_chain_interface_profiles.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Saved chain interface profiles")
    
    return interface_summary


def generate_hotspots(interface_summary, output_dir, pdbname, threshold=1.0):
    """Identify and visualize interface hotspots."""
    hotspots = interface_summary[abs(interface_summary['mean_ddg']) > threshold].copy()
    hotspots = hotspots.sort_values('mean_ddg', ascending=False)
    
    if len(hotspots) > 0:
        fig, ax = plt.subplots(figsize=(12, max(8, len(hotspots)*0.3)))
        
        labels = hotspots['chain'] + ':' + hotspots['original_aa'] + hotspots['position'].astype(str)
        colors = ['red' if x > 0 else 'blue' for x in hotspots['mean_ddg']]
        
        ax.barh(range(len(hotspots)), hotspots['mean_ddg'], color=colors, alpha=0.7)
        ax.set_yticks(range(len(hotspots)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
        ax.set_xlabel(r'Mean FoldX $\Delta\Delta$G (kcal/mol)')
        ax.set_title(f'{pdbname} Interface Hotspots (|ΔΔG| > {threshold} kcal/mol)\nRed: Destabilizing on mutation | Blue: Stabilizing on mutation')
        ax.grid(True, alpha=0.3, axis='x')
        
        plt.tight_layout()
        out_png = os.path.join(output_dir, "4_interface_hotspots.png")
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[INFO] Saved interface hotspots")


def generate_top_mutations(interface_df, output_dir, pdbname, top_n=30):
    """Generate top destabilizing/stabilizing mutations plot."""
    top_destab = interface_df.nlargest(top_n, 'ddg')
    top_stab = interface_df.nsmallest(top_n, 'ddg')
    top_combined = pd.concat([top_destab, top_stab]).reset_index(drop=True)
    
    fig, ax = plt.subplots(figsize=(12, 12))
    labels = top_combined['chain'] + ':' + top_combined['mutation_label']
    colors = ['red' if x > 0 else 'blue' for x in top_combined['ddg']]
    ax.barh(range(len(top_combined)), top_combined['ddg'], color=colors, alpha=0.7)
    ax.set_yticks(range(len(top_combined)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel(r'FoldX $\Delta\Delta$G (kcal/mol)')
    ax.set_title(f'{pdbname} Top {top_n} Destabilizing (Red) and Stabilizing (Blue) Mutations')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    out_png = os.path.join(output_dir, "5_top_mutations_bar.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Saved top mutations bar plot")


def generate_distributions(interface_df, folding_df, output_dir):
    """Generate distribution plots."""
    fig, axes = plt.subplots(1, 2 if folding_df is not None else 1, figsize=(15, 5))
    
    if folding_df is not None and not folding_df.empty:
        # Use viridis-inspired colors
        axes[0].hist(interface_df['ddg'], bins=50, alpha=0.7, color='#440154', edgecolor='black')
        axes[0].set_xlabel(r'FoldX $\Delta\Delta$G (kcal/mol)')
        axes[0].set_ylabel('Frequency')
        axes[0].set_title('Interface Stability Distribution')
        axes[0].axvline(x=interface_df['ddg'].mean(), color='#FDE725', linestyle='--', linewidth=2, label=f"Mean: {interface_df['ddg'].mean():.2f}")
        axes[0].axvline(x=0, color='black', linestyle='-', linewidth=1)
        axes[0].legend()
        
        axes[1].hist(folding_df['ddg'], bins=50, alpha=0.7, color='#21918C', edgecolor='black')
        axes[1].set_xlabel(r'FoldX $\Delta\Delta$G (kcal/mol)')
        axes[1].set_ylabel('Frequency')
        axes[1].set_title('Folding Stability Distribution')
        axes[1].axvline(x=folding_df['ddg'].mean(), color='#FDE725', linestyle='--', linewidth=2, label=f"Mean: {folding_df['ddg'].mean():.2f}")
        axes[1].axvline(x=0, color='black', linestyle='-', linewidth=1)
        axes[1].legend()
    else:
        if hasattr(axes, '__iter__'):
            ax = axes[0]
        else:
            ax = axes
        ax.hist(interface_df['ddg'], bins=50, alpha=0.7, color='#440154', edgecolor='black')
        ax.set_xlabel(r'FoldX $\Delta\Delta$G (kcal/mol)')
        ax.set_ylabel('Frequency')
        ax.set_title('Interface Stability Distribution')
        ax.axvline(x=interface_df['ddg'].mean(), color='#FDE725', linestyle='--', linewidth=2, label=f"Mean: {interface_df['ddg'].mean():.2f}")
        ax.legend()
    
    plt.tight_layout()
    out_png = os.path.join(output_dir, "6_distribution_plot.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Saved distribution plot")


def generate_chain_comparison(interface_df, output_dir):
    """Generate chain contribution comparison."""
    chain_comparison = interface_df.groupby('chain').agg({
        'ddg': ['mean', 'std', 'sum', 'count']
    }).reset_index()
    chain_comparison.columns = ['chain', 'mean_ddg', 'std_ddg', 'sum_ddg', 'count']
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Use plasma-inspired colors
    axes[0].bar(chain_comparison['chain'], chain_comparison['mean_ddg'], 
                yerr=chain_comparison['std_ddg'], capsize=5, color=['#0D0887', '#F0F921'][:len(chain_comparison)])
    axes[0].set_xlabel('Chain')
    axes[0].set_ylabel(r'Mean FoldX $\Delta\Delta$G (kcal/mol)')
    axes[0].set_title('Average Interface Effect by Chain')
    axes[0].axhline(y=0, color='black', linestyle='--')
    
    axes[1].bar(chain_comparison['chain'], chain_comparison['count'] / 20,
                color=['#0D0887', '#F0F921'][:len(chain_comparison)])
    axes[1].set_xlabel('Chain')
    axes[1].set_ylabel('Number of Interface Residues')
    axes[1].set_title('Interface Residues per Chain')
    
    plt.tight_layout()
    out_png = os.path.join(output_dir, "7_chain_comparison.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Saved chain comparison")


def generate_interface_surface_heatmap(interface_df, output_dir, pdbname, top_n=50):
    """Generate a heatmap showing the most promising interacting surface between chains.
    
    This creates a matrix where rows are top interface residues from one chain and
    columns are from another chain, with values showing their combined binding importance.
    """
    chains = sorted(interface_df['chain'].unique())
    
    if len(chains) < 2:
        print("[WARNING] Need at least 2 chains for interface surface heatmap")
        return
    
    # Get chain pair information to identify which residues interact across chains
    chain_pairs = interface_df['chain_pair'].unique()
    
    # Calculate mean ΔΔG per position (across all mutations)
    position_summary = interface_df.groupby(['chain', 'position', 'original_aa', 'position_label']).agg({
        'ddg': ['mean', 'max', 'std']
    }).reset_index()
    position_summary.columns = ['chain', 'position', 'original_aa', 'position_label', 'mean_ddg', 'max_ddg', 'std_ddg']
    
    # Calculate an importance score (absolute mean ΔΔG indicates binding importance)
    position_summary['importance'] = abs(position_summary['mean_ddg'])
    
    # Get top residues from each chain
    chain_a = chains[0]
    chain_b = chains[1]
    
    top_a = position_summary[position_summary['chain'] == chain_a].nlargest(top_n, 'importance')
    top_b = position_summary[position_summary['chain'] == chain_b].nlargest(top_n, 'importance')
    
    if len(top_a) == 0 or len(top_b) == 0:
        print("[WARNING] Insufficient data for interface surface heatmap")
        return
    
    # Create position labels
    top_a = top_a.sort_values('position')
    top_b = top_b.sort_values('position')
    
    labels_a = (top_a['original_aa'] + top_a['position'].astype(str)).tolist()
    labels_b = (top_b['original_aa'] + top_b['position'].astype(str)).tolist()
    
    # Create interaction matrix based on combined importance scores
    # Each cell (i,j) shows the geometric mean of importance of residue i from chain A
    # and residue j from chain B - representing their combined interface contribution
    interaction_matrix = np.zeros((len(labels_a), len(labels_b)))
    
    importance_a = top_a['importance'].values
    importance_b = top_b['importance'].values
    mean_ddg_a = top_a['mean_ddg'].values
    mean_ddg_b = top_b['mean_ddg'].values
    
    for i in range(len(labels_a)):
        for j in range(len(labels_b)):
            # Combined importance: geometric mean of absolute effects
            combined = np.sqrt(importance_a[i] * importance_b[j])
            # Sign based on whether both destabilize (positive) or stabilize (negative)
            if mean_ddg_a[i] > 0 and mean_ddg_b[j] > 0:
                interaction_matrix[i, j] = combined  # Both destabilizing
            elif mean_ddg_a[i] < 0 and mean_ddg_b[j] < 0:
                interaction_matrix[i, j] = -combined  # Both stabilizing
            else:
                interaction_matrix[i, j] = combined * 0.5  # Mixed effect
    
    # Create the heatmap with square cells
    cell_size = 0.5
    fig_width = len(labels_b) * cell_size + 4
    fig_height = len(labels_a) * cell_size + 3
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # Create DataFrame for seaborn
    interaction_df = pd.DataFrame(interaction_matrix, index=labels_a, columns=labels_b)
    
    sns.heatmap(interaction_df, cmap=HEATMAP_CMAP, center=0, annot=False,
                cbar_kws={'label': r'Combined Importance Score', 'shrink': 0.8},
                xticklabels=True, yticklabels=True, ax=ax, square=True)
    
    ax.set_xlabel(f'Chain {chain_b} Interface Residues', fontsize=11, fontweight='bold')
    ax.set_ylabel(f'Chain {chain_a} Interface Residues', fontsize=11, fontweight='bold')
    ax.set_title(f'{pdbname}\nInter-Chain Interface Hotspot Map\nTop {min(top_n, len(top_a))} × {min(top_n, len(top_b))} Most Significant Residues', fontsize=12, fontweight='bold')
    
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()
    
    out_png = os.path.join(output_dir, "8_interface_surface_heatmap.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Saved interface surface heatmap")
    
    # Also create a focused view with only the most critical residues (top 20)
    top_critical = 20
    if len(top_a) >= 10 and len(top_b) >= 10:
        top_a_critical = position_summary[position_summary['chain'] == chain_a].nlargest(top_critical, 'importance').sort_values('position')
        top_b_critical = position_summary[position_summary['chain'] == chain_b].nlargest(top_critical, 'importance').sort_values('position')
        
        labels_a_c = (top_a_critical['original_aa'] + top_a_critical['position'].astype(str)).tolist()
        labels_b_c = (top_b_critical['original_aa'] + top_b_critical['position'].astype(str)).tolist()
        
        importance_a_c = top_a_critical['importance'].values
        importance_b_c = top_b_critical['importance'].values
        mean_ddg_a_c = top_a_critical['mean_ddg'].values
        mean_ddg_b_c = top_b_critical['mean_ddg'].values
        
        critical_matrix = np.zeros((len(labels_a_c), len(labels_b_c)))
        for i in range(len(labels_a_c)):
            for j in range(len(labels_b_c)):
                combined = np.sqrt(importance_a_c[i] * importance_b_c[j])
                if mean_ddg_a_c[i] > 0 and mean_ddg_b_c[j] > 0:
                    critical_matrix[i, j] = combined
                elif mean_ddg_a_c[i] < 0 and mean_ddg_b_c[j] < 0:
                    critical_matrix[i, j] = -combined
                else:
                    critical_matrix[i, j] = combined * 0.5
        
        # Square cells for critical contacts heatmap
        cell_size = 0.6
        fig_width = len(labels_b_c) * cell_size + 4
        fig_height = len(labels_a_c) * cell_size + 3
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        critical_df = pd.DataFrame(critical_matrix, index=labels_a_c, columns=labels_b_c)
        
        sns.heatmap(critical_df, cmap=HEATMAP_CMAP, center=0, annot=True, fmt='.2f',
                    annot_kws={'size': 7},
                    cbar_kws={'label': r'Combined Importance Score', 'shrink': 0.8},
                    xticklabels=True, yticklabels=True, ax=ax, square=True)
        
        ax.set_xlabel(f'Chain {chain_b} Critical Residues', fontsize=11, fontweight='bold')
        ax.set_ylabel(f'Chain {chain_a} Critical Residues', fontsize=11, fontweight='bold')
        ax.set_title(f'{pdbname}\nCritical Inter-Chain Contacts\nTop {len(labels_a_c)} × {len(labels_b_c)} Binding Hotspots', fontsize=12, fontweight='bold')
        
        plt.xticks(rotation=45, ha='right', fontsize=9)
        plt.yticks(rotation=0, fontsize=9)
        plt.tight_layout()
        
        out_png = os.path.join(output_dir, "9_critical_interface_contacts.png")
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[INFO] Saved critical interface contacts heatmap")


def generate_selfmutation_heatmap(df, output_dir, pdbname, data_type='interface'):
    """Generate a bar chart showing selfmutation ΔΔG values across all residues.
    
    For selfmutation data, we show residue positions with ΔΔG values as bars.
    """
    if df.empty:
        return
    
    for chain in df['chain'].unique():
        chain_data = df[df['chain'] == chain].copy()
        chain_data = chain_data.sort_values('position')
        chain_data['pos_label'] = chain_data['original_aa'] + chain_data['position'].astype(str)
        
        # Create bar chart showing ΔΔG for each position
        max_cols = 100
        positions = chain_data['pos_label'].tolist()
        n_chunks = (len(positions) + max_cols - 1) // max_cols
        
        for chunk_idx in range(n_chunks):
            start_idx = chunk_idx * max_cols
            end_idx = min((chunk_idx + 1) * max_cols, len(positions))
            chunk_data = chain_data.iloc[start_idx:end_idx]
            
            # Auto-adjust figure width based on number of residues
            bar_width = 0.15
            fig_width = max(12, len(chunk_data) * bar_width + 2)
            fig_height = 6
            
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            
            colors = ['#d62728' if x > 0.5 else '#1f77b4' if x < -0.5 else '#7f7f7f' 
                      for x in chunk_data['ddg']]
            
            bars = ax.bar(range(len(chunk_data)), chunk_data['ddg'], 
                         yerr=chunk_data['ddg_std'], capsize=1, alpha=0.8,
                         color=colors, edgecolor='black', linewidth=0.3)
            
            ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
            ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='Destabilizing (>0.5)')
            ax.axhline(y=-0.5, color='blue', linestyle='--', alpha=0.5, label='Stabilizing (<-0.5)')
            
            ax.set_xticks(range(len(chunk_data)))
            ax.set_xticklabels(chunk_data['pos_label'], rotation=90, fontsize=6)
            ax.set_ylabel(r'ΔΔG (kcal/mol)', fontsize=10, fontweight='bold')
            ax.set_xlabel('Residue Position (AA + Chain + Number)', fontsize=10, fontweight='bold')
            
            prefix = "Interface Binding" if data_type == 'interface' else "Protein Folding"
            part_info = f" (Part {chunk_idx+1}/{n_chunks})" if n_chunks > 1 else ""
            ax.set_title(f'{pdbname}\nChain {chain} - {prefix} Stability Profile{part_info}', 
                        fontsize=12, fontweight='bold')
            ax.legend(loc='upper right', fontsize=8)
            ax.grid(True, alpha=0.3, axis='y')
            
            plt.tight_layout()
            
            part_suffix = f"_part{chunk_idx+1}" if n_chunks > 1 else ""
            type_prefix = "1_interface" if data_type == 'interface' else "2_folding"
            out_png = os.path.join(output_dir, f"{type_prefix}_selfmut_chain_{chain}{part_suffix}.png")
            plt.savefig(out_png, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"[INFO] Saved {data_type} selfmutation plot for chain {chain}{part_suffix}")


def generate_selfmutation_chain_comparison(interface_df, folding_df, output_dir, pdbname):
    """Generate comparison of interface vs folding effects by chain."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    chains = sorted(interface_df['chain'].unique())
    
    # Mean ΔΔG by chain
    ax = axes[0]
    x = np.arange(len(chains))
    width = 0.35
    
    interface_means = [interface_df[interface_df['chain'] == c]['ddg'].mean() for c in chains]
    interface_stds = [interface_df[interface_df['chain'] == c]['ddg'].std() for c in chains]
    
    bars1 = ax.bar(x - width/2, interface_means, width, yerr=interface_stds, 
                   label='Interface Binding', color='#440154', capsize=5)
    
    if folding_df is not None and not folding_df.empty:
        folding_means = [folding_df[folding_df['chain'] == c]['ddg'].mean() for c in chains]
        folding_stds = [folding_df[folding_df['chain'] == c]['ddg'].std() for c in chains]
        bars2 = ax.bar(x + width/2, folding_means, width, yerr=folding_stds,
                       label='Protein Folding', color='#21918C', capsize=5)
    
    ax.set_ylabel(r'Mean ΔΔG (kcal/mol)', fontsize=10, fontweight='bold')
    ax.set_xlabel('Protein Chain', fontsize=10, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(chains, fontsize=11)
    ax.set_title(f'{pdbname}\nMean Stability Effect by Chain', fontsize=12, fontweight='bold')
    ax.axhline(y=0, color='black', linestyle='--')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Residue count by chain
    ax = axes[1]
    interface_counts = [len(interface_df[interface_df['chain'] == c]) for c in chains]
    
    bars = ax.bar(x, interface_counts, color=['#0D0887', '#F0F921'][:len(chains)])
    ax.set_ylabel('Number of Residues', fontsize=10, fontweight='bold')
    ax.set_xlabel('Protein Chain', fontsize=10, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(chains, fontsize=11)
    ax.set_title(f'{pdbname}\nResidue Count per Chain', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add count labels on bars
    for i, (bar, count) in enumerate(zip(bars, interface_counts)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, 
                str(count), ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    out_png = os.path.join(output_dir, "3_chain_comparison.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Saved chain comparison plot")


def generate_selfmutation_distributions(interface_df, folding_df, output_dir, pdbname):
    """Generate distribution plots for selfmutation data."""
    n_plots = 2 if (folding_df is not None and not folding_df.empty) else 1
    fig, axes = plt.subplots(1, n_plots, figsize=(7*n_plots, 5))
    
    if n_plots == 1:
        axes = [axes]
    
    # Interface distribution
    ax = axes[0]
    ax.hist(interface_df['ddg'], bins=50, alpha=0.7, color='#440154', edgecolor='black')
    ax.axvline(x=interface_df['ddg'].mean(), color='#FDE725', linestyle='--', linewidth=2,
               label=f"Mean: {interface_df['ddg'].mean():.4f} kcal/mol")
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_xlabel(r'ΔΔG (kcal/mol)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Frequency (Number of Residues)', fontsize=10, fontweight='bold')
    ax.set_title(f'{pdbname}\nInterface Binding Energy Distribution', fontsize=12, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    
    # Folding distribution
    if n_plots == 2:
        ax = axes[1]
        ax.hist(folding_df['ddg'], bins=50, alpha=0.7, color='#21918C', edgecolor='black')
        ax.axvline(x=folding_df['ddg'].mean(), color='#FDE725', linestyle='--', linewidth=2,
                   label=f"Mean: {folding_df['ddg'].mean():.4f} kcal/mol")
        ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
        ax.set_xlabel(r'ΔΔG (kcal/mol)', fontsize=10, fontweight='bold')
        ax.set_ylabel('Frequency (Number of Residues)', fontsize=10, fontweight='bold')
        ax.set_title(f'{pdbname}\nFolding Stability Energy Distribution', fontsize=12, fontweight='bold')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    out_png = os.path.join(output_dir, "4_distribution.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Saved distribution plot")


def generate_selfmutation_hotspots(df, output_dir, pdbname, threshold=0.1, data_type='interface'):
    """Identify and plot residues with significant ΔΔG values."""
    # For selfmutation, most values should be ~0, so highlight deviations
    significant = df[abs(df['ddg']) > threshold].copy()
    significant = significant.sort_values('ddg', ascending=True)
    
    if len(significant) == 0:
        print(f"[INFO] No significant {data_type} hotspots found (threshold: {threshold})")
        return
    
    # Limit to top 50 for readability
    if len(significant) > 50:
        top_destab = significant.nlargest(25, 'ddg')
        top_stab = significant.nsmallest(25, 'ddg')
        significant = pd.concat([top_stab, top_destab])
    
    # Auto-adjust figure height based on number of residues
    bar_height = 0.3
    fig_height = max(8, len(significant) * bar_height + 2)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    
    labels = significant['chain'] + ':' + significant['original_aa'] + significant['position'].astype(str)
    colors = ['#d62728' if x > 0 else '#1f77b4' for x in significant['ddg']]
    
    ax.barh(range(len(significant)), significant['ddg'], 
            xerr=significant['ddg_std'], color=colors, alpha=0.7, capsize=3)
    ax.set_yticks(range(len(significant)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.set_xlabel(r'ΔΔG (kcal/mol)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Residue (Chain:AA+Position)', fontsize=10, fontweight='bold')
    
    prefix = "Interface Binding" if data_type == 'interface' else "Protein Folding"
    ax.set_title(f'{pdbname}\n{prefix} Stability Hotspots (|ΔΔG| > {threshold} kcal/mol)\n'
                 f'Red = Destabilizing Effect | Blue = Stabilizing Effect', 
                 fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    type_prefix = "5_interface" if data_type == 'interface' else "6_folding"
    out_png = os.path.join(output_dir, f"{type_prefix}_hotspots.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Saved {data_type} hotspots plot")


def generate_interface_vs_folding_scatter(interface_df, folding_df, output_dir, pdbname):
    """Generate scatter plot comparing interface vs folding ΔΔG."""
    if folding_df is None or folding_df.empty:
        return
    
    # Merge on position and chain
    merged = interface_df.merge(folding_df, on=['chain', 'position', 'original_aa'],
                                 suffixes=('_interface', '_folding'))
    
    if merged.empty:
        print("[WARNING] No matching residues for interface vs folding comparison")
        return
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Color by chain instead of position for better interpretation
    chains = merged['chain'].unique()
    colors = plt.cm.Set1(np.linspace(0, 1, len(chains)))
    chain_colors = {c: colors[i] for i, c in enumerate(chains)}
    
    for chain in chains:
        chain_data = merged[merged['chain'] == chain]
        ax.scatter(chain_data['ddg_interface'], chain_data['ddg_folding'],
                  c=[chain_colors[chain]], alpha=0.6, s=50, label=f'Chain {chain}')
    
    # Add diagonal line
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', alpha=0.5, label='Equal Effect (x=y)')
    
    ax.set_xlabel(r'Interface Binding ΔΔG (kcal/mol)', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'Protein Folding ΔΔG (kcal/mol)', fontsize=11, fontweight='bold')
    ax.set_title(f'{pdbname}\nInterface Binding vs Folding Stability Correlation\n'
                 f'Each point = one residue position', fontsize=12, fontweight='bold')
    ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle=':', alpha=0.5)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Make axes equal for fair comparison
    ax.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    
    out_png = os.path.join(output_dir, "7_interface_vs_folding.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Saved interface vs folding scatter plot")


def generate_interpretation_guide(viz_dir, pdbname, is_selfmutation=True):
    """Generate a text file explaining how to interpret the visualization results."""
    
    guide_content = f"""================================================================================
MUTATEX VISUALIZATION INTERPRETATION GUIDE
================================================================================
Structure: {pdbname}
Generated by: generate_mutatex_visualizations.py
================================================================================

UNDERSTANDING ΔΔG (Delta-Delta-G) VALUES
----------------------------------------
ΔΔG represents the change in free energy upon mutation, measured in kcal/mol.

• ΔΔG > 0 (POSITIVE): DESTABILIZING mutation
  - The mutation weakens the protein structure or binding interface
  - Higher positive values = more destabilizing effect
  - Threshold: ΔΔG > 0.5 kcal/mol is typically considered significant

• ΔΔG < 0 (NEGATIVE): STABILIZING mutation  
  - The mutation strengthens the protein structure or binding interface
  - More negative values = more stabilizing effect
  - Threshold: ΔΔG < -0.5 kcal/mol is typically considered significant

• ΔΔG ≈ 0: NEUTRAL mutation
  - The mutation has minimal effect on stability

================================================================================
FIGURE DESCRIPTIONS AND INTERPRETATION
================================================================================

"""
    
    if is_selfmutation:
        guide_content += """
1. INTERFACE STABILITY PROFILE (1_interface_selfmut_chain_*.png)
----------------------------------------------------------------
Location: heatmaps/
What it shows: Bar chart of ΔΔG values for each residue position in the chain.
X-axis: Residue positions (format: amino acid letter + position number)
Y-axis: ΔΔG in kcal/mol
Color coding:
  • RED bars: Destabilizing effect (ΔΔG > 0.5)
  • BLUE bars: Stabilizing effect (ΔΔG < -0.5)  
  • GRAY bars: Neutral effect (-0.5 ≤ ΔΔG ≤ 0.5)
Error bars: Standard deviation across FoldX runs
How to interpret: 
  - Tall red bars indicate residues critical for interface binding
  - These are potential drug targets or mutation-sensitive sites
  - Multiple parts may exist if the chain has many residues


2. FOLDING STABILITY PROFILE (2_folding_selfmut_chain_*.png)
------------------------------------------------------------
Location: heatmaps/
What it shows: Bar chart of ΔΔG values for protein folding stability.
Same format as interface profile but measures effect on protein fold stability.
How to interpret:
  - Residues with high ΔΔG are structurally important for the fold
  - Compare with interface data to identify binding-specific vs structural residues


3. CHAIN COMPARISON (3_chain_comparison.png)
--------------------------------------------
Location: analysis_plots/
Left panel - Mean ΔΔG by Chain:
  • Compares average stability effects between chains
  • Error bars show variation across all residues
  • Purple = Interface binding effect, Teal = Folding stability effect
Right panel - Residue Count:
  • Shows number of residues analyzed per chain
  • Useful for understanding data coverage
How to interpret:
  - Chains with higher mean ΔΔG contribute more to binding/stability
  - Larger error bars indicate more heterogeneous contribution


4. ENERGY DISTRIBUTION (4_distribution.png)
-------------------------------------------
Location: analysis_plots/
What it shows: Histogram of all ΔΔG values.
X-axis: ΔΔG bins (kcal/mol)
Y-axis: Frequency (number of residues)
Yellow dashed line: Mean ΔΔG value
How to interpret:
  - Most residues should cluster near zero (neutral effect)
  - Long tails toward positive values indicate destabilizing hotspots
  - Distribution shape reveals overall stability landscape


5. INTERFACE HOTSPOTS (5_interface_hotspots.png)
------------------------------------------------
Location: analysis_plots/
What it shows: Horizontal bar chart of residues with significant ΔΔG.
Y-axis: Residue labels (Chain:AminoAcid+Position)
X-axis: ΔΔG value
Color coding:
  • RED: Destabilizing (positive ΔΔG)
  • BLUE: Stabilizing (negative ΔΔG)
How to interpret:
  - Top destabilizing residues are key binding hotspots
  - These residues are critical for protein-protein interaction
  - Potential targets for interface disruption or stabilization


6. FOLDING HOTSPOTS (6_folding_hotspots.png)
--------------------------------------------
Location: analysis_plots/
What it shows: Same as interface hotspots but for folding stability.
How to interpret:
  - Identifies residues critical for protein structural integrity
  - Compare with interface hotspots to distinguish binding vs structural roles


7. INTERFACE VS FOLDING CORRELATION (7_interface_vs_folding.png)
----------------------------------------------------------------
Location: analysis_plots/
What it shows: Scatter plot comparing interface vs folding ΔΔG for each residue.
X-axis: Interface binding ΔΔG
Y-axis: Protein folding ΔΔG
Diagonal dashed line: Equal effect (x=y)
Points colored by chain
How to interpret:
  - Points above diagonal: More effect on folding than binding
  - Points below diagonal: More effect on binding than folding
  - Points near origin: Neutral residues
  - Upper-right quadrant: Destabilizing for both
  - Lower-left quadrant: Stabilizing for both
  - Off-diagonal points: Binding-specific or folding-specific effects

"""
    else:
        guide_content += """
1. INTERFACE ΔΔG HEATMAP (1_interface_ddg_heatmap_chain_*.png)
--------------------------------------------------------------
Location: heatmaps/interface/
What it shows: Matrix of mutation effects for each position.
X-axis: Wild-type residue positions (format: AA letter + position)
Y-axis: All 20 amino acid substitutions (in standard order)
Color scale: Purple-Yellow (plasma colormap)
  • Dark purple/blue: Low ΔΔG (stabilizing or neutral)
  • Yellow/orange: High ΔΔG (destabilizing)
How to interpret:
  - Each column shows all possible mutations at that position
  - Hot (yellow) columns indicate mutation-sensitive positions
  - Cool (purple) columns indicate mutation-tolerant positions
  - Horizontal patterns reveal amino acid-specific effects


2. FOLDING ΔΔG HEATMAP (2_folding_ddg_heatmap_chain_*.png)
----------------------------------------------------------
Location: heatmaps/folding/
Same format as interface heatmap but for folding stability.
How to interpret:
  - Compare with interface heatmap to distinguish binding vs structural roles


3. CHAIN INTERFACE PROFILES (3_chain_interface_profiles.png)
------------------------------------------------------------
Location: analysis_plots/
What it shows: Mean ΔΔG per position across all mutations.
How to interpret:
  - Peaks indicate mutation-sensitive positions
  - Valleys indicate mutation-tolerant positions


4. INTERFACE HOTSPOTS (4_interface_hotspots.png)
------------------------------------------------
Location: analysis_plots/
What it shows: Ranked list of mutation-sensitive positions.
How to interpret:
  - Top positions are critical for binding


5. TOP MUTATIONS (5_top_mutations_bar.png)
------------------------------------------
Location: analysis_plots/
What it shows: Most destabilizing and stabilizing specific mutations.
How to interpret:
  - Red bars: Mutations to avoid (destabilizing)
  - Blue bars: Potential beneficial mutations (stabilizing)


6. ENERGY DISTRIBUTIONS (6_distribution_plot.png)
-------------------------------------------------
Location: analysis_plots/
What it shows: Overall distribution of mutation effects.


7. CHAIN COMPARISON (7_chain_comparison.png)
--------------------------------------------
Location: analysis_plots/
What it shows: Comparison of stability contributions between chains.


8-9. INTERFACE SURFACE HEATMAPS (8_interface_surface_heatmap.png, 9_critical_interface_contacts.png)
----------------------------------------------------------------------------------------------------
Location: analysis_plots/
What it shows: Combined importance of residue pairs across the interface.
How to interpret:
  - Identifies which residue pairs form critical contacts
  - High values indicate important inter-chain interactions

"""

    guide_content += """
================================================================================
DATA FILES
================================================================================
Location: data/

CSV files contain the raw numerical data used to generate the plots.
Columns typically include:
  • original_aa: Wild-type amino acid
  • position: Residue number
  • chain: Protein chain identifier
  • ddg: ΔΔG value (kcal/mol)
  • ddg_std: Standard deviation
  • ddg_min, ddg_max: Range of values

Use these files for:
  - Custom analysis and filtering
  - Integration with structural visualization tools
  - Statistical analysis
  - Comparison with other datasets


================================================================================
PRACTICAL APPLICATIONS
================================================================================

1. DRUG DESIGN:
   - Target interface hotspots (high positive ΔΔG)
   - Design peptides mimicking stabilizing interactions

2. PROTEIN ENGINEERING:
   - Avoid mutating structurally critical residues
   - Introduce stabilizing mutations (negative ΔΔG)

3. DISEASE MUTATION ANALYSIS:
   - Check if disease mutations hit hotspot residues
   - Predict severity based on ΔΔG magnitude

4. BINDING SITE IDENTIFICATION:
   - Interface hotspots define the binding epitope
   - Compare chains to identify dominant binding partner


================================================================================
REFERENCES
================================================================================
• FoldX force field: Schymkowitz et al., Nucleic Acids Res, 2005
• MutateX pipeline: Delgado et al., Bioinformatics, 2019
• ΔΔG interpretation: Tokuriki & Tawfik, Curr Opin Struct Biol, 2009

================================================================================
"""

    guide_path = os.path.join(viz_dir, "INTERPRETATION_GUIDE.txt")
    with open(guide_path, 'w') as f:
        f.write(guide_content)
    
    print(f"[INFO] Saved interpretation guide")
    return guide_path


def process_single_result(results_dir, pdbname=None):
    """Process a single results directory and generate visualizations."""
    if pdbname is None:
        pdbname = os.path.basename(results_dir)
    
    print(f"\n{'='*60}")
    print(f"[INFO] Processing: {pdbname}")
    print(f"[INFO] Results directory: {results_dir}")
    print(f"{'='*60}")
    
    # Create visualization output directory structure
    viz_dir = os.path.join(results_dir, "visualizations")
    heatmaps_dir = os.path.join(viz_dir, "heatmaps")
    analysis_dir = os.path.join(viz_dir, "analysis_plots")
    data_dir = os.path.join(viz_dir, "data")
    
    os.makedirs(heatmaps_dir, exist_ok=True)
    os.makedirs(analysis_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)
    
    # Try selfmutation data format first (the actual data format present)
    print("[INFO] Looking for selfmutation data...")
    interface_df, folding_df = parse_selfmutation_results(results_dir)
    
    if not interface_df.empty:
        print(f"[INFO] Found selfmutation data: {len(interface_df)} interface residues")
        if not folding_df.empty:
            print(f"[INFO] Found {len(folding_df)} folding residues")
        
        # Set style
        sns.set_context("talk")
        
        # Generate selfmutation-specific visualizations
        print("[INFO] Generating interface heatmaps...")
        generate_selfmutation_heatmap(interface_df, heatmaps_dir, pdbname, 'interface')
        
        if not folding_df.empty:
            print("[INFO] Generating folding heatmaps...")
            generate_selfmutation_heatmap(folding_df, heatmaps_dir, pdbname, 'folding')
        
        print("[INFO] Generating chain comparison...")
        generate_selfmutation_chain_comparison(interface_df, folding_df, analysis_dir, pdbname)
        
        print("[INFO] Generating distributions...")
        generate_selfmutation_distributions(interface_df, folding_df, analysis_dir, pdbname)
        
        print("[INFO] Generating interface hotspots...")
        generate_selfmutation_hotspots(interface_df, analysis_dir, pdbname, threshold=0.01, data_type='interface')
        
        if not folding_df.empty:
            print("[INFO] Generating folding hotspots...")
            generate_selfmutation_hotspots(folding_df, analysis_dir, pdbname, threshold=0.01, data_type='folding')
            
            print("[INFO] Generating interface vs folding comparison...")
            generate_interface_vs_folding_scatter(interface_df, folding_df, analysis_dir, pdbname)
        
        # Save CSV summaries
        print("[INFO] Saving CSV summaries...")
        interface_csv = os.path.join(data_dir, "interface_selfmutation_summary.csv")
        interface_df.to_csv(interface_csv, index=False)
        
        if not folding_df.empty:
            folding_csv = os.path.join(data_dir, "folding_selfmutation_summary.csv")
            folding_df.to_csv(folding_csv, index=False)
        
        # Generate interpretation guide
        print("[INFO] Generating interpretation guide...")
        generate_interpretation_guide(viz_dir, pdbname, is_selfmutation=True)
        
        # Print summary
        n_heatmaps = len([f for f in os.listdir(heatmaps_dir) if f.endswith('.png')])
        n_analysis = len([f for f in os.listdir(analysis_dir) if f.endswith('.png')])
        n_data = len([f for f in os.listdir(data_dir) if f.endswith('.csv')])
        
        print(f"[SUCCESS] Visualizations generated for {pdbname}!")
        print(f"  Output: {viz_dir}/")
        print(f"    ├── heatmaps/       ({n_heatmaps} files)")
        print(f"    ├── analysis_plots/ ({n_analysis} files)")
        print(f"    └── data/           ({n_data} files)")
        
        return True
    
    # Fall back to full mutation scan format
    interface_dir = os.path.join(results_dir, "results", "interface_ddgs", "final_averages")
    folding_dir = os.path.join(results_dir, "results", "mutation_ddgs", "final_averages")
    
    if os.path.isdir(interface_dir):
        print("[INFO] Parsing full mutation scan data...")
        interface_df = parse_interface_data(interface_dir)
        
        if not interface_df.empty:
            print(f"[INFO] Parsed {len(interface_df)} mutations")
            
            folding_df = None
            if os.path.isdir(folding_dir):
                folding_df = parse_folding_data(folding_dir)
            
            # Set style
            sns.set_context("talk")
            
            # Create subdirectories
            heatmaps_interface_dir = os.path.join(heatmaps_dir, "interface")
            heatmaps_folding_dir = os.path.join(heatmaps_dir, "folding")
            os.makedirs(heatmaps_interface_dir, exist_ok=True)
            os.makedirs(heatmaps_folding_dir, exist_ok=True)
            
            # Generate full mutation visualizations
            print("[INFO] Generating interface heatmaps...")
            generate_interface_heatmaps(interface_df, heatmaps_interface_dir, pdbname)
            
            if folding_df is not None and not folding_df.empty:
                print("[INFO] Generating folding heatmaps...")
                generate_folding_heatmaps(folding_df, heatmaps_folding_dir, pdbname)
            
            print("[INFO] Generating chain interface profiles...")
            interface_summary = generate_chain_interface_profiles(interface_df, analysis_dir, pdbname)
            
            print("[INFO] Generating hotspots...")
            generate_hotspots(interface_summary, analysis_dir, pdbname)
            
            print("[INFO] Generating top mutations...")
            generate_top_mutations(interface_df, analysis_dir, pdbname)
            
            print("[INFO] Generating distribution plots...")
            generate_distributions(interface_df, folding_df, analysis_dir)
            
            print("[INFO] Generating chain comparison...")
            generate_chain_comparison(interface_df, analysis_dir)
            
            print("[INFO] Generating interface surface heatmap...")
            generate_interface_surface_heatmap(interface_df, analysis_dir, pdbname)
            
            # Save CSV summaries
            print("[INFO] Saving CSV summaries...")
            interface_csv = os.path.join(data_dir, "interface_mutations_summary.csv")
            interface_df.to_csv(interface_csv, index=False)
            
            if folding_df is not None and not folding_df.empty:
                folding_csv = os.path.join(data_dir, "folding_mutations_summary.csv")
                folding_df.to_csv(folding_csv, index=False)
            
            # Generate interpretation guide
            print("[INFO] Generating interpretation guide...")
            generate_interpretation_guide(viz_dir, pdbname, is_selfmutation=False)
            
            print(f"[SUCCESS] Visualizations generated for {pdbname}!")
            return True
    
    print(f"[WARNING] No valid data found in {results_dir}")
    return False


def discover_results_directories(base_dir):
    """Discover all result directories in the run_results folder."""
    results_dirs = []
    
    if not os.path.isdir(base_dir):
        return results_dirs
    
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path):
            # Check if it looks like a results directory
            # (has results_selfmutation or results subdirectory, or has .completed marker)
            if (os.path.isdir(os.path.join(item_path, "results_selfmutation")) or
                os.path.isdir(os.path.join(item_path, "results")) or
                os.path.isfile(os.path.join(item_path, ".completed"))):
                results_dirs.append(item_path)
    
    return sorted(results_dirs)


def main():
    # Determine script directory for default paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_run_results = os.path.join(script_dir, DEFAULT_RUN_RESULTS_DIR)
    
    if len(sys.argv) >= 3:
        # Legacy mode: specific results_dir and pdb_name
        results_dir = sys.argv[1]
        pdbname = sys.argv[2]
        process_single_result(results_dir, pdbname)
    
    elif len(sys.argv) == 2:
        # Single results directory specified
        results_dir = sys.argv[1]
        if os.path.isdir(results_dir):
            process_single_result(results_dir)
        else:
            print(f"[ERROR] Directory not found: {results_dir}")
            sys.exit(1)
    
    else:
        # Auto-discover mode: process all results in run_results/
        print(f"[INFO] Auto-discovery mode: looking for results in {default_run_results}")
        
        if not os.path.isdir(default_run_results):
            print(f"[ERROR] Run results directory not found: {default_run_results}")
            print("\nUsage:")
            print("  python generate_mutatex_visualizations.py                    # Process all in run_results/")
            print("  python generate_mutatex_visualizations.py <results_dir>      # Process single directory")
            print("  python generate_mutatex_visualizations.py <results_dir> <name>  # Legacy mode")
            sys.exit(1)
        
        results_dirs = discover_results_directories(default_run_results)
        
        if not results_dirs:
            print(f"[ERROR] No result directories found in {default_run_results}")
            sys.exit(1)
        
        print(f"[INFO] Found {len(results_dirs)} result directories:")
        for rd in results_dirs:
            print(f"  - {os.path.basename(rd)}")
        
        successful = 0
        failed = 0
        
        for results_dir in results_dirs:
            try:
                if process_single_result(results_dir):
                    successful += 1
                else:
                    failed += 1
            except Exception as e:
                print(f"[ERROR] Failed to process {results_dir}: {e}")
                failed += 1
        
        print(f"\n{'='*60}")
        print(f"[SUMMARY] Processed {successful + failed} directories:")
        print(f"  - Successful: {successful}")
        print(f"  - Failed: {failed}")
        print(f"{'='*60}")


if __name__ == "__main__":
    main()
