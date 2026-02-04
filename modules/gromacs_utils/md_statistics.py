#!/usr/bin/env python3
"""
MD Statistics Module for GROMACS Analysis

Extracts and calculates statistics from MD simulation outputs.
"""

import os
import json
import csv
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
import numpy as np

from .xvg_parser import parse_xvg, calculate_statistics


@dataclass
class MDStatistics:
    """Container for MD simulation statistics."""
    
    structure: str = ""
    simulation_type: str = "MD"
    timestep_fs: float = 2.0
    temperature_K: float = 300.0
    duration_ps: float = 0.0
    
    # Structure metrics
    rmsd: Dict[str, float] = field(default_factory=dict)
    rmsf: Dict[str, float] = field(default_factory=dict)
    gyration: Dict[str, float] = field(default_factory=dict)
    sasa: Dict[str, float] = field(default_factory=dict)
    
    # Interface metrics
    hbonds: Dict[str, float] = field(default_factory=dict)
    min_distance: Dict[str, float] = field(default_factory=dict)
    contacts: Dict[str, float] = field(default_factory=dict)
    
    # Flexible residues
    flexible_residues: List[int] = field(default_factory=list)
    
    # Raw timeseries (for plotting)
    timeseries: Dict[str, Dict] = field(default_factory=dict)


def extract_md_statistics(workdir: Path) -> MDStatistics:
    """
    Extract MD statistics from GROMACS output files.
    
    Args:
        workdir: Working directory containing XVG files
        
    Returns:
        MDStatistics object
    """
    stats = MDStatistics(structure=workdir.name)
    
    # Parse RMSD
    rmsd_data, _ = parse_xvg(workdir / 'rmsd.xvg')
    if rmsd_data:
        times = [d[0] for d in rmsd_data]
        values = [d[1] for d in rmsd_data]
        stats.rmsd = calculate_statistics(values)
        stats.timeseries['rmsd'] = {'time_ps': times, 'value_nm': values}
        stats.duration_ps = times[-1] if times else 0
    
    # Parse radius of gyration
    gyrate_data, _ = parse_xvg(workdir / 'gyrate.xvg')
    if gyrate_data:
        values = [d[1] for d in gyrate_data]
        stats.gyration = calculate_statistics(values)
        stats.timeseries['gyration'] = {
            'time_ps': [d[0] for d in gyrate_data],
            'value_nm': values
        }
    
    # Parse RMSF
    rmsf_data, _ = parse_xvg(workdir / 'rmsf.xvg')
    if rmsf_data:
        residues = [int(d[0]) for d in rmsf_data]
        values = [d[1] for d in rmsf_data]
        stats.rmsf = calculate_statistics(values)
        
        # Find flexible regions (top 10% RMSF)
        if values:
            threshold = np.percentile(values, 90)
            stats.flexible_residues = [r for r, v in zip(residues, values) if v > threshold][:20]
    
    # Parse SASA
    sasa_data, _ = parse_xvg(workdir / 'sasa.xvg')
    if sasa_data:
        values = [d[1] for d in sasa_data if len(d) > 1]
        if values:
            stats.sasa = calculate_statistics(values)
    
    # Parse H-bonds
    hbonds_data, _ = parse_xvg(workdir / 'hbonds.xvg')
    if hbonds_data:
        values = [d[1] for d in hbonds_data if len(d) > 1]
        if values:
            stats.hbonds = calculate_statistics(values)
    
    # Parse minimum distance
    mindist_data, _ = parse_xvg(workdir / 'mindist.xvg')
    if mindist_data:
        values = [d[1] for d in mindist_data if len(d) > 1]
        if values:
            stats.min_distance = calculate_statistics(values)
    
    # Parse contacts
    contacts_data, _ = parse_xvg(workdir / 'numcont.xvg')
    if contacts_data:
        values = [d[1] for d in contacts_data if len(d) > 1]
        if values:
            stats.contacts = calculate_statistics(values)
    
    return stats


def save_statistics_json(stats: MDStatistics, output_file: Path,
                         include_timeseries: bool = False) -> None:
    """
    Save MD statistics to JSON file.
    
    Args:
        stats: MDStatistics object
        output_file: Path to output file
        include_timeseries: Whether to include raw timeseries data
    """
    data = {
        'structure': stats.structure,
        'simulation': {
            'type': stats.simulation_type,
            'timestep_fs': stats.timestep_fs,
            'temperature_K': stats.temperature_K,
            'duration_ps': stats.duration_ps
        },
        'structure_metrics': {
            'rmsd_nm': stats.rmsd,
            'rmsf_nm': stats.rmsf,
            'gyration_nm': stats.gyration,
            'sasa_nm2': stats.sasa,
            'flexible_residues': stats.flexible_residues
        },
        'interface': {
            'hbonds': stats.hbonds,
            'min_distance_nm': stats.min_distance,
            'contacts': stats.contacts
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    # Save timeseries separately if requested
    if include_timeseries and stats.timeseries:
        ts_file = output_file.parent / 'timeseries_data.json'
        with open(ts_file, 'w') as f:
            json.dump(stats.timeseries, f)


def save_statistics_csv(stats: MDStatistics, output_file: Path) -> None:
    """
    Save MD statistics to CSV file.
    
    Args:
        stats: MDStatistics object
        output_file: Path to output file
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'mean', 'std', 'min', 'max', 'unit'])
        
        # Structure metrics
        for name, data, unit in [
            ('rmsd', stats.rmsd, 'nm'),
            ('rmsf', stats.rmsf, 'nm'),
            ('gyration', stats.gyration, 'nm'),
            ('sasa', stats.sasa, 'nm²')
        ]:
            if data and 'mean' in data:
                writer.writerow([
                    f'structure.{name}',
                    f"{data['mean']:.4f}",
                    f"{data['std']:.4f}",
                    f"{data['min']:.4f}",
                    f"{data['max']:.4f}",
                    unit
                ])
        
        # Interface metrics
        for name, data, unit in [
            ('hbonds', stats.hbonds, 'count'),
            ('min_distance', stats.min_distance, 'nm'),
            ('contacts', stats.contacts, 'count')
        ]:
            if data and 'mean' in data:
                writer.writerow([
                    f'interface.{name}',
                    f"{data['mean']:.4f}",
                    f"{data['std']:.4f}",
                    f"{data['min']:.4f}",
                    f"{data['max']:.4f}",
                    unit
                ])


def generate_md_statistics(workdir: Path) -> Dict[str, Path]:
    """
    Generate all MD statistics files.
    
    Args:
        workdir: Working directory
        
    Returns:
        Dictionary mapping file type to path
    """
    workdir = Path(workdir)
    stats_dir = workdir / 'statistics'
    stats_dir.mkdir(exist_ok=True)
    
    # Extract statistics
    stats = extract_md_statistics(workdir)
    
    # Save files
    files = {}
    
    json_file = stats_dir / 'md_statistics.json'
    save_statistics_json(stats, json_file, include_timeseries=True)
    files['json'] = json_file
    
    csv_file = stats_dir / 'md_statistics.csv'
    save_statistics_csv(stats, csv_file)
    files['csv'] = csv_file
    
    print(f"Created: {json_file}")
    print(f"Created: {csv_file}")
    if stats.timeseries:
        print(f"Created: {stats_dir / 'timeseries_data.json'}")
        files['timeseries'] = stats_dir / 'timeseries_data.json'
    
    return files


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate MD statistics")
    parser.add_argument("--workdir", required=True, help="Working directory")
    args = parser.parse_args()
    
    files = generate_md_statistics(Path(args.workdir))
    
    print("\nStatistics files generated:")
    for name, path in files.items():
        print(f"  {name}: {path}")
