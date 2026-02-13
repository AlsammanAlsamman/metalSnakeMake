#!/usr/bin/env python3
"""
Run METAL meta-analysis for a study combination.
Generates METAL configuration file and executes METAL.
"""

import sys
import os
import argparse
import subprocess
from pathlib import Path


def create_metal_config(
    combination_name,
    datasets,
    input_dir,
    output_dir,
    scheme,
    log_file
):
    """
    Create METAL configuration file for a study combination.
    
    Args:
        combination_name: Name of the study combination
        datasets: List of dataset names to meta-analyze
        input_dir: Directory containing input TSV files
        output_dir: Directory for METAL output
        scheme: METAL scheme (STDERR or SAMPLESIZE)
        log_file: Path to log file
    
    Returns:
        Path to generated METAL config file
    """
    config_file = Path(output_dir) / f"{combination_name}.metal"
    
    print(f"\n{'='*60}")
    print(f"Creating METAL configuration: {combination_name}")
    print(f"{'='*60}\n")
    print(f"Configuration file: {config_file}")
    print(f"Study combination: {combination_name}")
    print(f"Datasets ({len(datasets)}): {', '.join(datasets)}\n")
    
    with open(config_file, 'w') as f:
        # Write header
        f.write(f"# METAL analysis configuration - {combination_name}\n")
        f.write(f"SCHEME {scheme}\n")
        f.write(f"# GENOMICCONTROL ON\n")
        f.write(f"SEPARATOR TAB\n\n")
        
        # Write dataset blocks
        for dataset in datasets:
            input_file = Path(input_dir) / f"{dataset}.tsv"
            
            if not input_file.exists():
                print(f"  WARNING: Input file not found: {input_file}")
                print(f"           Skipping dataset: {dataset}\n")
                continue
            
            print(f"  Adding dataset: {dataset}")
            print(f"    Input: {input_file}")
            
            f.write(f"# Process {dataset} dataset\n")
            f.write(f"MARKER   markername\n")
            f.write(f"ALLELE   nea ea\n")
            f.write(f"FREQ     eaf\n")
            f.write(f"EFFECT   beta\n")
            f.write(f"STDERR   se\n")
            f.write(f"PVAL     p\n")
            f.write(f"WEIGHT   n\n")
            f.write(f"\n")
            f.write(f"PROCESS {input_file}\n\n")
        
        # Write output and analysis commands
        output_prefix = Path(output_dir) / combination_name
        f.write(f"OUTFILE {output_prefix} .tbl\n")
        f.write(f"ANALYZE HETEROGENEITY\n")
        f.write(f"QUIT\n")
    
    print(f"\n  Created METAL config: {config_file}")
    return config_file


def run_metal(metal_binary, config_file, log_file):
    """
    Execute METAL with the generated configuration file.
    
    Args:
        metal_binary: Path to METAL binary
        config_file: Path to METAL configuration file
        log_file: Path to log file
    
    Returns:
        Exit code from METAL
    """
    print(f"\n{'='*60}")
    print(f"Running METAL meta-analysis")
    print(f"{'='*60}\n")
    print(f"METAL binary: {metal_binary}")
    print(f"Config file: {config_file}\n")
    
    # Run METAL and capture output
    cmd = [metal_binary, config_file]
    
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False
        )
        
        # Write output to log
        with open(log_file, 'a') as f:
            f.write(result.stdout)
        
        # Also print to console
        print(result.stdout)
        
        if result.returncode != 0:
            print(f"\nERROR: METAL failed with exit code {result.returncode}")
            return result.returncode
        
        print(f"\nMETAL completed successfully")
        return 0
        
    except Exception as e:
        error_msg = f"ERROR: Failed to run METAL: {e}"
        print(f"\n{error_msg}")
        with open(log_file, 'a') as f:
            f.write(f"\n{error_msg}\n")
        return 1


def main():
    parser = argparse.ArgumentParser(
        description="Run METAL meta-analysis for a study combination",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--combination', required=True,
                        help='Study combination name')
    parser.add_argument('--datasets', required=True, nargs='+',
                        help='List of dataset names to meta-analyze')
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing input TSV files')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for METAL results')
    parser.add_argument('--done', required=True,
                        help='Done marker file path')
    parser.add_argument('--log', required=True,
                        help='Log file path')
    parser.add_argument('--scheme', default='STDERR',
                        help='METAL scheme (STDERR or SAMPLESIZE)')
    parser.add_argument('--metal_bin', required=True,
                        help='Path to METAL binary')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Initialize log file
    with open(args.log, 'w') as f:
        f.write(f"{'='*60}\n")
        f.write(f"METAL Meta-Analysis: {args.combination}\n")
        f.write(f"{'='*60}\n\n")
        f.write(f"Combination: {args.combination}\n")
        f.write(f"Datasets: {', '.join(args.datasets)}\n")
        f.write(f"Scheme: {args.scheme}\n")
        f.write(f"Input directory: {args.input_dir}\n")
        f.write(f"Output directory: {args.output_dir}\n\n")
    
    try:
        # Create METAL configuration
        config_file = create_metal_config(
            combination_name=args.combination,
            datasets=args.datasets,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            scheme=args.scheme,
            log_file=args.log
        )
        
        # Run METAL
        exit_code = run_metal(
            metal_binary=args.metal_bin,
            config_file=config_file,
            log_file=args.log
        )
        
        if exit_code != 0:
            sys.exit(exit_code)
        
        # Check output files
        output_file = Path(args.output_dir) / f"{args.combination}1.tbl"
        if not output_file.exists():
            print(f"\nERROR: Expected output file not found: {output_file}")
            sys.exit(1)
        
        print(f"\nOutput files created:")
        for f in Path(args.output_dir).glob(f"{args.combination}*"):
            print(f"  {f}")
        
        # Write .done marker
        with open(args.done, 'w') as f:
            f.write(f"METAL meta-analysis completed for {args.combination}\n")
            f.write(f"Combination: {args.combination}\n")
            f.write(f"Datasets: {', '.join(args.datasets)}\n")
            f.write(f"Output directory: {args.output_dir}\n")
            f.write(f"Config file: {config_file}\n")
        
        print(f"\nCreated done marker: {args.done}")
        print(f"\n{'='*60}")
        print(f"METAL meta-analysis completed successfully")
        print(f"{'='*60}\n")
        
    except Exception as e:
        error_msg = f"\nERROR: {type(e).__name__}: {e}"
        print(error_msg)
        with open(args.log, 'a') as f:
            f.write(f"{error_msg}\n")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
