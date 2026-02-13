#!/usr/bin/env python3
"""
Standardize column names and extract required columns from GWAS summary statistics.

This script:
1. Reads column mapping from bioconfigme configuration
2. Parses variant_id (chrom:pos format) if present
3. Extracts and renames columns to standardized names
4. Creates n column from n_cases + n_controls if missing
5. Reports missing columns without failing
6. Outputs standardized TSV file

Usage:
    python standardize_columns.py --dataset CHI --output results/00_standardized/CHI.tsv --done results/00_standardized/CHI.done --log results/log/CHI.log
"""

import sys
import argparse
import pandas as pd
from pathlib import Path

# Add utils to path for bioconfigme import
sys.path.append("utils")
from bioconfigme import (
    get_dataset_config,
    get_column_mapping,
    has_variant_id_column,
    has_separate_chrom_pos,
    get_standardized_columns
)


def parse_variant_id(df, variant_id_col):
    """
    Parse variant_id column in 'chrom:pos' format into separate chrom and pos columns.
    Handles formats like "chrom:pos" or "chrom:pos:ref:alt".
    
    Args:
        df: DataFrame with variant_id column
        variant_id_col: Name of variant_id column
    
    Returns:
        DataFrame with added chrom and pos columns
    """
    print(f"Parsing variant_id column '{variant_id_col}' (format: chrom:pos)")
    
    # Split variant_id on ':' - extract first two fields (chrom and pos)
    split_data = df[variant_id_col].astype(str).str.split(':', expand=True)
    
    if split_data.shape[1] < 2:
        raise ValueError(f"variant_id column '{variant_id_col}' does not contain ':' separator")
    
    # First field is chrom, second field is pos
    df['chrom'] = split_data[0]
    df['pos'] = split_data[1].astype(int)
    
    print(f"  Extracted chrom and pos from {len(df)} variants")
    return df


def is_column_all_na(series):
    """
    Check if a pandas Series contains all NA/null/empty values.
    
    Args:
        series: pandas Series to check
    
    Returns:
        True if all values are NA/null/empty, False otherwise
    """
    if series is None or len(series) == 0:
        return True
    
    # Check for various types of missing/empty values
    is_na = series.isna()
    is_empty_string = series.astype(str).str.strip().isin(['', 'NA', 'na', 'NaN', 'nan'])
    
    # Consider column all NA if all values are NA or empty
    all_missing = (is_na | is_empty_string).all()
    
    return all_missing


def extract_column_with_split(df, column_spec):
    """
    Extract a column and optionally split it to get specific element.
    
    Handles column specifications like:
      - "MarkerName" -> just extract the column as-is
      - "MarkerName:3" -> split MarkerName by ':' and take 3rd element (1-based indexing)
      - "variant_id:2" -> split variant_id by ':' and take 2nd element
    
    Args:
        df: DataFrame containing the column
        column_spec: Column specification (e.g., "MarkerName" or "MarkerName:3")
    
    Returns:
        tuple: (pandas Series with extracted data, description string for logging)
    """
    if ':' not in column_spec:
        # Simple column extraction - no splitting
        if column_spec in df.columns:
            return df[column_spec], f"from '{column_spec}'"
        else:
            return None, f"'{column_spec}' NOT FOUND"
    
    # Parse column_name:index format
    parts = column_spec.split(':')
    if len(parts) != 2:
        raise ValueError(f"Invalid column specification: '{column_spec}'. Expected 'column_name' or 'column_name:index'")
    
    col_name, index_str = parts
    
    if col_name not in df.columns:
        return None, f"'{col_name}' NOT FOUND"
    
    try:
        # Convert to 0-based index (user specifies 1-based)
        split_index = int(index_str) - 1
        
        # Split the column values by ':' and extract the specified element
        split_data = df[col_name].astype(str).str.split(':', expand=True)
        
        if split_index < 0 or split_index >= split_data.shape[1]:
            raise ValueError(f"Index {index_str} out of range for column '{col_name}' (has {split_data.shape[1]} parts)")
        
        extracted_series = split_data[split_index]
        description = f"from '{col_name}' (split by ':', element {index_str})"
        
        return extracted_series, description
        
    except ValueError as e:
        raise ValueError(f"Invalid index in column specification '{column_spec}': {e}")


def standardize_dataset(dataset_name, output_file, done_file, log_file):
    """
    Standardize column names for a single dataset.
    
    Args:
        dataset_name: Name of dataset (e.g., 'CHI', 'EUR')
        output_file: Path to output TSV file
        done_file: Path to .done marker file
        log_file: Path to log file
    """
    print(f"\n{'='*60}")
    print(f"Standardizing dataset: {dataset_name}")
    print(f"{'='*60}\n")
    
    # Get dataset configuration
    dataset_config = get_dataset_config(dataset_name)
    input_file = dataset_config['file']
    n_cases = dataset_config.get('n_cases', 0)
    n_controls = dataset_config.get('n_controls', 0)
    
    print(f"Input file: {input_file}")
    print(f"Sample size: {n_cases} cases + {n_controls} controls = {n_cases + n_controls} total")
    
    # Get column mapping
    col_mapping = get_column_mapping(dataset_name)
    print(f"\nColumn mapping (standardized -> original):")
    for std_name, orig_name in col_mapping.items():
        # Highlight if this uses column splitting
        if ':' in orig_name and not orig_name.startswith('variant_id'):
            print(f"  {std_name:15s} <- {orig_name} [SPLIT MODE]")
        else:
            print(f"  {std_name:15s} <- {orig_name}")
    
    # Read input file
    print(f"\nReading input file...")
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
        print(f"  Loaded {len(df):,} variants with {len(df.columns)} columns")
    except Exception as e:
        raise RuntimeError(f"Failed to read input file: {e}")
    
    # Parse variant_id if present
    if has_variant_id_column(dataset_name) and not has_separate_chrom_pos(dataset_name):
        variant_id_col = col_mapping.get('variant_id')
        if variant_id_col and variant_id_col in df.columns:
            df = parse_variant_id(df, variant_id_col)
            print(f"  Added chrom and pos columns from variant_id")
    
    # Get standardized column names
    std_columns = get_standardized_columns()
    
    # Check for beta and or columns
    has_beta = 'beta' in col_mapping and col_mapping.get('beta') in df.columns
    has_or = 'or' in col_mapping and col_mapping.get('or') in df.columns
    
    # Add 'or' to standard columns if beta is missing but or is available
    if not has_beta and has_or and 'or' not in std_columns:
        std_columns = std_columns.copy()
        # Insert 'or' after 'neaf' (before 'beta' position)
        beta_idx = std_columns.index('beta') if 'beta' in std_columns else len(std_columns)
        std_columns.insert(beta_idx, 'or')
        print(f"\n  NOTE: 'beta' not available but 'or' is - including 'or' column for later conversion")
    
    # Map columns and track missing
    print(f"\nExtracting and renaming columns...")
    mapped_data = {}
    missing_columns = []
    all_na_columns = []
    
    for std_name in std_columns:
        # First, check if this column was created from variant_id parsing
        if std_name in ['chrom', 'pos'] and std_name in df.columns:
            # Check if column is all NA
            if is_column_all_na(df[std_name]):
                all_na_columns.append(f"{std_name} (from variant_id)")
                print(f"  ⊘ {std_name:10s} (from variant_id) [ALL NA - skipped]")
            else:
                mapped_data[std_name] = df[std_name]
                print(f"  ✓ {std_name:10s} (from variant_id parsing)")
        elif std_name in col_mapping:
            # Get the column specification (may include split instruction)
            column_spec = col_mapping[std_name]
            
            # Extract column (with optional splitting)
            extracted_series, description = extract_column_with_split(df, column_spec)
            
            if extracted_series is not None:
                # Check if column is all NA before adding
                if is_column_all_na(extracted_series):
                    all_na_columns.append(f"{std_name} ({description})")
                    print(f"  ⊘ {std_name:10s} {description} [ALL NA - skipped]")
                else:
                    mapped_data[std_name] = extracted_series
                    print(f"  ✓ {std_name:10s} {description}")
            else:
                # Column specified in config but not in file
                missing_columns.append(f"{std_name} ({description})")
                print(f"  ✗ {std_name:10s} {description} [MISSING - will skip]")
        else:
            # Column not in mapping - check if it's an optional column
            if std_name in ['eaf', 'neaf', 'n', 'or']:
                missing_columns.append(f"{std_name} (optional)")
                print(f"  - {std_name:10s} [NOT MAPPED - optional]")
            else:
                missing_columns.append(f"{std_name} (required)")
                print(f"  ✗ {std_name:10s} [NOT MAPPED - required]")
    
    # Create n column if missing but we have case/control counts
    if 'n' not in mapped_data and (n_cases > 0 or n_controls > 0):
        total_n = n_cases + n_controls
        mapped_data['n'] = total_n
        print(f"\n  ✓ Created 'n' column: {total_n} (cases: {n_cases} + controls: {n_controls})")
    
    # Create standardized dataframe
    std_df = pd.DataFrame(mapped_data)
    
    # Report results
    if all_na_columns:
        print(f"\n  Columns skipped (all NA) ({len(all_na_columns)}):")
        for col in all_na_columns:
            print(f"    - {col}")
    
    # Verify required columns are present
    required_cols = ['snpid', 'chrom', 'pos', 'ea', 'nea', 'se', 'p']
    missing_required = [col for col in required_cols if col not in std_df.columns]
    
    if missing_required:
        print(f"\n  WARNING: Missing required columns: {', '.join(missing_required)}")
        print(f"  Dataset may not be suitable for meta-analysis without these columns.")
    
    # Check for beta or OR
    has_beta_in_output = 'beta' in std_df.columns
    has_or_in_output = 'or' in std_df.columns
    
    if not has_beta_in_output and not has_or_in_output:
        print(f"\n  WARNING: Neither 'beta' nor 'or' column found in output!")
        print(f"  Effect size calculation may fail in downstream steps.")
    elif not has_beta_in_output and has_or_in_output:
        print(f"\n  NOTE: Dataset has 'or' but not 'beta' - BETA will be calculated in Step 3")
        print(f"        (beta = log(or))s.")
    
    # Check for beta or OR
    has_beta_in_output = 'beta' in std_df.columns
    has_or_in_output = 'or' in std_df.columns
    
    if not has_beta_in_output and not has_or_in_output:
        print(f"\n  WARNING: Neither 'beta' nor 'or' column found in output!")
        print(f"  Effect size calculation may fail in downstream steps.")
    elif not has_beta_in_output and has_or_in_output:
        print(f"\n  NOTE: Dataset has 'or' but not 'beta' - BETA will be calculated in Step 3")
        print(f"        (beta = log(or))")
    
    # Write output
    print(f"\nWriting standardized output...")
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Check for split columns and report
    split_columns_used = [(std, orig) for std, orig in col_mapping.items() 
                          if ':' in orig and std in std_df.columns and not orig.startswith('variant_id')]
    
    if split_columns_used:
        print(f"\n  Columns extracted using split mode:")
        for std_col, orig_spec in split_columns_used:
            print(f"    - {std_col}: {orig_spec}")
    
    # Write TSV file
    std_df.to_csv(output_file, sep='\t', index=False, na_rep='NA')
    print(f"  Wrote {len(std_df):,} variants to: {output_file}")
    
    # Write .done marker
    split_columns_info = [f"{std}: {orig}" for std, orig in col_mapping.items() if ':' in orig and std in std_df.columns]
    
    with open(done_file, 'w') as f:
        f.write(f"Standardization completed for {dataset_name}\n")
        f.write(f"Input: {input_file}\n")
        f.write(f"Output: {output_file}\n")
        f.write(f"Variants: {len(std_df):,}\n")
        f.write(f"Columns: {', '.join(std_df.columns)}\n")
        if split_columns_info:
            f.write(f"Split columns: {', '.join(split_columns_info)}\n")
        if missing_columns:
            f.write(f"Missing: {', '.join(missing_columns)}\n")
        if all_na_columns:
            f.write(f"Skipped (all NA): {', '.join(all_na_columns)}\n")
    
    print(f"  Created done marker: {done_file}")
    print(f"\n{'='*60}")
    print(f"Standardization completed successfully for {dataset_name}")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Standardize GWAS summary statistics column names",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python standardize_columns.py --dataset CHI --output results/00_standardized/CHI.tsv --done results/00_standardized/CHI.done --log results/log/CHI.log
  python standardize_columns.py --dataset EUR --output results/00_standardized/EUR.tsv --done results/00_standardized/EUR.done --log results/log/EUR.log

Required columns (extracted from config):
  snpid, chrom, pos, ea, nea, beta/or, se, p
  
Optional columns:
  eaf, neaf, n (created from n_cases + n_controls if missing)
"""
    )
    
    parser.add_argument('--dataset', required=True,
                        help='Dataset name (e.g., CHI, EUR, MVP)')
    parser.add_argument('--output', required=True,
                        help='Output TSV file path')
    parser.add_argument('--done', required=True,
                        help='Done marker file path')
    parser.add_argument('--log', required=True,
                        help='Log file path')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.dataset:
        parser.error("--dataset is required")
    
    try:
        standardize_dataset(args.dataset, args.output, args.done, args.log)
        sys.exit(0)
    except Exception as e:
        print(f"\nERROR: Standardization failed for {args.dataset}", file=sys.stderr)
        print(f"  {type(e).__name__}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
