"""
bioconfigme - Configuration management for meta-analysis pipeline

Provides centralized access to YAML configuration files for analysis
and software settings. All Snakemake rules and scripts should use this
module instead of reading YAML files directly.

Usage in Snakemake rules:
    import sys
    sys.path.append("utils")
    from bioconfigme import get_results_dir, get_software_module, get_analysis_value
"""

import yaml
import os
import sys
from pathlib import Path
from typing import Any, Optional, Sequence, Dict, List, Union


# Global config cache
_CONFIGS = {
    'analysis': None,
    'software': None
}


def _load_config(config_type: str) -> Dict:
    """
    Load configuration file from configs/ directory.
    
    Args:
        config_type: Either 'analysis' or 'software'
    
    Returns:
        Parsed YAML configuration dictionary
    
    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If YAML is malformed
    """
    if _CONFIGS[config_type] is not None:
        return _CONFIGS[config_type]
    
    config_file = Path(f"configs/{config_type}.yml")
    
    if not config_file.exists():
        raise FileNotFoundError(
            f"Configuration file not found: {config_file}\n"
            f"Expected location: configs/{config_type}.yml"
        )
    
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
            _CONFIGS[config_type] = config
            return config
    except yaml.YAMLError as e:
        raise ValueError(f"Failed to parse {config_file}: {e}")


def get_results_dir() -> str:
    """
    Get the absolute path to the results directory.
    
    Returns:
        Absolute path to results directory from analysis.yml
    
    Raises:
        ValueError: If results_dir is not defined in config
    """
    config = _load_config('analysis')
    
    if 'results_dir' not in config:
        raise ValueError(
            "Missing required key 'results_dir' in configs/analysis.yml"
        )
    
    results_dir = config['results_dir']
    
    if not results_dir or results_dir == "/path/to/results":
        raise ValueError(
            "results_dir in configs/analysis.yml must be set to an absolute path.\n"
            "Current value is placeholder: /path/to/results"
        )
    
    return str(results_dir)


def get_analysis_value(
    path: Sequence[str],
    *,
    required: bool = True,
    default: Optional[Any] = None
) -> Any:
    """
    Get a value from analysis configuration using nested key path.
    
    Args:
        path: Sequence of keys to navigate (e.g., ['metal', 'scheme'])
        required: If True, raise error when key is missing
        default: Default value if key is missing and required=False
    
    Returns:
        Configuration value at the specified path
    
    Raises:
        ValueError: If required=True and key path doesn't exist
    
    Examples:
        >>> get_analysis_value(['target_build'])
        'hg38'
        >>> get_analysis_value(['metal', 'scheme'])
        'STDERR'
        >>> get_analysis_value(['datasets', 'CHI', 'build'])
        'hg19'
    """
    config = _load_config('analysis')
    
    current = config
    for key in path:
        if not isinstance(current, dict) or key not in current:
            if required:
                path_str = ' -> '.join(path)
                raise ValueError(
                    f"Required configuration key not found: {path_str}\n"
                    f"Missing key: {key}\n"
                    f"Check configs/analysis.yml"
                )
            return default
        current = current[key]
    
    return current


def get_software_module(tool_name: str) -> str:
    """
    Get the module load command or path for a software tool.
    
    Checks for 'module' key first, then 'path' key in software.yml.
    Returns the appropriate string for loading the tool.
    
    Args:
        tool_name: Name of the tool (e.g., 'metal', 'python', 'liftover')
    
    Returns:
        Module name or absolute path to the tool
    
    Raises:
        ValueError: If tool is not found or has neither module nor path
    
    Examples:
        >>> get_software_module('metal')
        'metal/2020-05-05'  # or '/usr/local/bin/metal'
    """
    config = _load_config('software')
    
    if tool_name not in config:
        raise ValueError(
            f"Tool '{tool_name}' not found in configs/software.yml\n"
            f"Available tools: {', '.join(config.keys())}"
        )
    
    tool_config = config[tool_name]
    
    # Check for module first, then path
    if 'module' in tool_config and tool_config['module']:
        return tool_config['module']
    elif 'path' in tool_config and tool_config['path']:
        return tool_config['path']
    else:
        raise ValueError(
            f"Tool '{tool_name}' must have either 'module' or 'path' defined in configs/software.yml"
        )


def get_dataset_list() -> List[str]:
    """
    Get list of all dataset names configured for meta-analysis.
    
    Returns:
        List of dataset names (e.g., ['CHI', 'EUR', 'MEX'])
    """
    config = _load_config('analysis')
    
    if 'datasets' not in config:
        raise ValueError("No 'datasets' section found in configs/analysis.yml")
    
    return list(config['datasets'].keys())


def get_dataset_config(dataset_name: str) -> Dict:
    """
    Get complete configuration for a specific dataset.
    
    Args:
        dataset_name: Name of the dataset (e.g., 'CHI', 'EUR')
    
    Returns:
        Dictionary with dataset configuration including file, build, columns
    
    Raises:
        ValueError: If dataset not found in configuration
    """
    datasets = get_analysis_value(['datasets'])
    
    if dataset_name not in datasets:
        available = ', '.join(datasets.keys())
        raise ValueError(
            f"Dataset '{dataset_name}' not found in configs/analysis.yml\n"
            f"Available datasets: {available}"
        )
    
    return datasets[dataset_name]


def get_target_build() -> str:
    """
    Get the target genome build for the meta-analysis.
    
    Returns:
        Target build (e.g., 'hg38', 'hg19')
    """
    return get_analysis_value(['target_build'])


def get_liftover_chain(from_build: str, to_build: str) -> str:
    """
    Get the path to liftOver chain file for build conversion.
    
    Args:
        from_build: Source genome build (e.g., 'hg19')
        to_build: Target genome build (e.g., 'hg38')
    
    Returns:
        Path to chain file
    
    Raises:
        ValueError: If chain file path not configured
    """
    chain_key = f"{from_build}_to_{to_build}"
    
    try:
        chain_path = get_analysis_value(
            ['resources', 'liftover_chains', chain_key],
            required=True
        )
        return chain_path
    except ValueError:
        raise ValueError(
            f"No liftOver chain file configured for {from_build} -> {to_build}\n"
            f"Expected key: resources.liftover_chains.{chain_key}\n"
            f"Check configs/analysis.yml"
        )


def get_qc_thresholds() -> Dict:
    """
    Get quality control thresholds for filtering variants.
    
    Returns:
        Dictionary with QC parameters (min_info, min_maf, etc.)
    """
    return get_analysis_value(['qc'], default={})


def get_liftover_binary() -> str:
    """
    Get the liftOver binary name or path.
    
    Returns liftOver binary name from software.yml if specified,
    otherwise defaults to 'liftOver' (assumes it's in PATH after module load).
    
    Returns:
        Binary name or path for liftOver executable
    
    Examples:
        >>> get_liftover_binary()
        'liftOver'  # or '/usr/local/bin/liftOver' if path specified
    """
    config = _load_config('software')
    
    if 'liftover' not in config:
        return 'liftOver'  # Default binary name
    
    tool_config = config['liftover']
    
    # If path is specified, return it
    if 'path' in tool_config and tool_config['path']:
        return tool_config['path']
    
    # If binary name is specified, return it
    if 'binary' in tool_config and tool_config['binary']:
        return tool_config['binary']
    
    # Default to 'liftOver'
    return 'liftOver'


def get_metal_binary() -> str:
    """
    Get the METAL binary name or path.
    
    Returns METAL binary name from software.yml if specified,
    otherwise defaults to 'metal' (assumes it's in PATH after module load).
    
    Returns:
        Binary name or path for METAL executable
    
    Examples:
        >>> get_metal_binary()
        'metal'  # or '/usr/local/bin/metal' if path specified
    """
    config = _load_config('software')
    
    if 'metal' not in config:
        return 'metal'  # Default binary name
    
    tool_config = config['metal']
    
    # If path is specified, return it
    if 'path' in tool_config and tool_config['path']:
        return tool_config['path']
    
    # If binary name is specified, return it
    if 'binary' in tool_config and tool_config['binary']:
        return tool_config['binary']
    
    # Default to 'metal'
    return 'metal'


def get_metal_config() -> Dict:
    """
    Get METAL meta-analysis configuration.
    
    Returns:
        Dictionary with METAL parameters (scheme, genomiccontrol, etc.)
    """
    return get_analysis_value(['metal'], default={})


def get_metal_combinations() -> Dict[str, List[str]]:
    """
    Get METAL study combinations for meta-analysis.
    
    Returns dictionary mapping combination names to lists of dataset names.
    Each combination will be analyzed separately by METAL.
    
    Returns:
        Dictionary with combination_name: [dataset_list]
    
    Example:
        >>> get_metal_combinations()
        {'hisp_mvp': ['MEX', 'MEX123', 'LAMR', 'CLM', 'MVP']}
    
    Raises:
        ValueError: If study_combinations not defined or invalid
    """
    combinations = get_analysis_value(['metal', 'study_combinations'], default={})
    
    if not combinations:
        raise ValueError(
            "No study combinations defined in configs/analysis.yml\n"
            "Add 'metal.study_combinations' with combination_name: [dataset_list]"
        )
    
    # Validate that all datasets in combinations exist
    all_datasets = get_dataset_list()
    for combo_name, datasets in combinations.items():
        if not isinstance(datasets, list):
            raise ValueError(
                f"Study combination '{combo_name}' must be a list of dataset names"
            )
        
        invalid_datasets = [d for d in datasets if d not in all_datasets]
        if invalid_datasets:
            raise ValueError(
                f"Study combination '{combo_name}' contains invalid datasets: {invalid_datasets}\n"
                f"Valid datasets: {all_datasets}"
            )
    
    return combinations


def get_standardized_columns() -> List[str]:
    """
    Get the list of standardized column names used across the pipeline.
    
    These are the target column names after Step 0 (standardization).
    All datasets will be transformed to use these column names.
    Note: 'or' may be added dynamically if beta is not available.
    
    Returns:
        List of standardized column names in order
    
    Example:
        >>> get_standardized_columns()
        ['snpid', 'chrom', 'pos', 'ea', 'nea', 'eaf', 'neaf', 'beta', 'se', 'p', 'n']
    """
    return ['snpid', 'chrom', 'pos', 'ea', 'nea', 'eaf', 'neaf', 'beta', 'se', 'p', 'n']


def get_required_columns() -> List[str]:
    """
    Get the minimal required columns for meta-analysis.
    
    These columns MUST be present after standardization and calculation steps.
    
    Returns:
        List of required column names
    
    Example:
        >>> get_required_columns()
        ['snpid', 'chrom', 'pos', 'ea', 'nea', 'beta', 'se', 'p']
    """
    return ['snpid', 'chrom', 'pos', 'ea', 'nea', 'beta', 'se', 'p']


def get_column_mapping(dataset_name: str) -> Dict[str, str]:
    """
    Get the column name mapping for a specific dataset.
    
    Maps standardized names (keys) to original column names (values).
    
    Args:
        dataset_name: Name of the dataset (e.g., 'CHI', 'EUR')
    
    Returns:
        Dictionary mapping standardized_name -> original_name
    
    Example:
        >>> get_column_mapping('CHI')
        {'snpid': 'rsid', 'chrom': 'chrom', 'pos': 'pos', ...}
    """
    dataset_config = get_dataset_config(dataset_name)
    
    if 'columns' not in dataset_config:
        raise ValueError(
            f"No column mapping found for dataset '{dataset_name}'\n"
            f"Check configs/analysis.yml datasets.{dataset_name}.columns"
        )
    
    return dataset_config['columns']


def has_variant_id_column(dataset_name: str) -> bool:
    """
    Check if a dataset has variant_id column (in chrom:pos format).
    
    Files with variant_id need to be parsed to extract chrom and pos.
    
    Args:
        dataset_name: Name of the dataset (e.g., 'CLM', 'MEX')
    
    Returns:
        True if dataset has variant_id column, False otherwise
    
    Example:
        >>> has_variant_id_column('CLM')
        True
        >>> has_variant_id_column('CHI')
        False
    """
    col_mapping = get_column_mapping(dataset_name)
    return 'variant_id' in col_mapping


def has_separate_chrom_pos(dataset_name: str) -> bool:
    """
    Check if a dataset has separate chrom and pos columns.
    
    Args:
        dataset_name: Name of the dataset (e.g., 'CHI', 'EUR')
    
    Returns:
        True if dataset has both chrom and pos columns, False otherwise
    """
    col_mapping = get_column_mapping(dataset_name)
    return 'chrom' in col_mapping and 'pos' in col_mapping


# Helper function for scripts that need to validate config on import
def validate_config():
    """
    Validate that essential configuration keys are present.
    Call this at the start of pipeline execution to fail fast.
    
    Raises:
        ValueError: If any essential configuration is missing or invalid
    """
    errors = []
    
    # Check results_dir
    try:
        results_dir = get_results_dir()
        if not os.path.isabs(results_dir):
            errors.append("results_dir must be an absolute path")
    except Exception as e:
        errors.append(f"results_dir error: {e}")
    
    # Check target_build
    try:
        target_build = get_target_build()
        if target_build not in ['hg19', 'hg38', 'hg18']:
            errors.append(f"Invalid target_build: {target_build}")
    except Exception as e:
        errors.append(f"target_build error: {e}")
    
    # Check datasets exist
    try:
        datasets = get_dataset_list()
        if not datasets:
            errors.append("No datasets configured")
    except Exception as e:
        errors.append(f"datasets error: {e}")
    
    if errors:
        raise ValueError(
            "Configuration validation failed:\n" + 
            "\n".join(f"  - {err}" for err in errors)
        )


if __name__ == "__main__":
    # Basic validation when run directly
    print("Validating bioconfigme configuration...")
    try:
        validate_config()
        print("✓ Configuration valid")
        print(f"✓ Results directory: {get_results_dir()}")
        print(f"✓ Target build: {get_target_build()}")
        print(f"✓ Datasets: {', '.join(get_dataset_list())}")
    except Exception as e:
        print(f"✗ Configuration error: {e}", file=sys.stderr)
        sys.exit(1)
