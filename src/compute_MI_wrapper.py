"""
MIIC Helper - Python Wrapper for R MIIC Feature Selection

This module provides a Python interface to run MIIC (Multivariate Information-based 
Inductive Causation) feature selection using an R backend. It facilitates the selection
of relevant features from genomic data by computing mutual information scores.

Author: Ali Chemkhi
Date: June 2025
"""

import os
import subprocess

def run_miic_selection(count_matrix_path, metadata_path, variables_of_interest, selection_pool, output_path, r_script_path=None):
    """
    Python wrapper function to call the R MIIC feature selection script.
    
    Parameters:
    -----------
    count_matrix_path : str
        Path to the count matrix CSV file
    metadata_path : str
        Path to the metadata CSV file
    variables_of_interest : list or str
        List of variables of interest (genes or metadata columns), or comma-separated string
    selection_pool : list or str
        List of genes for selection pool, or comma-separated string
    output_path : str
        Path for the output CSV file
    r_script_path : str, optional
        Path to the R script. If None, uses default path
    
    Returns:
    --------
    str
        Path to the output file containing mutual information scores
    
    Raises:
    -------
    subprocess.CalledProcessError
        If the R script execution fails
    FileNotFoundError
        If input files or R script are not found
    """
    variables_of_interest = list(variables_of_interest)
    selection_pool = list(selection_pool)
    
    # Default R script path
    if r_script_path is None:
        r_script_path = '/Users/alichemkhi/Desktop/myProjects/miic_helper/src/run_miic_select.R'
    
    # Convert lists to comma-separated strings if needed
    if isinstance(variables_of_interest, list):
        variables_str = ','.join(variables_of_interest)
    else:
        variables_str = str(variables_of_interest)
    
    if isinstance(selection_pool, list):
        selection_pool_str = ','.join(selection_pool)
    else:
        selection_pool_str = str(selection_pool)
    
    # Validate input files exist
    if not os.path.exists(count_matrix_path):
        raise FileNotFoundError(f"Count matrix file not found: {count_matrix_path}")
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"R script not found: {r_script_path}")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Build the R script command
    cmd = [
        'Rscript',
        r_script_path,
        '--matrix', count_matrix_path,
        '--metadata', metadata_path,
        '--variables_of_interest', variables_str,
        '--selection_pool', selection_pool_str,
        '--output', output_path
    ]


    
    print(f"Running MIIC selection with command:")
    print(' '.join(cmd))
    print(f"Variables of interest: {variables_str}")
    print(f"Selection pool: {selection_pool_str}")
    print(f"Output will be saved to: {output_path}")
    
    try:
        # Run the R script and wait for it to finish
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Print R script output
        if result.stdout:
            print("\nR script output:")
            print(result.stdout)

        print(f"\nMIIC selection completed successfully!")
        print(f"Results saved to: {output_path}")
        
        return output_path
        
    except subprocess.CalledProcessError as e:
        print(f"\nError running R script:")
        print(f"Return code: {e.returncode}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        raise
    except Exception as e:
        print(f"\nUnexpected error: {str(e)}")
        raise