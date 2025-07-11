"""
MIIC Helper - Python Wrapper for R MIIC Feature Selection

This module provides a Python interface to run MIIC (Multivariate Information-based 
Inductive Causation) feature selection using an R backend. 

Author: Ali Chemkhi
Date: June 2025
"""

import os
import subprocess
# import sys
# sys.path.append('/Users/alichemkhi/Desktop/myProjects/miic_helper/src') # Adjust this path as needed
from process_anndata_object import save_anndata_to_csv 


def run_miic_selection(adata, variables_of_interest, selection_pool, output_path, compute_meta_mi_scores=False, r_script_path=None, run_script=True):
    """
    Python wrapper function to call the R MIIC feature selection script.
    
    This function takes an AnnData object, saves it as tmp CSV files, and calls an R script
    to perform MIIC feature selection. 
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object containing the count matrix and metadata
    variables_of_interest : list or str
        List of variables of interest (genes or metadata columns), or comma-separated string.
        Cannot be empty.
    selection_pool : list, str, or None
        List of genes for selection pool, comma-separated string, or None for no selection pool
    output_path : str
        Path for the output CSV file
    compute_meta_mi_scores : bool, optional
        Whether to compute mutual information scores for metadata variables (default: False)
    r_script_path : str, optional
        Path to the R script. If None, uses default path
    run_script : bool, optional
        If True, runs the R script; if False, just prints the command without executing (default: True)
    
    Returns:
    --------
    str
        Path to the output file containing mutual information scores
    
    Raises:
    -------
    ValueError
        If variables_of_interest is empty
    subprocess.CalledProcessError
        If the R script execution fails
    FileNotFoundError
        If input files or R script are not found
    """

    # Save the AnnData object to tmp CSV files so that the R script can read them
    # Note: This assumes that the adata object has been preprocessed and contains the necessary data
    count_matrix_path, metadata_path = save_anndata_to_csv(
        adata, 
        f"/tmp/tmp_annadata_raw_matrix",
        save_counts=True,
        save_metadata=True,
        transpose_counts=False
    )


    # Validate inputs
    if not variables_of_interest:
        raise ValueError("variables_of_interest cannot be empty")

    # Convert variables_of_interest to list if it's a string
    if isinstance(variables_of_interest, str):
        variables_of_interest = [variables_of_interest]
    elif not isinstance(variables_of_interest, list):
        variables_of_interest = list(variables_of_interest)
    
    # Handle selection_pool parameter
    if selection_pool:
        if isinstance(selection_pool, str):
            selection_pool_str = selection_pool
        else:
            selection_pool_list = list(selection_pool)
            selection_pool_str = ','.join(selection_pool_list)
    else:
        selection_pool_str = ''

    # Default R script path
    if r_script_path is None:
        r_script_path = '/Users/alichemkhi/Desktop/myProjects/miic_helper/src/run_miic_select.R'
    
    # Convert variables_of_interest to comma-separated string
    variables_str = ','.join(variables_of_interest)
    
    
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
    
    # Build the R script command (multi-line with \)

    
    cmd = [
        'Rscript', r_script_path, \
        '--matrix', count_matrix_path, \
        '--metadata', metadata_path, \
        '--variables_of_interest', variables_str, \
        '--output', output_path
    ]

    if selection_pool_str:  # Only add selection_pool if it's not empty
        cmd.extend(['--selection_pool', selection_pool_str])

    if compute_meta_mi_scores:
        cmd.extend(['--compute_meta_mi_scores', str(compute_meta_mi_scores).upper()])

    print(f"Running MIIC selection with command:")
    print(' '.join(cmd))
    print(f"Variables of interest: {variables_str}")
    print(f"Selection pool: {selection_pool_str}")
    print(f"Output will be saved to: {output_path}")
    
    if  run_script:
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
    else:
        # print the command without executing in a way that it can be copied and run later
        print("\nCommand to run MIIC selection (not executed):")
        cmd_str = " \\\n    ".join(cmd)
        print(cmd_str)
        return output_path
    