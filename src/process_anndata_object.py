"""
process_anndata_object.py

Author: Ali Chemkhi

This module provides utility functions for working with AnnData objects, specifically:
- Creating AnnData objects from pandas DataFrames (count matrix + metadata)
- Saving AnnData objects back to CSV files for external analysis
- Handling format conversions between different data representations

Designed for single-cell RNA-seq data processing workflows.
"""

import pandas as pd
import anndata as ad
import scipy.sparse as sp
import os
import numpy as np
import scanpy as sc

from scipy import io
from scipy.sparse import csr_matrix




def generate_miic_files(
    adata,
    selection,
    metadata,
    foi,
    outdir,
    prefix_name,
    groups={},
    assert_consequence=[],
    contextual=[]
):
    """
    Generate MIIC input files (matrix and metadata) from AnnData object and selection.

    Parameters:
        adata: AnnData object
        selection: List of selected gene/feature names (strings)
        metadata: Dict[str, Optional[List[str]]], where value is None (continuous) or a list of levels (categorical)
        foi: List of features of interest (reference variables)
        outdir: Output directory for files
        prefix_name: Prefix for output files
        groups: Optional dict mapping group names to lists of variable names
        assert_consequence: Optional list of variable names to mark as consequences
        contextual: Optional list of variable names to mark as contextual
    """

    raw_count_matrix = adata.raw.X if adata.raw is not None else adata.X
    gene_names = adata.raw.var_names if adata.raw is not None else adata.var_names
    selected_genes = list(set(selection))
    print("Unique selected genes size (removed metadata):", len(selected_genes))
 
    # Remove metadata from selected genes
    foi_genes = [f for f in foi if f not in metadata] 
    selected_genes = list(set(selected_genes + foi_genes))
    nb_genes = len(selected_genes) # final

    for f in foi:
        assert f in selected_genes or f in metadata, f"Feature of interest '{f}' not in selected genes or metadata."

    # Build matrix
    matrix_use = pd.DataFrame(
        raw_count_matrix[:, [gene_names.get_loc(g) for g in selected_genes]].todense()
        if hasattr(raw_count_matrix, 'todense')
        else raw_count_matrix[:, [gene_names.get_loc(g) for g in selected_genes]],
        columns=selected_genes,
        index=adata.obs_names
    )

    for meta_key, levels in metadata.items():
        col_data = adata.obs[meta_key]
        assert((matrix_use.index == col_data.index).all()), "Index mismatch between matrix_use and adata.obs_names"
        matrix_use[meta_key] = pd.Categorical(col_data, categories=levels) if levels is not None else col_data

    # Save CSV matrix
    outfile_miic = os.path.join(outdir, f"{prefix_name}.{nb_genes}.csv")
    matrix_use.to_csv(outfile_miic, index=False)

    # ----------------- State order of file -----------------

    outfile_contextual = os.path.join(outdir, f"{prefix_name}.{nb_genes}.st.txt")
    category_order_use = pd.DataFrame({
        "var_names": matrix_use.columns,
        "var_type": 1,
        "levels_increasing_order": [np.nan] * matrix_use.shape[1],
        "group": "gene",
        "is_contextual": 0,
        "is_consequence": 0
    })

    for meta_key, levels in metadata.items():
        if levels is not None and meta_key in matrix_use.columns:
            category_order_use.loc[category_order_use.var_names == meta_key, "var_type"] = 0
            category_order_use.loc[category_order_use.var_names == meta_key, "levels_increasing_order"] = ','.join(levels)

    for var in contextual:
        category_order_use.loc[category_order_use.var_names == var, "is_contextual"] = 1

    for var in assert_consequence:
        category_order_use.loc[category_order_use.var_names == var, "is_consequence"] = 1

    # default metadata group assignment
    category_order_use.loc[category_order_use.var_names.isin(metadata.keys()), "group"] = "metadata"

    for group_name, group_vars in groups.items():
        category_order_use.loc[category_order_use.var_names.isin(group_vars), "group"] = group_name

    category_order_use.to_csv(outfile_contextual, index=False, sep='\t')

    print("MIIC files generated:")
    print("  - Matrix:", outfile_miic)
    print("  - Structure:", outfile_contextual)



def load_anndata_from_folder(folder_path):
    """
    Load a Scanpy-compatible dataset saved by `save_seurat_files_mtx()` in R.

    Parameters
    ----------
    folder_path : str
        Path to the folder containing the files:
        - raw_counts.mtx
        - genes.txt
        - cells.txt (optional, can be inferred from metadata)
        - metadata.csv

    Returns
    -------
    adata : scanpy.AnnData
        AnnData object with raw counts, metadata, and gene names.
    """
    import os

    mtx_path = os.path.join(folder_path, "raw_counts.mtx")
    genes_path = os.path.join(folder_path, "genes.txt")
    metadata_path = os.path.join(folder_path, "metadata.csv")

    # Load sparse matrix
    X = io.mmread(mtx_path).tocsr()

    # Load gene names
    var = pd.read_csv(genes_path, header=None, names=["gene_symbols"])
    var.index = var["gene_symbols"]
 

    # Load cell metadata
    obs = pd.read_csv(metadata_path, index_col=0)
    # Construct AnnData object
    adata = sc.AnnData(X=X, obs=obs, var=var)

    return adata





def create_anndata_from_dataframes(count_matrix_df, metadata_df):
    """
    Create an AnnData object from count matrix and metadata dataframes.
    
    Parameters:
    -----------
    count_matrix_df : pd.DataFrame
        Count matrix dataframe (genes x cells or cells x genes)
    metadata_df : pd.DataFrame
        Metadata dataframe (cells x metadata fields)
    
    Returns:
    --------
    adata : AnnData
        AnnData object with count matrix as X and metadata as obs
    """
    try:
        # Make copies to avoid modifying original dataframes
        count_matrix = count_matrix_df.copy()
        metadata = metadata_df.copy()
        
        print(f"Count matrix shape: {count_matrix.shape}")
        print(f"Metadata shape: {metadata.shape}")
        
        # Check if we need to transpose the count matrix
        # AnnData expects cells x genes format
        if count_matrix.shape[1] == metadata.shape[0]:
            print("Transposing count matrix to cells x genes format")
            count_matrix = count_matrix.T
        elif count_matrix.shape[0] != metadata.shape[0]:
            raise ValueError(f"Dimension mismatch: count matrix has {count_matrix.shape[0]} rows, metadata has {metadata.shape[0]} rows")
        
        
        # Ensure indices match
        common_indices = count_matrix.index.intersection(metadata.index)
        if len(common_indices) == 0:
            raise ValueError("No common indices found between count matrix and metadata")
        
        print(f"\nFound {len(common_indices)} common indices")
        
        # Filter to common indices
        count_matrix = count_matrix.loc[common_indices]
        metadata = metadata.loc[common_indices]
        
        # Create AnnData object
        adata = ad.AnnData(X=count_matrix)
        adata.obs = metadata
        
        # Add basic information
        print(f"\nAnnData object created successfully!")
        print(f"Shape: {adata.shape} (cells x genes)")
        print(f"Available metadata columns: {list(adata.obs.columns)}")
        
        return adata
        
    except Exception as e:
        print(f"Error creating AnnData: {str(e)}")
        raise


def save_anndata_to_csv(adata, output_prefix, save_counts=True, save_metadata=True, transpose_counts=False):
    """
    Save both count matrix and metadata from AnnData object to CSV files.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object containing the count matrix and metadata
    output_prefix : str
        Prefix for output files (will add '_counts.csv' and '_metadata.csv')
    save_counts : bool
        Whether to save the count matrix (default: True)
    save_metadata : bool
        Whether to save the metadata (default: True)
    transpose_counts : bool
        If True, transpose the count matrix (genes as rows, cells as columns)
        If False (default), keep AnnData format (cells as rows, genes as columns)
    
    Returns:
    --------
    tuple
        (counts_df, metadata_df) or subset based on what was saved
    """
    
    # Save count matrix if requested
    if save_counts:
        # Check if X is sparse matrix and convert to dense
        if sp.issparse(adata.X):
            counts_matrix = adata.X.toarray()
        else:
            counts_matrix = adata.X
        
        # Create DataFrame with proper indices and column names
        counts_df = pd.DataFrame(
            counts_matrix,
            index=adata.obs.index,  # Cell names as row indices
            columns=adata.var.index  # Gene names as column names
        )
        
        # Transpose if requested
        if transpose_counts:
            counts_df = counts_df.T
        
        counts_path = f"{output_prefix}_counts.csv"
        counts_df.to_csv(counts_path)
        
        print(f"Count matrix saved to: {counts_path}")
        print(f"Shape: {counts_df.shape} ({'genes x cells' if transpose_counts else 'cells x genes'})")
    
    # Save metadata if requested
    if save_metadata:
        # Get metadata from obs
        metadata_df = adata.obs.copy()
        
        metadata_path = f"{output_prefix}_metadata.csv"
        metadata_df.to_csv(metadata_path)
        
        print(f"Metadata saved to: {metadata_path}")
        print(f"Shape: {metadata_df.shape} (cells x metadata_columns)")
        print(f"Metadata columns: {list(metadata_df.columns)}")
    
    return counts_path,metadata_path