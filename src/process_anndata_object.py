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
        
        # Save to CSV
        counts_path = f"{output_prefix}_counts.csv"
        counts_df.to_csv(counts_path)
        
        print(f"Count matrix saved to: {counts_path}")
        print(f"Shape: {counts_df.shape} ({'genes x cells' if transpose_counts else 'cells x genes'})")
    
    # Save metadata if requested
    if save_metadata:
        # Get metadata from obs
        metadata_df = adata.obs.copy()
        
        # Save to CSV
        metadata_path = f"{output_prefix}_metadata.csv"
        metadata_df.to_csv(metadata_path)
        
        
        print(f"Metadata saved to: {metadata_path}")
        print(f"Shape: {metadata_df.shape} (cells x metadata_columns)")
        print(f"Metadata columns: {list(metadata_df.columns)}")
    
    return counts_path,metadata_path