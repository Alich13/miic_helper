# MIIC Feature Selection Script
# Author: Ali Chemkhi
# 
# This script runs MIIC feature selection to identify variables that are most informative
# about specified variables of interest, given metadata conditions.
#
# Parameters:
# --matrix: Path to CSV file containing gene expression count matrix (genes as columns)
# --metadata: Path to CSV file containing sample metadata
# --variables_of_interest: Comma-separated list of variables of interest (genes or metadata columns, e.g., "Il10,Phgdh,condition,mtMean")
# --selection_pool: [OPTIONAL] Comma-separated list of genes to consider for selection (subset of matrix columns)
# --output: Path for output CSV file containing mutual information scores
# --compute_meta_mi_scores: Whether to compute MI scores between metadata variables (TRUE/FALSE, default: FALSE)
#
# Example usage:
# Rscript /Users/alichemkhi/Desktop/myProjects/miic_helper/src/run_miic_select.R \
#   --matrix  /Users/alichemkhi/Desktop/myProjects/miic_helper/data/wnt_data_counts.csv \
#   --metadata /Users/alichemkhi/Desktop/myProjects/miic_helper/data/wnt_data_metadata.csv \
#   --variables_of_interest T,condition,time \
#   --output /Users/alichemkhi/Desktop/myProjects/output/v0.2/miic_results.csv \
#   --compute_meta_mi_scores TRUE
#
# Or with selection pool:
# Rscript /Users/alichemkhi/Desktop/myProjects/miic_helper/src/run_miic_select.R \
#   --matrix  /Users/alichemkhi/Desktop/myProjects/miic_helper/data/wnt_data_counts.csv \
#   --metadata /Users/alichemkhi/Desktop/myProjects/miic_helper/data/wnt_data_metadata.csv \
#   --variables_of_interest T,condition,time \
#   --selection_pool T,WNT3,BMP4 \
#   --output /Users/alichemkhi/Desktop/myProjects/output/v0.2/miic_results.csv \
#   --compute_meta_mi_scores TRUE


# Load libraries
#library(miic)
library(miic,lib.loc = "/Users/alichemkhi/Desktop/code/miic_v210_branck_lib") 
library(Matrix)
print("----> miic version:")
print(packageVersion("miic"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: Rscript run_miic_select.R --matrix matrix.csv --metadata metadata.csv --variables_of_interest variables_csv_string --output output.csv [--selection_pool selection_pool_csv_string] [--compute_meta_mi_scores TRUE/FALSE]")
}

# Parse named arguments
matrix_file <- args[which(args == "--matrix") + 1] # Path to the gene expression matrix CSV file
metadata_file <- args[which(args == "--metadata") + 1]  # Path to the metadata CSV file
variables_str <- args[which(args == "--variables_of_interest") + 1] # Comma-separated list of variables of interest
output_file <- args[which(args == "--output") + 1]  # Path for output CSV file

# Parse optional selection_pool parameter
selection_pool_str <- NULL
if ("--selection_pool" %in% args) {
  selection_pool_str <- args[which(args == "--selection_pool") + 1]
}

# Parse optional compute_meta_mi_scores parameter with default value FALSE
compute_meta_mi_scores <- FALSE
if ("--compute_meta_mi_scores" %in% args) {
  compute_meta_mi_scores_str <- args[which(args == "--compute_meta_mi_scores") + 1]
  compute_meta_mi_scores <- as.logical(compute_meta_mi_scores_str)
}

# Convert comma-separated strings to vectors
all_variables <- unlist(strsplit(variables_str, ","))
selection_pool <- if (!is.null(selection_pool_str)) unlist(strsplit(selection_pool_str, ",")) else NULL


# Read matrix and metadata
matrix_df <- read.csv(matrix_file,check.names = FALSE, row.names = 1)
meta_df   <- read.csv(metadata_file,check.names = FALSE, row.names = 1)

# Classify variables into genes (matrix columns) and metadata columns
var_of_interest_names <- c()  # genes in matrix
metadata_names <- c()  # metadata columns
var_of_interest_values <- data.frame()  # metadata values



for (var in all_variables) {
  if (var %in% colnames(matrix_df)) {
    # Variable is a gene in the matrix
    var_of_interest_names <- c(var_of_interest_names, var)
    #cat("Found variable '", var, "' in matrix columns\n")
  } else if (var %in% colnames(meta_df)) {
    # Variable is a metadata column
    metadata_names <- c(metadata_names, var)
    #cat("Found variable '", var, "' in metadata columns\n")
  } else {
    # Variable not found in either
    stop("Warning: Variable '", var, "' not found in matrix columns or metadata columns\n")
  }
}

# Extract metadata values for the found metadata columns
if (length(metadata_names) > 0) {
  var_of_interest_values <- meta_df[, metadata_names, drop = FALSE]
}

cat("Gene variables:", paste(var_of_interest_names, collapse = ", "), "\n")
cat("Metadata variables:", paste(metadata_names, collapse = ", "), "\n")
cat("Selection pool:", if (!is.null(selection_pool)) paste(selection_pool, collapse = ", ") else "None (using all genes)", "\n")
cat("Compute metadata MI scores:", compute_meta_mi_scores, "\n")

# Filter matrix to keep only genes of interest and selection pool genes
if (!is.null(selection_pool) && length(selection_pool) > 0 && !is.na(selection_pool[1]) && selection_pool[1] != "") {
  # Combine genes of interest and selection pool genes
  genes_to_keep <- unique(c(var_of_interest_names, selection_pool))
  
  # Check which genes are actually present in the matrix
  available_genes <- intersect(genes_to_keep, colnames(matrix_df))
  missing_genes <- setdiff(genes_to_keep, colnames(matrix_df))
  
  if (length(missing_genes) > 0) {
    cat("Warning: The following genes are not found in the matrix:", paste(missing_genes, collapse = ", "), "\n")
  }
  
  if (length(available_genes) == 0) {
    stop("Error: None of the specified genes are found in the matrix columns")
  }
  
  # Filter the matrix to keep only available genes
  cat("Filtering matrix to keep", length(available_genes), "genes:", paste(available_genes, collapse = ", "), "\n")
  matrix_df <- matrix_df[, available_genes, drop = FALSE]
  
  cat("Matrix filtered. New dimensions:", nrow(matrix_df), "x", ncol(matrix_df), "\n")
} else {
  cat("No selection pool specified. Using all genes in the matrix.\n")
}


# Run MIIC feature selection
cat("Computing mutual information scores...\n")
# NOTE : MI is not computed between metadata features, only between genes and metadata variables 
mi_scores <- miic:::compute_mi_batch(
  matrix_df,
  var_of_interest_values = var_of_interest_values, # metadata subset - a data frame
  var_of_interest_names = var_of_interest_names,   # genes of interest - should be in matrix_df columns
  n_threads = 6,
  verbose = 3,
  corrected = FALSE,
  unit="bits"
)
# Save result
write.csv(mi_scores, output_file, row.names = TRUE)
cat("output path : ", output_file)

# return mi scores for metadata against each others
if (compute_meta_mi_scores) {
  cat("Computing mutual information scores for metadata variables...\n")
  if (length(metadata_names) > 0) {
    mi_scores_meta <- miic:::compute_mi_batch(
      meta_df,
      var_of_interest_names = metadata_names,
      n_threads = 6,
      verbose = 3,
    corrected = FALSE,
    unit="bits"
  )
  } else {
    mi_scores_meta <- data.frame()  # Empty data frame if no metadata variables
    cat("No metadata variables specified for MI computation.\n")
  }
  # Save metadata MI scores
  meta_output_file <- sub("\\.csv$", "_meta_mi_scores.csv", output_file)
  write.csv(mi_scores_meta, meta_output_file, row.names = TRUE)
  cat("Metadata MI scores saved to:", meta_output_file, "\n")
}






