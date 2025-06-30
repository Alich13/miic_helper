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
# --selection_pool: Comma-separated list of genes to consider for selection (subset of matrix columns)
# --output: Path for output CSV file containing mutual information scores
#
# Example usage:
# Rscript /Users/alichemkhi/Desktop/myProjects/miic_helper/src/run_miic_select.R \
#   --matrix  /Users/alichemkhi/Desktop/myProjects/miic_helper/data/wnt_data_counts.csv \
#   --metadata /Users/alichemkhi/Desktop/myProjects/miic_helper/data/wnt_data_metadata.csv \
#   --variables_of_interest T,condition,time \
#   --selection_pool T,WNT3,BMP4 \
#   --output /Users/alichemkhi/Desktop/myProjects/output/v0.2/miic_results.csv


# Load libraries
#library(miic)
library(miic,lib.loc = "/Users/alichemkhi/Desktop/code/miic_v210_branck_lib") 
library(Matrix)
print("----> miic version:")
print(packageVersion("miic"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop("Usage: Rscript run_miic_select.R --matrix matrix.csv --metadata metadata.csv --variables_of_interest variables_csv_string --selection_pool selection_pool_csv_string --output output.csv")
}

# Parse named arguments
matrix_file <- args[which(args == "--matrix") + 1] # Path to the gene expression matrix CSV file
metadata_file <- args[which(args == "--metadata") + 1]  # Path to the metadata CSV file
variables_str <- args[which(args == "--variables_of_interest") + 1] # Comma-separated list of variables of interest
selection_pool_str <- args[which(args == "--selection_pool") + 1] # Comma-separated list of genes to consider for selection
output_file <- args[which(args == "--output") + 1]  # Path for output CSV file

# Convert comma-separated strings to vectors
all_variables <- unlist(strsplit(variables_str, ","))
selection_pool <- unlist(strsplit(selection_pool_str, ","))


# Read matrix and metadata
matrix_df <- read.csv(matrix_file,check.names = FALSE, row.names = 1)
meta_df   <- read.csv(metadata_file,check.names = FALSE, row.names = 1)

# Classify variables into genes (matrix columns) and metadata columns
var_of_interest_names <- c()  # genes in matrix
var_of_interest_value_names <- c()  # metadata columns
var_of_interest_values <- data.frame()  # metadata values



for (var in all_variables) {
  if (var %in% colnames(matrix_df)) {
    # Variable is a gene in the matrix
    var_of_interest_names <- c(var_of_interest_names, var)
    #cat("Found variable '", var, "' in matrix columns\n")
  } else if (var %in% colnames(meta_df)) {
    # Variable is a metadata column
    var_of_interest_value_names <- c(var_of_interest_value_names, var)
    #cat("Found variable '", var, "' in metadata columns\n")
  } else {
    # Variable not found in either
    stop("Warning: Variable '", var, "' not found in matrix columns or metadata columns\n")
  }
}

# Extract metadata values for the found metadata columns
if (length(var_of_interest_value_names) > 0) {
  var_of_interest_values <- meta_df[, var_of_interest_value_names, drop = FALSE]
}

cat("Gene variables:", paste(var_of_interest_names, collapse = ", "), "\n")
cat("Metadata variables:", paste(var_of_interest_value_names, collapse = ", "), "\n")
cat("Selection pool:", paste(selection_pool, collapse = ", "), "\n")

# Filter matrix to keep only genes of interest and selection pool genes
if (length(selection_pool) > 0 && !is.na(selection_pool[1]) && selection_pool[1] != "") {
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


# Run MIIC
mi_scores <- miic:::compute_mi_batch(
  matrix_df,
  var_of_interest_values = var_of_interest_values,
  var_of_interest_names = var_of_interest_names,
  n_threads = 6,
  verbose = 3,
  corrected = FALSE,
  unit="bits"
)

# Save result
cat("output path : ", output_file)
write.csv(mi_scores, output_file, row.names = TRUE)




