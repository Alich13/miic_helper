# convert from Seurat to anndata format https://smorabit.github.io/blog/2021/velocyto/

library(Seurat)
# CD to the miic folder
#R CMD INSTALL --library=/Users/alichemkhi/Desktop/code/miic_lib_210 ./
library(miic ,lib.loc="/Users/alichemkhi/Desktop/code/miic_lib_210") 
# check miic version
packageVersion("miic")
library(reshape2)

#' @title wrap_selection
#' @description This function selects features based on mutual information scores.
#' @param dataset Seurat object containing the dataset
#' @param subset_vars Vector of variable names to be used for feature selection (selection pool)
#' @param var_of_interest Variable of interest for feature selection
#' @param threads Number of threads to use for parallel processing
#' @param run_outdir Output directory for saving results
#' @param name Prefix for output files
#' @return A data frame with variables, feature of interest, and MI scores
#' @examples
#' # Not run:
#' wrap_MI_compute(dataset, c("gene1", "gene2"), "outcome", 4, "/path/to/dir", "results")
#' # End(Not run)
wrap_MI_compute <- function(dataset,subset_vars,var_of_interest,threads,run_outdir,name) {
  
  matrix_ <-dataset@assays$RNA$counts
  print("Matrix loaded ...")
  # if subset_vars is not NULL, subset the matrix
  if (!is.null(subset_vars)) {
    matrix_ <- matrix_[subset_vars,]
    print("Matrix subsetted ...")
  }

  # add the metadata of interest to the matrix
  metadata_of_interest<- intersect(var_of_interest, colnames(dataset@meta.data))
  matrix_ <-cbind(t(matrix_), dataset@meta.data[metadata_of_interest])

  print("Computing MI scores ...")
  mi_scores<-miic::selectFeatures(
      matrix_,
      n_features=0,
      var_of_interest_names=unique(var_of_interest),
      n_threads=threads, verbose=3, plot=F)
  print("MI selection done ... saving results")

  matrix_melt<-reshape2::melt(mi_scores$mis) 
  names(matrix_melt) <-  c("variables","foi","value")
  # save selection as tsv
  write.table(matrix_melt,
              file=file.path(run_outdir,paste0(name,"_MI_table.csv")),sep="\t",quote = F,row.names = F)

  return (matrix_melt)
}




#' @title average_gene_count
#' @description Computes the average count of genes matching a given pattern across all cells in the Seurat object, adds the information as metadata, and removes the original genes.
#' @param pattern A regex pattern to match gene names.
#' @param Object A Seurat object.
#' @param name The name of the metadata column to be added.
#' @return The Seurat object with the new metadata column and without the original genes.
#' @examples
#' # Not run:
#' seurat_obj <- average_gene_count("MT-", seurat_obj, "mito_avg_count")
#' # End(Not run)
average_gene_count <- function(pattern,Object,name){
  
  pattern_genes <- rownames(Object)[grepl(pattern, rownames(Object), ignore.case = TRUE)]
  pattern_matrix <- GetAssayData(Object, layer = "count")[pattern_genes,]
  mean_expression <- colMeans(pattern_matrix)
  
  Object <- AddMetaData(
    object = Object,
    metadata = mean_expression, # Mean expression 
    col.name = name )
  
  #Remove the original ribosomal genes
  remaining_genes <- setdiff(rownames(Object@assays$RNA$counts), c(pattern_genes))
  Object <- subset(Object, features = remaining_genes)
  
  return(Object)
}


#' @title load_seurat_data
#' @description Loads 10X Genomics formatted data and creates a Seurat object. If metadata is available, it is added to the Seurat object.
#' @param data_dir Directory containing the 10X Genomics data files: matrix.mtx, features.tsv.gz, barcodes.tsv.gz, and optionally metadata.txt.
#' @return A Seurat object containing the loaded data.
#' @examples
#' # Not run:
#' seurat_obj <- load_seurat_data("/path/to/data_dir")
#' # End(Not run)
load_seurat_data <- function(data_dir) {
  
  # Load expression matrix
  counts <- readMM(file.path(data_dir,"matrix.mtx"))
  rownames(counts) <- readLines(file.path(data_dir,"features.tsv.gz"))  # Assign gene names
  colnames(counts) <- readLines(file.path(data_dir,"barcodes.tsv.gz"))  # Assign cell barcodes
  
  # Check if metadata is provided
  metadata_file<-file.path(data_dir,"metadata.txt")
  if (!is.null(metadata_file)) {
    metadata <- read.table(metadata_file, sep="\t", header=TRUE, row.names=1)
    seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
  } else {
    seurat_obj <- CreateSeuratObject(counts = counts)
  }
  
  seurat_obj@meta.data$orig.ident=seurat_obj@meta.data$myIdent
  seurat_obj<-SetIdent(seurat_obj, value = "orig.ident")
  
  return(seurat_obj)
}



#' @title myParseUnsplicedSO
#' @description Parses unspliced RNA count data and creates a Seurat object.
#' @param dir Directory containing the unspliced data files: unspliced.mtx, unspliced.barcodes.txt, unspliced.genes.txt.
#' @param sample_name Name of the sample, used as the project name in Seurat.
#' @return A Seurat object containing the unspliced RNA data.
#' @examples
#' # Not run:
#' seurat_obj <- myParseUnsplicedSO("/path/to/dir", "sample_name")
#' # End(Not run)
myParseUnsplicedSO <- function(dir, sample_name){
  # Load the matrix
  matrix_data <- readMM(file.path(dir,"unspliced.mtx"))
  # Load the barcodes
  barcodes <- readLines(file.path(dir,"unspliced.barcodes.txt"))
  # Load the gene names
  genes <- read.table(file.path(dir,"unspliced.genes.txt"), header = FALSE, stringsAsFactors = FALSE)
  
  # Assign row and column names
  colnames(matrix_data) <- genes$V1
  rownames(matrix_data) <- barcodes
  
  # Create the Seurat object
  seurat_obj <- CreateSeuratObject(counts = t(matrix_data), project = sample_name)
  
  # Print summary
  return(seurat_obj)
}


#' @title myParseSO
#' @description Parses spliced RNA count data and creates a Seurat object.
#' @param dir Directory containing the spliced data files: spliced.mtx, spliced.barcodes.txt, spliced.genes.txt.
#' @param sample_name Name of the sample, used as the project name in Seurat.
#' @return A Seurat object containing the spliced RNA data.
#' @examples
#' # Not run:
#' seurat_obj <- myParseSO("/path/to/dir", "sample_name")
#' # End(Not run)
myParseSO <- function(dir, sample_name){
  # Load the matrix
  matrix_data <- readMM(file.path(dir,"spliced.mtx"))
  # Load the barcodes
  barcodes <- readLines(file.path(dir,"spliced.barcodes.txt"))
  # Load the gene names
  genes <- read.table(file.path(dir,"spliced.genes.txt"), header = FALSE, stringsAsFactors = FALSE)
  
  # Assign row and column names
  colnames(matrix_data) <- genes$V1
  rownames(matrix_data) <- barcodes
  
  # Create the Seurat object
  seurat_obj <- CreateSeuratObject(counts = t(matrix_data), project = sample_name)
  
  # Print summary
  return(seurat_obj)
}


#' @title myStrainingAmbRNA
#' @description Corrects for ambient RNA using the SoupX algorithm and creates a Seurat object.
#' @param exp Path to the experiment directory containing the raw and filtered matrices.
#' @param sample Sample name, used for naming the Seurat object.
#' @return A Seurat object with corrected RNA counts.
#' @examples
#' # Not run:
#' seurat_obj <- myStrainingAmbRNA("/path/to/exp", "sample_name")
#' # End(Not run)
myStrainingAmbRNA <- function(exp, sample){
  
  raw_matrix <- file.path(exp,sample, "outs/raw_feature_bc_matrix")
  filtered_matrix <- file.path(exp,sample,"outs/filtered_feature_bc_matrix")
  
  message("Matrices are being read ...")
  toc <- Seurat::Read10X(filtered_matrix)
  tod <- Seurat::Read10X(raw_matrix)
  
  srat <- CreateSeuratObject(counts = toc) # maybe set assay next time ?
  soup.channel <- SoupChannel(tod, toc)
  
  # Find Clusters for SoupX
  message("Clustering for SoupeX algorithm ...")
  srat    <- SCTransform(srat, verbose = F)
  srat    <- RunPCA(srat, verbose = F)
  srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat    <- FindClusters(srat, verbose = T) # res 0.8 by default
  
  # Saving metadata for sanity check
  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  
  # Calculating ambient RNA profile
  message("Calculating ambiant RNA profile ...")
  soup.channel <- autoEstCont(soup.channel)
  
  # SoupX matrix correction
  message("SoupX matrix correction ...")
  adj.matrix  <- adjustCounts(soup.channel, roundToInt = T) #roundToInt to make sure we output an integer matrix
  
  # Reincorporating Antibody Capture assay to the matrix
  data <- list("Gene Expression" = adj.matrix)
  
  
  message("Creation of Seurat object...")
  Object <- CreateSeuratObject(counts = data$`Gene Expression`, project = sample, assay = "RNA")
  
  return(Object)
}


#' @title myDoubletRemoval
#' @description Removes doublets from the Seurat object using the scDblFinder package.
#' @param Object A Seurat object.
#' @return A Seurat object with doublets removed.
#' @examples
#' # Not run:
#' seurat_obj <- myDoubletRemoval(seurat_obj)
#' # End(Not run)
myDoubletRemoval <- function(Object){
  # Doublet removal
  message("Doublet removal using scDblFinder ...")
  sce <- SingleCellExperiment(assays=list(counts=Object@assays$RNA@counts))
  message("SingleCellExperiment done")
  sce <- scDblFinder(sce)
  
  sce_name <- paste0(project.name,'_',day,'_sce')
  sce <- get(sce_name)
  
  Object<-Object[,sce@colData$scDblFinder.class == "singlet"]
  
  return(Object)
  
}

#' @title myQCAnnotation
#' @description Adds quality control metrics to the Seurat object, including mitochondrial, ribosomal protein S, and ribosomal protein L gene expression percentages.
#' @param Object A Seurat object.
#' @return A Seurat object with additional metadata columns for QC metrics.
#' @examples
#' # Not run:
#' seurat_obj <- myQCAnnotation(seurat_obj)
#' # End(Not run)
myQCAnnotation <- function(Object){
  message("Adding QC metadata to the seurat object")
  Object[["percent.mt"]] <- PercentageFeatureSet(Object, pattern = "^mt-")
  Object[["percent.Rps"]] <- PercentageFeatureSet(Object, pattern = "^Rps")
  Object[["percent.Rpl"]] <- PercentageFeatureSet(Object, pattern = "^Rpl")
  
  return(Object)
}



#' @title myNormalizeDatabySizeFactor
#' @description Normalizes the Seurat object using size factors estimated from the data.
#' @param SeuratObject A Seurat object.
#' @return A Seurat object with normalized data.
#' @examples
#' # Not run:
#' seurat_obj <- myNormalizeDatabySizeFactor(seurat_obj)
#' # End(Not run)
myNormalizeDatabySizeFactor <- function(SeuratObject){
  ## Normalisation dependant of the theoric size of cells (number of counts and features)
  sce.object <- as.SingleCellExperiment(SeuratObject)
  clusters <- quickCluster(sce.object)
  sce.object <- computeSumFactors(sce.object, clusters = clusters)
  sce.object <- logNormCounts(sce.object)
  
  # Reconverting the sce.object into a seurat object and redo the defaults idents and assays
  SeuratObject <- as.Seurat(sce.object,)
  DefaultAssay(SeuratObject) <- 'RNA'
  
  return(SeuratObject)
}


#' @title myEnhancedVolcano
#' @description Creates a volcano plot for visualizing differentially expressed genes (DEGs) between two groups.
#' @param markers A data frame containing the results of differential expression analysis.
#' @param pvalue_cutoff Adjusted p-value cutoff for significance.
#' @param FC_cutoff Log2 fold change cutoff for significance.
#' @param onepop Optional string to specify a population for labeling in the plot.
#' @param ident.1 Name of the first identity (group) in the Seurat object.
#' @param ident.2 Name of the second identity (group) in the Seurat object.
#' @return A ggplot object representing the volcano plot.
#' @examples
#' # Not run:
#' volcano_plot <- myEnhancedVolcano(markers, 0.05, 0.5, "B cells", "Group1", "Group2")
#' # End(Not run)
myEnhancedVolcano <- function(markers,
                              pvalue_cutoff = 0.05,
                              FC_cutoff = 0.5,
                              onepop = "",
                              ident.1,
                              ident.2
){
  
  
  keyvals.color.p_val <- ifelse(
    markers$avg_log2FC < -FC_cutoff & markers$p_val_adj < pvalue_cutoff, "royalblue" , 
    ifelse(markers$avg_log2FC > FC_cutoff & markers$p_val_adj < pvalue_cutoff, "red2", "grey30"))
  names(keyvals.color.p_val)[keyvals.color.p_val == "red2"] <- paste("Overexpressed genes in", ident.1, onepop)
  names(keyvals.color.p_val)[keyvals.color.p_val == "royalblue"] <- paste("Overexpressed genes in",ident.2, onepop)
  names(keyvals.color.p_val)[keyvals.color.p_val == "grey30"] <- "Not sig."
  
  
  
  plotDEG <- EnhancedVolcano(markers,
                             lab = rownames(markers),
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             legendPosition = 'right',
                             legendLabSize = 15,
                             legendIconSize = 5,
                             legendLabels = c(paste("Overexpressed genes in",ident.2, onepop), paste("Overexpressed genes in",ident.1, onepop), "Not sig."),
                             pointSize = 4.5,
                             colCustom = keyvals.color.p_val,
                             drawConnectors = T,
                             widthConnectors = 0.75,
                             pCutoff = pvalue_cutoff,
                             FCcutoff = FC_cutoff,
                             boxedLabels = T,
                             title = paste('DEG analysis of',ident.1, onepop, 'versus',ident.2, onepop, ' with wilcox test'),
                             subtitle = "Bonferroni adj p-value"
  )  +
    theme(legend.position = "bottom")
  
  return(plotDEG)
}



# generate interact_edges: Generate ligand-receptor interaction edges from a data frame
# Args:
#   L_to_R_df: Data frame with columns 'ligand_complex' and 'receptor_complex'
generate_interact_edges <- function(L_to_R_df) {
  
  interact_edges <- data.frame(ligands = character(), receptors = character())
  
  for (i in 1:nrow(L_to_R_df)) {
    oneligand <- L_to_R_df$ligand_complex[i]
    onereceptor <- L_to_R_df$receptor_complex[i]
    onereceptor <- unique(unlist(str_split(onereceptor, "_")))
    for (onerecp in onereceptor) {
      interact_edges[nrow(interact_edges) +1,] <- c(oneligand, onerecp)
    }
  }
  interact_edges <- interact_edges[!(duplicated(interact_edges)),]
  
  return(interact_edges)
}


#"https://gist.github.com/crazyhottommy/4e46298045a329b47669"
human_to_mouse=read.csv(file="/Users/alichemkhi/Desktop/data/human_mouse_1to1_orthologs.csv")


# convert_human_to_mouse: Convert a list of human gene names to mouse gene name format (capitalize first letter, rest lowercase)
# Args:
#   gene_list: Character vector of human gene names
convert_human_to_mouse <- function(gene_list) {
  # human_to_mouse[human_to_mouse$human %in% ligands,]$mouse
  # missing some genes in the table above 
  # so converting manually and making sure we don't lose genes in the process  
  # mouse - human consistence workaround
  # Convert to lowercase and capitalize the first letter
  # Convert Liana gene names to mouse format (SERPINE1 -> Serpine1)
  
  gene_list_cap <- unname(sapply(gene_list, function(x) {
    paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
  }))
  
  return(gene_list_cap)
}