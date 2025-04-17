

# convert from Seurat to anndata format https://smorabit.github.io/blog/2021/velocyto/

library(Seurat)
# CD to the miic folder
#R CMD INSTALL --library=/Users/alichemkhi/Desktop/code/miic_lib_210 ./
library(miic ,lib.loc="/Users/alichemkhi/Desktop/code/miic_lib_210") 
# check miic version
packageVersion("miic")
library(reshape2)

# doc
#' @title wrap_selection
#' @description This function selects features based on mutual information scores.
#' @param dataset Seurat object containing the dataset
#' @param subset_vars Vector of variable names to be used for feature selection (selection pool)
#' @param var_of_interest Variable of interest for feature selection
#' @param threads Number of threads to use for parallel processing
#' @param run_outdir Output directory for saving results
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



# Parse spliced RNA counts 
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


# Parse spliced RNA counts 
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


#Function straining
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


# Function Doublet removal
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

myQCAnnotation <- function(Object){
  message("Adding QC metadata to the seurat object")
  Object[["percent.mt"]] <- PercentageFeatureSet(Object, pattern = "^mt-")
  Object[["percent.Rps"]] <- PercentageFeatureSet(Object, pattern = "^Rps")
  Object[["percent.Rpl"]] <- PercentageFeatureSet(Object, pattern = "^Rpl")
  
  return(Object)
}



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