############################################################################################################################################################
# Author : achemkhi ali.chemkhi@curie.fr

# Helper functions to preprocess seurat object for miic analysis
############################################################################################################################################################


library(argparse)
library(Seurat)
library(dplyr)
library(grid)
library(gridExtra)
library(fs)
library(millefy)
library(ggplot2)
library(miic) #devtools::install_github("interactMIIC/causalCCC/miic_R_package", force = T, ref = "achemkhi/no-goi-limit-branch") #latest version of MIIC
source("/Users/alichemkhi/Desktop/myProjects/miic_helper/src/aux.r") # TODO add the function in the same script 


# Filter_seurat_object: Filter a Seurat object based on a list of filters (cell type, etc.)
# Args:
#   seurat_object: A Seurat object to filter.
#   filter_list: Named list of filters (e.g., cell types or clusters to keep).
Filter_seurat_object <- function(seurat_object, filter_list) {

  if (length(filter_list) > 0) {  # Apply filters only if filter_list is not empty
    for (name in names(filter_list)) {
      print(paste0("filtered ",name," colomn. keeping only ",paste(filter_list[[name]],collapse = " and ") ))
      Idents(seurat_object)<-name
      print(name)
      seurat_object <- subset(seurat_object, idents = filter_list[[name]])
    }
  }
  
  # A workaround to make MIselection function work. Because it needs a cell type filter 
  seurat_object$myIdentFilter<-"oneinteract"
  
  return(seurat_object)
}


# MI_heatmap: Plot and save a heatmap of mutual information values, highlighting selected genes
MI_heatmap <- function(heatmap_selection,highlight_genes,plot_name,output_path) {
  
  # Add a column to determine label colors based on highlight_genes
  heatmap_selection$label_color <- ifelse(heatmap_selection$variables %in% highlight_genes, "red", "black")
  
  breaks <- c(0, 1000, 2000, Inf)
  
  # Plot heatmap
  plot<-ggplot(heatmap_selection, aes(x = foi, y = variables, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("white","pink", "red", "darkred"), 
      name = "MI"
    )+
    theme_minimal() +
    labs(title = plot_name, x = "reference variables", y = "selected genes") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size=6), # Rotate x-axis labels
      axis.text.y = element_blank(), # Increase y-axis label size for readability
    )+
    geom_text(aes(label = variables, y = variables, x = 0, color = label_color), 
              inherit.aes = FALSE, hjust = 0,size=2) +
    scale_color_identity() # Use label_color as is without legend
  
  ggsave(output_path, plot = plot, width = 6, height = 10, bg = "white") 
  
  return(plot)
}



# handle_MI_Selection: Select top variables by MI, save selection, and plot heatmap
handle_MI_Selection <- function(mi_df,
                                foi,
                                topN,
                                outdir,
                                groupby=FALSE,
                                topn=10,
                                colorFeatures=c(""))
{
  # init to empty
  selection_df <- data.frame()
  if (groupby)
  {
     selection_df<-mi_df %>%
       ungroup() %>%
       filter(foi %in% myfoi) %>% 
       filter (!foi %in% c("RiboSMean","RiboLMean","mtMean")) %>%
       group_by(foi) %>% 
       arrange(desc(value)) %>%
       slice_head(n = topn) %>%
       ungroup() %>%
       arrange(desc(value)) %>%
       distinct(variables,.keep_all = TRUE) %>%
       slice_head(n = topN)     
  }
  else
  {
    selection_df<-mi_df %>%
      ungroup() %>%
      filter(foi %in% myfoi) %>% 
      filter (!foi %in% c("RiboSMean","RiboLMean","mtMean")) %>%
      arrange(desc(value)) %>%
      slice_head(n = topN) %>%
      distinct(variables,.keep_all = TRUE)
  }
  
  selection <- selection_df %>% pull(variables)
  write.csv(selection_df,file = file.path(outdir,paste0("Selection.",length(selection),factor(groupby),".csv")))
  # heatmap 
  file_name<-paste("heatmap",length(selection),factor(groupby),"highlighted.png",sep = ".")
  heatmap_selection<-mi_df[(mi_df$foi %in% myfoi & mi_df$variables %in% selection),] 
  outpath=paste0(outdir,"/",file_name)
  MI_heatmap(heatmap_selection,colorFeatures,"",outpath)
  
  return(selection_df)
}
    


# Generate_miic_files: Generate MIIC input files (matrix and metadata) from Seurat object and selection
# Args:
#   seurat_object: Seurat object containing data.
#   selection: Vector of selected gene/feature names.
#   metadata: Named list of metadata variables and their levels (or NULL for continuous).
#   foi: Vector of features of interest (reference variables).
#   OUTDIR: Output directory for files.
#   prefix_name: Prefix for output files.
#   groups: Named list of groupings for variables (optional).
#   assert_consequence: Vector of variables to mark as consequences (optional).
#   contextual: Vector of variables to mark as contextual (optional).
Generate_miic_files <- function(
    seurat_object,
    selection,
    metadata,
    foi,
    OUTDIR,
    prefix_name,
    groups=list(),
    assert_consequence=c(),
    contextual=c()
)
{
  raw.count.matrix=seurat_object@assays$RNA$counts
  selected_genes<-unique(selection) 
  metadata_selected<-selected_genes[selected_genes %in% c(names(metadata))]
  selected_genes<-selected_genes[!selected_genes %in% c(names(metadata))]
  foi_genes = foi[!foi%in% names(metadata)]
  # Add foi (genes) - make sure they are selected 
  print("Unique selected genes size (removed metadata):")
  print(length(selected_genes))
  print("Selected metadata :")
  print(metadata_selected)
  print("Selected genes INTER foi :")
  print(intersect(foi_genes,selected_genes))
  
  selected_genes<-unique(c(foi_genes,selected_genes))
  nb_genes=length(selected_genes)

  
  
  
  lapply(foi,function(foi_item){
    # foi variable must be either a metadata or a gene 
    stopifnot(foi_item %in% selected_genes | foi_item %in% names(metadata))
  })
  
  # -----------dump csv matrix-------------------------------------------------
  outfile_miic <- paste0(OUTDIR,"/",prefix_name,".",nb_genes,".csv")
  
  matrix.use <- as.data.frame(t(raw.count.matrix[selected_genes,]))
  
  # Add foi(metadata) and selected metadata
  for (metadata_item in names(metadata)) {
    if (metadata_item %in% foi | metadata_item %in% metadata_selected) {
      
      if (is.null(metadata[[metadata_item]])) { # discrete
        matrix.use[[metadata_item]] <- seurat_object@meta.data[[metadata_item]]
      } else { # categorical
        matrix.use[[metadata_item]] <- as.factor(seurat_object@meta.data[[metadata_item]])
      }
      
    }
  }
  
  write.csv(matrix.use,outfile_miic,row.names = FALSE)
  
  
  # --------- dump categorical order
  outfile_contextual <- paste0(OUTDIR,"/",prefix_name,".",nb_genes,".st.txt")
  category.order.use <- data.frame(var_names=colnames(matrix.use),
                                   var_type = rep(1),
                                   levels_increasing_order = rep(NA),
                                   group = "gene",
                                   is_contextual = rep(0),
                                   is_consequence = rep(0))
  

  #define categorical vars and orders
  for (metadata_item in names(metadata)) {
    if (!is.null(metadata[[metadata_item]]) & metadata_item %in% colnames(matrix.use)) {
      category.order.use[category.order.use$var_names == metadata_item, "var_type"] <- 0
      category.order.use[category.order.use$var_names == metadata_item, "levels_increasing_order"] <- metadata[[metadata_item]]  # c("CTL,RES")
    }
  }
  print("category order use")
  #define contextual
  if (length(contextual)>0)
  {
    for (item in contextual) {
      category.order.use[category.order.use$var_names == item, "is_contextual"] <- 1
    }
  }
  
  #define consequences
  if (length(assert_consequence)>0)
  {
    category.order.use<-category.order.use %>%
      mutate(is_consequence = if_else(var_names %in% assert_consequence,1, is_consequence))
    
  }
  
  #define groups
  category.order.use<-category.order.use %>%
    mutate(group = if_else(var_names %in% names(metadata), "metadata", group))

  if (length(groups)>0)
  {
    for (item in names(groups)) {  
      category.order.use <- category.order.use %>%
        mutate(group = if_else(var_names %in% groups[[item]], item, group))
    }
  } 
  
  write.table(category.order.use,outfile_contextual,row.names = FALSE)
  
}