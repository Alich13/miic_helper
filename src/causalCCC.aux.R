# generate interact_edges 
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


"https://gist.github.com/crazyhottommy/4e46298045a329b47669"
human_to_mouse=read.csv(file="/Users/alichemkhi/Desktop/data/human_mouse_1to1_orthologs.csv")


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