# miic_helper

This repo contain helper functions to preprocess seurat objects and prepare for miic analysis . 

Planning to wrap up all the functions in a single command line (python) in the future .

- Read a seurat object (need to support adata files also )
- Filter according to the defined filters
- Run MI selection on our features of interest
- Handle the selection (group by foi or not )
- Dump miic preprocessed files  

This repo contains other functions that are useful in miic preprocessing like  `selectMI` from franck and `layout_network.ipynb` . 