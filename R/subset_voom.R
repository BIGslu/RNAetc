data <- kimma::example.voom

data$targets
data$genes
data$E
data$weights
data$design


# name weights columns
# 1. identify libraries
# provide name of column and list of levels
# alternative to subset by genes

# example filter: "virus == 'none' & asthma == 'healthy'"

subset_voom <- function(dat_voom, lib_filter = NULL, gene_filter=NULL){ # add library ID column in targets?

  dat_voom_sub<-dat_voom
  # Add gene & sample ID's to weights matrix
  # NOTE: may want to add checks for dims
  rownames(dat_voom_sub$weights)<-rownames(dat_voom$E)
  colnames(dat_voom_sub$weights)<-colnames(dat_voom$E)

  #identify libraries to subset to

  if(is.null(lib_filter)){
    libs_sub <- colnames(dat_voom_sub$E) # all libraries if no sublist provided
  } else if(any(lib_filter %in% colnames(dat_voom_sub$E))){
    libs_sub <- intersect(lib_filter, colnames(dat_voom_sub$E)) # add warning for if not ALL libIds from input in voom
  } else{
    # warning that your list is nonsense
  }


  if(is.null(gene_filter)){
    genes_sub <- rownames(dat_voom_sub$E) # all genes if no sublist provided
  } else if(any(gene_filter %in% rownames(dat_voom_sub$E))){
    genes_sub <- intersect(genes_filter, rownames(dat_voom_sub$E)) # add warning for if not ALL gene names from input in voom

    # if we want to include HGNC symbols, think about how to do that
  } else{
    # warning that your list is nonsense
  }

}
