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

subset_voom <- function(dat_voom, lib_keep = NULL, lib_remove = NULL, lib_filter = NULL, gene_keep=NULL, libraryID = "libID"){ # add library ID column in targets?
  ### Checks
  # check which library filter parameters are NULL, only want either one or none to be specified
  if(sum(is.null(lib_keep),is.null(lib_remove),is.null(lib_filter)) <= 1){
    stop("You specified more than one way to filter libraries, please only specify either lib_keep, lib_remove, or lib_filter")
  }


  #create subset object
  dat_voom_sub<-dat_voom

  # Add gene & sample ID's to weights matrix
  # first check dims are same
  if(!identical(dim(dat_voom$E),dim(dat_voom$weights))){
    stop("The expression matrix and weights matrix in your voom object have different dimensions...")
  }
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
    genes_sub <- intersect(gene_filter, rownames(dat_voom_sub$E)) # add warning for if not ALL gene names from input in voom

    # if we want to include HGNC symbols, think about how to do that
  } else{
    # warning that your list is nonsense
  }

  #print message for user to see how many libraries and genes are being included in subset
  message(paste0("Subsetting to ",length(libs_sub)," of ",ncol(dat_voom$E)," libraries"))
  message(paste0("Subsetting to ",length(genes_sub)," of ",nrow(dat_voom$E)," genes"))

  dat_voom_sub$targets <- dat_voom_sub$targets[libs_sub,]
  dat_voom_sub$genes <- dat_voom_sub$genes[genes_sub,]
  dat_voom_sub$E <- dat_voom_sub$E[genes_sub,libs_sub]
  dat_voom_sub$weights <- dat_voom_sub$weights[genes_sub,libs_sub]
  dat_voom_sub$design <- dat_voom_sub$design[libs_sub,]

  message("Done!")

  return(dat_voom_sub)
}
