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

subset_voom <- function(dat_voom, lib_filter =NULL, gene_filter=NULL){

  dat_voom_sub<-dat_voom
  # Add gene & sample ID's to weights matrix
  # NOTE: may want to add checks for dims
  rownames(dat_voom_sub$weights)<-rownames(dat_voom$E)
  colnames(dat_voom_sub$weights)<-colnames(dat_voom$E)

  #identify libraries to subset to

  if(is.null(gene_filter))

  if(is.null(gene_filter))

}
