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

subset_voom <- function(dat_voom,lib_filter,gene_filter){
  #identify libraries to subset to
}
