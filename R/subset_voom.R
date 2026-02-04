#' Subset voom object
#'
#' @param dat limma EList object to subset
#' @param lib_keep Character vector of library IDs to keep (Default NULL)
#' @param lib_remove Character vector of library IDs to remove (Default NULL)
#' @param lib_filter Character string to use for filtering libraries to keep (Default NULL)
#' @param gene_keep Character vector of genes to keep (Default NULL)
#' @param libraryID Character string specifying the name of the column with library IDs (Default "libID")
#'
#' @return limma EList object
#' @export
#'
#' @examples
#' dat.voom <- RNAetc::example.voom
#'
#' subset_voom(dat.voom, lib_keep = c("lib1","lib2"))
#' subset_voom(dat.voom, lib_remove = c("lib1","lib2"))
#' subset_voom(dat.voom, lib_filter = "asthma == 'healthy' & virus == 'none'")
#' subset_voom(dat.voom, gene_keep = c("ENSG00000000460", "ENSG00000001460"))
#' subset_voom(dat.voom, lib_keep = c("lib1","lib2"),
#'            gene_keep = c("ENSG00000000460", "ENSG00000001460"))

subset_voom <- function(dat,
                        lib_keep = NULL, lib_remove = NULL, lib_filter = NULL,
                        gene_keep=NULL,
                        libraryID = "libID"){ # add library ID column in targets?
  ### Checks
  # check which library filter parameters are NULL, only want either one or none to be specified
  if(sum(is.null(lib_keep),is.null(lib_remove),is.null(lib_filter)) <= 1){
    stop("You specified more than one way to filter libraries, please only specify either lib_keep, lib_remove, or lib_filter")
  }

  #create subset object
  dat_voom_sub<-dat

  # Add gene & sample ID's to weights matrix
  # first check dims are same
  if(!is.null(dat$weights)){
    if(!identical(dim(dat$E),dim(dat$weights))){
      stop("The expression matrix and weights matrix in your voom object have different dimensions...")
    }
    rownames(dat_voom_sub$weights)<-rownames(dat$E)
    colnames(dat_voom_sub$weights)<-colnames(dat$E)
  }

  #identify libraries to subset to:
  if(is.null(lib_keep) & is.null(lib_remove) & is.null(lib_filter)){
    libs_sub <- colnames(dat_voom_sub$E) # all libraries if no library filtering parameter provided
  } else if(!is.null(lib_keep)){
    if(all(lib_keep %in% colnames(dat_voom_sub$E))){
      libs_sub <- lib_keep # just the listed libraries (if they exist in the voom object)
    } else{
      missing <- lib_keep[!(lib_keep %in% colnames(dat_voom_sub$E))]
      stop("I didn't find ",length(missing)," of the libraries you want to keep in the voom object: ",paste0(missing,collapse = ", "))
    }
  } else if(!is.null(lib_remove)){
    if(all(lib_remove %in% colnames(dat_voom_sub$E))){
      libs_sub <- colnames(dat_voom_sub$E)[!(colnames(dat_voom_sub$E) %in% lib_remove)] # all libraries minus the ones to remove (if they exist in the voom object)
    } else{
      missing <- lib_remove[!(lib_remove %in% colnames(dat_voom_sub$E))]
      stop("I didn't find ",length(missing)," of the libraries you want to remove in the voom object: ",paste0(missing,collapse = ", "))
    }
  } else if(!is.null(lib_filter)){
    libs_sub <- dat_voom_sub$targets %>%
      dplyr::filter(!!rlang::parse_expr(lib_filter)) %>%
      dplyr::pull(!!libraryID)
    if(length(libs_sub) == 0){
      stop("Your filter statement didn't find any libraries matching the criteria, double check your statement")
    }
  }

  # identify genes to keep
  if(is.null(gene_keep)){
    genes_sub <- rownames(dat_voom_sub$E) # all genes if no sublist provided
  } else if(all(gene_keep %in% rownames(dat_voom_sub$E))){
    genes_sub <- gene_keep # just the listed libraries (if they exist in the voom object)
  } else{
    missing <- gene_keep[!(gene_keep %in% rownames(dat_voom_sub$E))]
    stop("I didn't find ",length(missing)," of the genes you want to keep in the voom object: ",paste0(missing,collapse = ", "))
  }

  #print message for user to see how many libraries and genes are being included in subset
  message(paste0("Subsetting to ",length(libs_sub)," of ",ncol(dat$E)," libraries"))
  message(paste0("Subsetting to ",length(genes_sub)," of ",nrow(dat$E)," genes"))

  dat_voom_sub$targets <- dat_voom_sub$targets %>% dplyr::filter({{libraryID}} %in% libs_sub)
  dat_voom_sub$genes <- dat_voom_sub$genes[genes_sub,]
  dat_voom_sub$E <- dat_voom_sub$E[genes_sub,libs_sub]
  if(!is.null(dat$weights)){ dat_voom_sub$weights <- dat_voom_sub$weights[genes_sub,libs_sub] }
  if(!is.null(dat$design)){ dat_voom_sub$design <- dat_voom_sub$design[libs_sub,] }

  message("Done!")

  return(dat_voom_sub)
}
