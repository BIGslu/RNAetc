#' Combine 2 voom objects into 1
#'
#' @param dat1 EList voom object resulting from voom() or voomWithQualityWeights()
#' @param dat2 EList voom object resulting from voom() or voomWithQualityWeights()
#' @param libraryID Character. Column name in dat$targets that matches column names in dat$E
#'
#' @return data frame
#' @export
#'
#' @examples
#' example.voom <- example.voom
#' example.voom$design <- NULL
#' dat1 <- subset_voom(example.voom, lib_filter="virus=='none'")
#' dat2 <- subset_voom(example.voom,
#'                lib_keep=c("lib1","lib2","lib4","lib6","lib8","lib10","lib12"))
#' combined <- combine_voom(dat1, dat2)

combine_voom <- function(dat1, dat2, libraryID = "libID"){
  message(paste("Combining genes:",nrow(dat1$E),"from dat1 and", nrow(dat2$E),"from dat2."))
  message(paste("Combining libraries:",ncol(dat1$E),"from dat1 and", ncol(dat2$E),"from dat2."))

  dat.combined <- dat1

  # check duplicate library ID
  names1 <- colnames(dat1)
  names2 <- colnames(dat2)
  if(any(names1 %in% names2)){
    dups <- names1[names1 %in% names2]
    suppressMessages(dat1_dups <- RNAetc::subset_voom(dat1, lib_keep=dups))
    suppressMessages(dat2_dups <- RNAetc::subset_voom(dat2, lib_keep=dups))
    if(!identical(dat1_dups, dat2_dups)){
      stop(paste("    ", paste(dups,collapse=", "), "exist in both voom objects and do not match. Please correct."))
    } else{
      message(paste("    ", paste(dups,collapse=", "), "exist in both voom objects. Only 1 copy will be retained."))
    }
  } else{
    dups <- NULL
  }

  #Counts
  dat.combined$E <- as.data.frame(dat1$E) %>%
    tibble::rownames_to_column() %>%
    dplyr::full_join(as.data.frame(dat2$E) %>%
                       tibble::rownames_to_column() %>%
                       dplyr::select(-dups),
                     by = "rowname") %>%
    tibble::column_to_rownames()

  #Targets
  dat.combined$targets <- dat1$targets %>%
    dplyr::bind_rows(dat2$targets %>% dplyr::filter(!!dplyr::sym(libraryID) %in% dups))

  #Genes
  genes_total <- sum(!is.null(dat1$genes), !is.null(dat2$genes))
  if(genes_total==2){
    suppressMessages(dat.combined$genes <- dplyr::full_join(dat1$genes, dat2$genes))
  } else if(genes_total==1){
    stop("Genes only found in one voom object. Cannot combine.")
  }else {
      dat.combined$genes <- NULL
    }

  #Weights
  weights_total <- sum(!is.null(dat1$weights), !is.null(dat2$weights))
  if(weights_total==2){
    rownames(dat1$weights)<-rownames(dat1$E)
    colnames(dat1$weights)<-colnames(dat1$E)
    rownames(dat2$weights)<-rownames(dat2$E)
    colnames(dat2$weights)<-colnames(dat2$E)

    dat.combined$weights <- as.data.frame(dat1$weights) %>%
      tibble::rownames_to_column() %>%
      dplyr::full_join(as.data.frame(dat2$weights) %>%
                         tibble::rownames_to_column() %>%
                         dplyr::select(-dups),
                       by = "rowname") %>%
      tibble::column_to_rownames()

  } else if(weights_total==1){
    stop("Weights only found in one voom object. Cannot combine.")
  } else{
    dat.combined$weights <- NULL
  }

  #Design - remove as is not valid for combined object
  dat.combined$design <- NULL

  message(paste("\nReturning voom object with", ncol(dat.combined$E),
                "libraries and", nrow(dat.combined$E), "genes."))
  return(dat.combined)
}
