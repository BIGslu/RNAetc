#' Combine all parts of a voom object into a single dataframe
#'
#' @param dat EList voom object resulting from voom() or voomWithQualityWeights()
#' @param libraryID Character. Column name in dat$targets that matches column names in dat$E
#' @param geneID Character. Column name in dat$genes that matches row names in dat$E
#' @param include_weights Logical. If gene-level weights should be included in the result. Default is FALSE
#'
#' @return data frame
#' @export
#'
#' @examples
#' combined <- collapse_voom(example.voom)
#'

collapse_voom <- function(dat,
                          libraryID = "libID", geneID = "geneName",
                          include_weights = FALSE){
  combined <- as.data.frame(dat$E) %>%
    tibble::rownames_to_column(geneID) %>%
    tidyr::pivot_longer(-geneID, names_to = libraryID) %>%
    dplyr::left_join(dat$targets, by=libraryID)

  if(!is.null(dat$genes)){
    combined <- dplyr::left_join(combined, dat$genes, by=geneID)
  }

  if(include_weights){
    colnames(dat$weights) <- colnames(dat$E)
    rownames(dat$weights) <- rownames(dat$E)

    combined <- combined %>%
      dplyr::left_join(as.data.frame(dat$weights) %>%
                         tibble::rownames_to_column(geneID) %>%
                         tidyr::pivot_longer(-geneID, names_to = libraryID,
                                             values_to = "weights"),
                       by=c(libraryID, geneID))
  }
  return(combined)
}
