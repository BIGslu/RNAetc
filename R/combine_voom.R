#' Combine data frame elements of a voom object
#'
#' @param dat EList voom object resulting from voom() or voomWithQualityWeights()
#' @param libraryID Character. Column name in dat$targets that matches column names in dat$E
#' @param geneID Character. Column name in dat$genes that matches row names in dat$E
#'
#' @return data frame
#' @export
#'
#' @examples
#' combined <- combine_voom(example.voom)

combine_voom <- function(dat, libraryID = "libID", geneID = "geneName"){

  combined <- as.data.frame(dat$E) %>%
    tibble::rownames_to_column(geneID) %>%
    tidyr::pivot_longer(-geneID, names_to = libraryID, values_to = "log2cpm") %>%
    dplyr::left_join(dat$targets, by = libraryID)

  if(!is.null(dat$genes)){
    combined <- combined %>%
      dplyr::left_join(dat$genes, by = geneID)
  }

  return(combined)
}
