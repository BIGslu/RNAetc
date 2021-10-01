#' Read in counts table from Subread featureCounts
#'
#' @param file Character string for file path of counts table
#'
#' @return data frame
#' @export

read_fcounts <- function(file){
  Geneid <- Chr <- Length <- NULL

  count <- readr::read_tsv(file, skip = 1, show_col_types = FALSE) %>%
    dplyr::select(-c(Chr:Length)) %>%
    dplyr::rename_all(~gsub("_filter.bam","", basename(.))) %>%
    dplyr::rename(ensemblID = Geneid)

  return(count)
}
