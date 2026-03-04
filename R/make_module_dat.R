#' Create dat object with mean module/gene set expression
#'
#' @param module_df Data frame with module names in one column and gene membership in another column
#' @param dat limma EList object to alter
#' @param genecol String showing name of the column in module_df containing the gene IDs. Be sure these IDs match the row names in the input dat object. (Default "ensembl_gene")
#' @param modcol String showing name of the column in module_df containing the module IDs. (Default "gs_name")
#' @param lib String showing the name of the column containing sample IDs in dat$targets. This will become the column names for the new collapsed expression matrix. (Default "libID")
#'
#' @return limma EList object
#' @export
#'
#' @examples
#' dat.voom <- RNAetc::example.voom
#' mod_df_ex <- data.frame("gs_name" = c(rep("M1", 5),
#'                                       rep("M2", 5),
#'                                       rep("M3", 5)),
#'                         "ensembl_gene" = c(rownames(dat.voom$E)[1:15]))
#'
#' make_module_dat(module_df = mod_df_ex, dat = dat.voom)

make_module_dat <- function(module_df = NULL,
                            dat = NULL,
                            genecol = "ensembl_gene",
                            modcol = "gs_name",
                            lib = "libID"){

  . <- module <- gene <- NULL

  # format module df
  conmod_df <- data.frame("gene" = module_df[[genecol]],
                          "module" = module_df[[modcol]])

  # get mod names
  conmods <- unique(conmod_df$module)

  # get mod genes
  conmodgenes <- list()
  for(i in conmods){
    conmodgenes[[i]] <- conmod_df %>%
      dplyr::filter(module == i) %>%
      dplyr::pull(gene)
  }

  # grab libIDs to relate expression to metadata
  base_df <- data.frame("libID" = dat[["targets"]][[lib]])

  for(i in 1:length(conmodgenes)){
    # get gene ids for module
    gl_name <- names(conmodgenes)[i]
    gl_temp <- conmodgenes[[i]]
    gl_temp <- gl_temp[which(gl_temp %in% rownames(dat[["E"]]))]

    # subset to genes that exist in input dat
    etemp <- base::t(dat[["E"]]) %>%
      as.data.frame() %>%
      dplyr::select(dplyr::all_of(gl_temp))

    # sample-wise average expression over module
    jdf <- data.frame("libID" = rownames(etemp))
    jdf[[gl_name]] <- rowSums(etemp)/length(conmodgenes[[i]])

    # join to libID df
    base_df <- base_df %>%
      dplyr::left_join(jdf)
  }

  dat_new <- dat # targets stays the same

  # clean and transpose module-wise expression matrix
  dat_new$E <- base_df %>%
    dplyr::select(dplyr::all_of(c("libID",names(conmodgenes)))) %>%
    tibble::column_to_rownames(var = "libID") %>%
    base::t()

  # make genes (modules) df
  dat_new$genes <- data.frame("gene" = names(conmodgenes))
  dat_new$genes[[modcol]] <- names(conmodgenes)
  rownames(dat_new$genes) <- names(conmodgenes)

  # weights are irrelevant
  dat_new$weights <- NULL

  return(dat_new)
}
