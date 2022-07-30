#' Construct WGCNA modules and associated data
#'
#' Make WGCNA modules from gene expression data. Also outputs mean or eigenvalue module expression and DAVID formatted gene lists
#'
#' @param fit List object output by fit_modules( )
#' @param Rsq.min Numeric minimum R-squared for soft threshold selection. If set, sft.value is ignored
#' @param sft.value Numeric soft threshold. Set when minimum R-squared is ignored
#' @param minModuleSize Numeric minimum module size. Default is 20
#' @param maxBlockSize Integer giving maximum block size for module detection. Default is 500
#' @param deepSplit Integer value between 0 and 4. Provides a simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive
#' @param networkType Character string. Network type. One of "unsigned", "signed", "signed hybrid". Default is "signed"
#' @param TOMType Character string. One of "none", "unsigned", "signed", "signed Nowick", "unsigned 2", "signed 2", or "signed Nowick 2". If "none", adjacency will be used for clustering. See TOMsimilarityFromExpr for details. Default is "signed"
#' @param mods.mean Logical if include mean module expression in output
#' @param mods.eigen Logical if include eigenvalue module expression in output
#' @param david Logical if include DAVID formatted genes in modules in output
#' @param nThread Integer for number of threads to use
#'
#' @return List including:
#' \itemize{
#'   \item{genes} Character vector of genes used in module building
#'   \item{sft} Data frame with soft thresholding selected for module building. Includes power, minimum R-squared, and connectivity
#'   \item{top.plot} ggplot object of soft thresholding topology
#'   \item{connect.plot} ggplot object of soft thresholding connectivity
#'   \item{mods} Data frame of genes in modules
#'   \item{mods.mean} Data frame of mean module expression for each library
#'   \item{mods.eigen} Data frame of module eigenvalue expression for each library
#'   \item{david} DAVID formatted data frame of genes in modules
#' }
#' @export
#'
#' @examples
#' fit <- fit_modules(dat = example.voom)
#' dat.mods <- make_modules(fit = fit, sft.value = 14,
#'     mods.mean = TRUE, mods.eigen = TRUE, david = TRUE)

make_modules <- function(fit,
                         Rsq.min = NULL, sft.value = NULL,
                         minModuleSize = 20, maxBlockSize=500, deepSplit = 3,
                         networkType="signed", TOMType="signed",
                         mods.mean = FALSE, mods.eigen = FALSE, david = FALSE,
                         nThread=2){
  SFT.R.sq <- Power <- module <- geneName <- module.char <- NULL
  ##### Soft-thresholding #####
  #Select threshold
  if(!is.null(Rsq.min)){
    sft.select <- dplyr::filter(fit$sft, SFT.R.sq  >= Rsq.min)
    if(nrow(sft.select)>0){
      power.t <- min(sft.select$Power)
    } else {
      stop("R-squared minimum not reached. Please input lower Rsq.min or set sft.value instead.")
    }
  } else if(!is.null(sft.value)){
    sft.select <- dplyr::filter(fit$sft, Power == sft.value)
    power.t <- unique(sft.select$Power)
  } else{
    stop("Please set Rsq.min or sft.value.")
  }

  ##### Build modules #####
  cor <- WGCNA::cor #Reassign cor fxn to prevent namespace error with stats
  mod.net <- WGCNA::blockwiseModules(t(fit$dat$E),
                                     power=power.t,
                                     networkType=networkType,
                                     TOMType=TOMType,
                                     maxBlockSize=maxBlockSize,
                                     minModuleSize=minModuleSize,
                                     deepSplit=deepSplit,
                                     numericLabels=TRUE,
                                     saveTOMFileBase="TOM-blockwise",
                                     nthreads=nThread)
  cor <- stats::cor #Revert cor fxn assignment

  mods <- as.data.frame(mod.net$colors) %>%
    tibble::rownames_to_column("geneName") %>%
    dplyr::rename(module = "mod.net$colors") %>%
    #add leading 0 to module names for correct sorting of factor
    dplyr::mutate(module.char = ifelse(module <= 9,
                                       paste("0", module, sep=""),
                                       module)) %>%
    #Add color var
    dplyr::mutate(mod.color = WGCNA::labels2colors(mod.net$colors))

  #Add gene info if present
  if(!is.null(fit$dat$genes)){
    mods <- mods %>%
      dplyr::left_join(fit$dat$genes, by="geneName")
  }

  ##### Combine results into list #####
  dat.mods <- list()
  dat.mods[["genes"]] <- fit$genes
  dat.mods[["mods"]] <- mods
  dat.mods[["sft"]] <- sft.select

  ##### plot #####
  dat.mods[["top.plot"]] <- fit$top.plot +
    ggplot2::geom_hline(ggplot2::aes(yintercept = sft.select$SFT.R.sq[1]), color="red") +
    ggplot2::geom_vline(xintercept = power.t, color="red")
  dat.mods[["connect.plot"]] <- fit$connect.plot +
    ggplot2::geom_hline(ggplot2::aes(yintercept = sft.select$mean.k.[1]), color="red") +
    ggplot2::geom_vline(xintercept = power.t, color="red")

  ##### Mean module expression #####
  if(mods.mean){
    mods.voom <- mods %>%
      #Combine count and module data
      dplyr::select(geneName, module.char) %>%
      dplyr::left_join(tibble::rownames_to_column(as.data.frame(fit$dat$E), "geneName"),
                       by="geneName") %>%
      #Calculate mean by module
      dplyr::group_by(module.char) %>%
      dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>%
      tibble::rownames_to_column() %>%
      #Make module names
      dplyr::mutate(rowname=paste("module", module.char, sep="_")) %>%
      tibble::column_to_rownames()

    dat.mods[["mods.mean"]] <- mods.voom
  }

  ##### Eigenvalue expression #####
  if(mods.eigen){
    mods.E <- mod.net$MEs %>%
      t() %>% as.data.frame() %>%
      tibble::rownames_to_column("module.char") %>%
      dplyr::mutate(module.char = as.numeric(gsub("ME","", module.char))) %>%
      dplyr::mutate(module.char = ifelse(module.char <= 9,
                                         paste("0", module.char, sep=""),
                                         module.char)) %>%
      dplyr::arrange(module.char)
    rownames(mods.E) <- paste0("module_", mods.E$module.char)
    dat.mods[["mods.eigen"]] <- mods.E
  }

  ##### DAVID format #####
  if(david){
    #list modules
    mod.names <- sort(rownames(mods.voom))
    #calculate maximum module gene list length
    max.mod <- max(table(mod.net$colors))
    #Create data frame with that number of rows
    david.df <- data.frame(rowname=1:max.mod)

    for (i in 1:length(mod.names)){
      #Create module name as "module#"
      mod.name <- mod.names[i]
      #Filter gene list to module of interest
      gene.list <- mods %>%
        dplyr::filter(module == i-1) %>%
        dplyr::select(geneName)
      #Calculate the total number of rows that need to be added for nrows to match
      add.genes <- max.mod - nrow(gene.list)
      #Add NAs to fill out gene.list and convert to data frame for merging
      gene.list <- c(gene.list$geneName, rep(NA, times=add.genes))
      gene.list <- as.data.frame(gene.list)
      colnames(gene.list) <- mod.name

      #Combine module lists and rename to mod.name
      david.df <- david.df %>%
        dplyr::bind_cols(gene.list)
    }

    dat.mods[["david"]] <- david.df
  }

  return(dat.mods)
}
