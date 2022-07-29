#' Compare WGCNA module builds
#'
#' Make WGCNA modules from gene expression data with dynamic soft threshold selection. Also outputs mean module expression and DAVID formatted gene lists
#'
#' @param fit List object output by fit_modules( )
#' @param Rsq.min Numeric vector minimum R-squared for soft threshold selection. If set, sft.value is ignored
#' @param sft.value Integer vector soft threshold. Set when minimum R-squared is ignored
#' @param minModuleSize Numeric vector minimum module size. Default is 20
#' @param maxBlockSize Integer vector giving maximum block size for module detection. Default is 500
#' @param deepSplit Integer vector between 0 and 4. Provides a simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive
#' @param networkType Character vector. Network type from "unsigned", "signed", "signed hybrid". Default is "signed"
#' @param TOMType Character vector from "none", "unsigned", "signed", "signed Nowick", "unsigned 2", "signed 2", or "signed Nowick 2". If "none", adjacency will be used for clustering. See TOMsimilarityFromExpr for details. Default is "signed"
#' @param nThread Integer for number of threads to use
#'
#' @return List including:
#' \itemize{
#'   \item{summary} Data frame with total modules and sizes for all builds
#'   \item{modules} List of data frames with module sizes for each build
#' }#'
#' @export
#'
#' @examples
#' fit <- fit_modules(dat = example.voom)
#' compare.mods <- compare_modules(fit = fit, sft.value = c(14, 20))

compare_modules <- function(fit,
                            Rsq.min = NULL, sft.value = NULL,
                            minModuleSize = 20, maxBlockSize=500, deepSplit = 3,
                            networkType="signed", TOMType="signed",
                            nThread=2){
  ##### Check inputs #####
  if(!is.null(Rsq.min) & !is.null(sft.value)) {
    stop("Only one of Rsq.min or sft.value may be used.")
  }
  SFT.R.sq <- Power <- module <- param <- NULL

  ##### Loop through all unique parameter sets #####
  if(!is.null(Rsq.min)){
    # get sft thresholds for Rsq values
    sft.value <- c()
    for(r in Rsq.min) {
      sft.select <- dplyr::filter(fit$sft, SFT.R.sq  >= r)
      if(nrow(sft.select)>0){
        sft.value <- c(sft.value, min(sft.select$Power))
      }
    }

    sft.value <- unique(sft.value)
    if(length(sft.value) == 0){
      stop("R-squared minimum not reached. Please input lower Rsq.min or set sft.value instead.")
    }
  }
  #all unique parameter combinations
  params <- expand.grid(sft.value, minModuleSize, maxBlockSize, deepSplit,
                        networkType, TOMType)
  colnames(params) <- c("sft.value", "minModuleSize", "maxBlockSize", "deepSplit",
                        "networkType", "TOMType")

  #Add actual R-square value
  Rsq.actual <- fit$sft %>%
    dplyr::select(Power, SFT.R.sq) %>%
    dplyr::rename(sft.value = Power)

  params <- dplyr::left_join(params, Rsq.actual, by = "sft.value") %>%
    dplyr::select(sft.value, SFT.R.sq, dplyr::everything())

  message(paste("Running", nrow(params), "module builds."))

  ##### Build modules #####
  compare.result <- data.frame()
  mod.ls <- list()
  # for(i in 1:nrow(params)){
  for(i in 1:2){
    param.OI <- params[i,]

    ##### Build modules #####
    cor <- WGCNA::cor #Reassign cor fxn to prevent namespace error with stats
    mod.temp <- WGCNA::blockwiseModules(t(fit$dat$E),
                                        power=param.OI$sft.value,
                                        networkType=param.OI$networkType,
                                        TOMType=param.OI$TOMType,
                                        maxBlockSize=param.OI$maxBlockSize,
                                        minModuleSize=param.OI$minModuleSize,
                                        deepSplit=param.OI$deepSplit,
                                        numericLabels=TRUE,
                                        saveTOMFileBase="TOM-blockwise",
                                        nthreads=nThread)
    cor <- stats::cor #Revert cor fxn assignment

    ##### Extract results #####
    mods <- as.data.frame(mod.temp$colors) %>%
      dplyr::rename(module = "mod.temp$colors") %>%
      dplyr::count(module)
    mod.ls[[paste0("param",i)]] <- mods

    compare.result <- param.OI %>%
      dplyr::bind_cols(data.frame(
        param = paste0("param",i),
        genes = length(mod.temp$colors),
        mod0 = mods[mods$module == 0,"n"],
        modMin = min(mods[mods$module != 0,"n"]),
        modMax = max(mods[mods$module != 0,"n"]),
        modN = max(mods[mods$module != 0,"module"]))) %>%
      dplyr::bind_rows(compare.result)
  }

  mod.ls.result <- list()
  mod.ls.result[["summary"]] <- compare.result %>%
    dplyr::arrange(param)
  mod.ls.result[["modules"]] <- mod.ls
  return(mod.ls.result)
}
