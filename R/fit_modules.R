#' Perform soft thresholding in preparation for WGCNA module building
#'
#' @param dat limma EList output by voom( ). Must include dat$E
#' @param genes Character vector of genes to used in module building. Must match rownames in dat. If not set, all genes in dat are used
#' @param powerVector Numeric vector of soft thresholding powers for which the scale free topology fit indices are to be calculated. Default c(1:30)
#' @param networkType Character string of network type. Allowed values are "unsigned", "signed", "signed hybrid"
#' @param nThread Integer for number of threads to use
#'
#' @return List including:
#' \itemize{
#'   \item{genes} Character vector of genes used in module building
#'   \item{sft} Data frame with soft thresholding results including R-squared, slope, and k metrics
#'   \item{top.plot} ggplot object of soft thresholding topology
#'   \item{connect.plot} ggplot object of soft thresholding connectivity
#' }
#' @export
#'
#' @examples
#' fit <- fit_modules(dat = example.voom)

fit_modules <- function(dat, genes=NULL,
                         powerVector=c(1:30), networkType = "signed",
                         nThread=2){

  ##### Check data format #####
  if(!is.numeric(dat$E) & is.numeric(dat$E[,1])){
    stop("Gene names must be rownames or the first column of dat.")
  }

  slope <- SFT.R.sq <- fit <- mean.k. <- Power <- value <- rowname <- signed.R.sq <- NULL

  ##### Format expression data #####
  #Move rownames if present
  dat.format <- dat
  if(is.numeric(dat$E)){
    dat.format$E <- as.data.frame(dat.format$E) %>%
      tibble::rownames_to_column()
  } else{
    dat.format$E <- as.data.frame(dat.format$E)
    colnames(dat.format$E)[1] <- "rowname"
  }

  ##### Subset significant genes #####
  dat.signif <- dat.format
  #If gene names as rownames
  if(!is.null(genes)){
    if(any(!(genes %in% dat.format$E[,1]))){ stop("Subset genes are not present in dat. Check that the same format is used, ie HGNC, ENSEMBL, etc") }

    dat.signif$E <- as.data.frame(dat.signif$E) %>%
      dplyr::filter(rowname %in% genes) %>%
      tibble::column_to_rownames()
  } else {
    dat.signif$E <- as.data.frame(dat.signif$E) %>%
      tibble::column_to_rownames()
  }

  ##### Soft-thresholding #####
  WGCNA::allowWGCNAThreads(nThreads=nThread)
  #Calculate
  sft <- WGCNA::pickSoftThreshold(t(dat.signif$E),
                                  powerVector = powerVector,
                                  networkType = networkType,
                                  verbose=0)
  sft.format <- as.data.frame(sft$fitIndices) %>%
    dplyr::mutate(signed.R.sq = -sign(slope)*SFT.R.sq)

  #Plot
  minR <- round(min(sft.format$signed.R.sq), digits=1)
  maxR <- round(max(sft.format$signed.R.sq), digits=1)
  sft.plot1 <- sft.format %>%
    ggplot2::ggplot(ggplot2::aes(x=Power, y=signed.R.sq, label=Power)) +
    ggplot2::geom_text(size=3) +
    ggplot2::theme_classic() +
    ggplot2::labs(y="Scale free topology model fit,signed R^2",
                  x="Soft threshold (power)") +
    ggplot2::scale_y_continuous(breaks =
                                  round(seq(minR, maxR, by = 0.1), digits = 1))

  sft.plot2 <- sft.format %>%
    ggplot2::ggplot(ggplot2::aes(x=Power, y=mean.k., label=Power)) +
    ggplot2::geom_text(size=3) +
    ggplot2::theme_classic() +
    ggplot2::labs(y="Mean connectivity",
                  x="Soft threshold (power)")

  ##### Combine results into list #####
  #Put rownames back so works with make_modules

  fit_result <- list()
  fit_result[["dat"]] <- dat.signif
  fit_result[["genes"]] <- rownames(dat.signif$E)
  fit_result[["sft"]] <- sft.format
  fit_result[["top.plot"]] <- sft.plot1
  fit_result[["connect.plot"]] <- sft.plot2
  return(fit_result)
}
