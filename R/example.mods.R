#' Example dataframe defining gene co-expression modules.
#'
#' @format A data frame with 1000 rows and 8 variables:
#' \describe{
#'   \item{geneName}{character. ENSEMBL gene ID}
#'   \item{module}{dbl. Module number defined by WGNCA}
#'   \item{module.char}{character. Module number as character defined by WGNCA}
#'   \item{mod.color}{character Module name as defined by WGNCA}
#'   \item{hgnc_symbol}{character. Current approved HGNC symbol.}
#'   \item{Previous symbols}{character. Previous HGNC symbols.}
#'   \item{Alias symbols}{character. Alias HGNC symbols.}
#'   \item{gene_biotype}{character. Gene product type. All = protein-coding.}
#'   }
#'
#' @source \url{https://github.com/altman-lab/P259_pDC_public}
#' @references Dill-McFarland et al. 2021. Eosinophil-mediated suppression and Anti-IL-5 enhancement of plasmacytoid dendritic cell interferon responses in asthma. J Allergy Clin Immunol. In revision
#' @description Example gene-coexpression modules defined in a dataframe.
#' @docType data
#' @name example.mods
#' @keywords datasets
"example.mods"
