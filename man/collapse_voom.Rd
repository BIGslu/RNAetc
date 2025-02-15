% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collapse_voom.R
\name{collapse_voom}
\alias{collapse_voom}
\title{Combine all parts of a voom object into a single dataframe}
\usage{
collapse_voom(
  dat,
  libraryID = "libID",
  geneID = "geneName",
  include_weights = FALSE
)
}
\arguments{
\item{dat}{EList voom object resulting from voom() or voomWithQualityWeights()}

\item{libraryID}{Character. Column name in dat$targets that matches column names in dat$E}

\item{geneID}{Character. Column name in dat$genes that matches row names in dat$E}

\item{include_weights}{Logical. If gene-level weights should be included in the result. Default is FALSE}
}
\value{
data frame
}
\description{
Combine all parts of a voom object into a single dataframe
}
\examples{
combined <- collapse_voom(example.voom)

}
