% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_metrics.R
\name{clean_metrics}
\alias{clean_metrics}
\title{Clean alignment metrics into a data frame}
\usage{
clean_metrics(dir = "data_raw/metrics/", flagstat = TRUE, picard = TRUE)
}
\arguments{
\item{dir}{Character. Path to directory containing metric files}

\item{flagstat}{Logical. Include samtools flagstat metrics?}

\item{picard}{Logical. Include Picard RNAseq metrics?}
}
\value{
Data frame with alignment quality metrics
}
\description{
Clean alignment metrics into a data frame
}
