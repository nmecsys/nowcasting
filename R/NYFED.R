#' @name  NYFED
#' @title Example of replication files in Giannone et al. 2008
#' @docType data
#' @format  A \code{list} with 4 elements: 
#' 
#' \itemize{
#' \item \code{base} is a \code{mts} with 25 series and 385 observations. There are missing values;
#' \item \code{legend} is a \code{data.frame} with specifications of the series in NYFED$base;
#' \item \code{Time} is a \code{vector} of length 385 with dates;
#' \item \code{blocks} is a \code{matrix} showing the groups of variables.
#' }
#' @description partial dataset used to replicate the results in the New York Fed Staff Nowcasting Report . 
#' @usage NYFED
#' @source This dataset is available in the following url: \url{https://github.com/FRBNY-TimeSeriesAnalysis/Nowcasting}
NULL
# 