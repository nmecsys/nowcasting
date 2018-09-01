#' @title Balanced panel
#' @description This function transforms the original monthly time series to its stationary representation following the user specification. The time series with more than 1/3 missing, i.e. NAs, are deleted and the remaining are modified such that the missings and outliers are replaced by an approximated value.
#'
#' The missings and outliers are “corrected” following the same method available in the replication files of Giannone et al. 2008. Outliers are defined as observations that lie more than 4 IQR from the median. All missings and outliers are replaced by the median. A centered moving average of degree **k** is calculated, forming a new panel. Then the missings and outliers are replaced by their equivalent observations on this new panel. We've made an important modification on the outlier_correction function found in the above mentioned files: Here the median of an even-sized sample is calculated by the mean of the two most central values, rather than using the largest of those numbers. Because of this modification the results obtained with the original replication files in (USGDP) are slightly different from those found here.
#' 
#' @param base A \code{mts} with the series to be transformed. 
#' @param trans A \code{vector} where each coordinate is a code for the transformation of the correspondent coordinate in the \code{base} argument. The transformation is specified by codes, as follows:
#' \itemize{
#'   \item{trans = 0: the original serie is preserved;}
#'   
#'   \item{trans = 1: monthly rate of change
#'   
#'   \deqn{\frac{x_{i,t} - x_{i,t-1}}{x_{i,t-1}}}}
#'   
#'   \item{trans = 2: monthly difference
#'   
#'   \deqn{x_{i,t} - x_{i,t-1}}}
#'   
#'   \item{trans = 3: monthly difference in year-over-year rate of change
#'   
#'   \deqn{\frac{x_{i,t} - x_{i,t-12}}{x_{i,t-12}}  -  \frac{x_{i,t-1} - x_{i,t-13}}{x_{i,t-13}}}}
#'   
#'   \item{trans = 4: monthly difference in year difference
#'   
#'   \deqn{(x_{i,t} - x_{i,t-12})  -  (x_{i,t-1} - x_{i,t-13})}}
#'   
#'   \item{trans = 5: yearly difference
#'   
#'   \deqn{(x_{i,t} - x_{i,t-12})}}
#'   
#'   \item{trans = 6: yearly rate of change
#'   
#'   \deqn{\frac{x_{i,t} - x_{i,t-12}}{x_{i,t-12}}}}
#'  } 
#' 
#' @param aggregate A \code{boolean} representing if you want aggregate the monthly variables to represent quarterly quantities. If \code{TRUE} the aggregation is made following the approximation of \emph{Mariano and Murasawsa 2003}.
#' @param k.ma A \code{numeric} representing the degree of the moving average correction.
#' @param na.prop A \code{numeric} representing the proportion of NA allowed. Default is 1/3.
#' @param h A \code{numeric} representing the number of steps ahead to forecasting. Default is 12.
#' @references Giannone, D., Reichlin, L., & Small, D. (2008). Nowcasting: The real-time informational content of macroeconomic data. Journal of Monetary Economics, 55(4), 665-676.<doi:10.1016/j.jmoneco.2008.05.010>
#' 
#' Mariano, R. S., & Murasawa, Y. (2003). A new coincident index of business cycles based on monthly and quarterly series. Journal of applied Econometrics, 18(4), 427-443.<doi:10.1002/jae.695>
#' @examples 
#' # Example from database BRGDP:
#' Bpanel(BRGDP,rep(3,ncol(BRGDP)))
#' @import zoo
#' @importFrom stats filter lag
#' @export


Bpanel <- function(base = NULL, trans = NULL, aggregate = F, k.ma = 3, na.prop = 1/3, h = 12){
  
  if(is.null(trans)){
    stop('trans can not to be NULL')
  }
  
  if(sum(is.na(trans)) != 0){
    stop('trans can not support missings values')
  }
  
  if(length(trans) != ncol(base)){
    stop('the number of elements in the vector must be equal to the number of columns of base')
  }
  
  if(sum(!names(table(trans)) %in% c(0:6)) != 0){
    stop('the only available transformations are 0, 1, 2, 3, 4, 5 and 6.')
  }
  
  if(na.prop <= 0 | na.prop >= 1){
    stop("na.prop must be between 0 and 1.")
  }
  
  # data transformation
  base1 <- base
  for(j in 1:ncol(base)){
    base1[,j] <- NA
    if(trans[j] == 1){  # monthly rate of change
      temp <- diff(base[,j]) / stats::lag(base[,j], -1)
      base1[-1,j] <- temp
    }else if(trans[j] == 2){ # monthly difference
      temp <- diff(base[,j])
      base1[-1,j] <- temp
    }else if(trans[j] == 3){ # monthly difference in year-over-year rate of change
      temp <- diff(diff(base[,j], 12) / stats::lag(base[,j], -12))
      base1[-c(1:13),j] <- temp
    }else if(trans[j] == 4){ # monthly difference in year difference
      temp <- diff(diff(base[,j],12))
      base1[-c(1:13),j] <- temp
    }else if(trans[j] == 5){ # yearly difference
      temp <- diff(base[,j],12)
      base1[-c(1:12),j] <- temp  
    }else if(trans[j] == 6){ # yearly rate of change
      temp <- base[,j] / stats::lag(base[,j],-12)
      base1[-c(1:12),j] <- temp  
    }else if(trans[j] == 0){ # no transformation
      base1[,j] <- base[,j]
    }
  }
  
  # quartly quantity transformation
  if(aggregate == T){
    base1 <- stats::filter(base1, c(1,2,3,2,1), sides = 1)
  }
  colnames(base1) <- colnames(base)
  
  # remove series whith more than na.prop missings values
  SerOk <- colSums(is.na(base1)) < (nrow(base1) * na.prop)
  base2 <- base1[, which(SerOk)]
  
  if(sum(SerOk) == 1){
    stop("the procedure can not be done with only one series available.")
  }
  
  if(sum(!SerOk) > 0){
    warning(paste(sum(!SerOk),'series ruled out due to lack in observations (more than', round(na.prop*100,2),'is NA).'))
  }
  
  seriesdeletadas <- colnames(base1[, which(!SerOk)])
  print(seriesdeletadas)
  
  
  # replacing missings and outliers 
  base3 <- base2 * NA
  
  for(i in 1:ncol(base2)){
    # ignoring the last missings values
    na <- is.na(base2[,i])
    na2 <- NULL
    for(j in 1:length(na)){
      na2[j] <- ifelse(sum(na[j:length(na)]) == length(j:length(na)), 1, 0)
    }
    na_position <- min(which(na2 == 1)) - 1
    if(length(which(na2 == 1)) == 0){ na_position <- nrow(base2)} 
    base3[,i] <- c(outliers_correction(base2[1:na_position,i], k.ma), rep(NA, nrow(base2) - na_position))
  }
  
  # add h lines to base
  base4 <- ts(rbind(base3, matrix(NA, nrow = h, ncol = ncol(base3))),
            start = start(base3), frequency = 12)
  
  # output
  return(base4)
}

