#' @title Quarterly to monthly transformation
#' @description It transforms a quarterly  time series in a monthly one. The values of the quarterly \code{ts} are set to the last month of the quarter.
#' @param x a \code{ts} or \code{mts} in quarterly frequency
#' @param reference_month a integer to define the position of a quarter value in a quarter. Default is 3. The options are 1, 2 or 3.
#' @param interpolation logical. The NA values can be estimated by linear interpolation (approx function from stats package). Default is FALSE.
#' @return The correpondent monthly transformation.
#' @examples 
#' # Selecting the quarterly GDP variable in BRGDP
#' brgdp <- month2qtr(BRGDP[,ncol(BRGDP)])
#' 
#' qtr2month(brgdp) 
#' 
#' @importFrom zoo as.Date yearmon
#' @importFrom stats lag approx is.ts
#' @export

qtr2month <- function(x, reference_month = 3, interpolation = FALSE){
  
  if(!reference_month %in% c(1,2,3)){
    stop("The reference_month should be 1,2 or 3")
  }
  
  if(!is.ts(x)){
    stop("x should be a ts object")
  }
  
  if(!is.null(dim(x))){
    stop("x should be a single ts object")
  }
  
  data_q <- zoo::as.Date(x)
  data_m <- seq(data_q[1], data_q[length(data_q)], by = 'months')
  out_x <- ts(rep(NA,length(data_m)),
              start =  as.numeric(c(substr(data_q[1],1,4), substr(data_q[1],6,7))),
              frequency = 12)
  out_x[data_m %in% data_q] <- x
  
  if(reference_month %in% c(2,3)){
    out_x <- stats::lag(out_x, -(reference_month-1))
    data_q <- zoo::as.Date(out_x)
    data_m <- seq(data_q[1], data_q[length(data_q)], by = 'months')
  }
  
  if(interpolation){
    xout <- zoo::as.yearmon(data_m)
    out_x <- stats::approx(out_x, xout = xout, method = "linear")
    out_x <- ts(out_x$y, start = out_x$x[1], end = out_x$x[length(out_x$x)], frequency = 12)
  }
  
  # output
  return(out_x)
}