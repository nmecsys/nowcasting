#' @title Quarterly to monthly transformation
#' @description It transforms a quarterly  time series in a monthly one. The values of the quarterly \code{ts} are set to the last month of the quarter.
#' @param x a \code{ts} or \code{mts} in quarterly frequency
#' @return The correpondent monthly transformation.
#' @examples 
#' # Selecting the quarterly GDP variable in BRGDP
#' brgdp <- month2qtr(BRGDP[,ncol(BRGDP)])
#' 
#' qtr2month(brgdp) 
#' 
#' @importFrom zoo as.Date
#' @importFrom stats lag
#' @export

qtr2month <- function(x){
  data_q <- zoo::as.Date(x)
  data_m <- seq(data_q[1], data_q[length(data_q)], by = 'months')
  out_x <- ts(rep(NA,length(data_m)),
              start =  as.numeric(c(substr(data_q[1],1,4), substr(data_q[1],6,7))),
              frequency = 12)
  out_x[data_m %in% data_q] <- x
  out_x <- stats::lag(out_x, -2)
  
  # output
  return(out_x)
}