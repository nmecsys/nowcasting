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
#' @import zoo
#' @export

qtr2month<-function(x){
  data<-zoo::as.Date(x)
  datas<-seq(data[1],data[length(data)],by = 'months')
  out_x<-ts(rep(NA,length(datas)),start=  as.numeric(c(substr(as.Date(x)[1],1,4),substr(as.Date(x)[1],6,7))),frequency = 12)
  out_x[datas %in% data] <- x
  out_x<-stats::lag(out_x,-2)
  return(out_x)
}