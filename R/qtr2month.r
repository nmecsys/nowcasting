#' @title Quarterly to monthly desagregation
#' @description It transforsms a quarterly time series in a monthly. To the last month of the monthly \code{ts} is set the value of the quarterly \code{ts}.
#' @param x Variable in quarterly frequency
#' @return The correpondent monthly transformation or agregation.
#' @examples 
#' # Selecting only last month of matrix time series BRGDP:
#' mestri_vintage<-month2qtr(BRGDP[,dim(BRGDP)[2]])
#' 
#' qtr2month(mestri_vintage) 
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