#' @title Pseudo Real Time Data Base
#' @description Create a pseudo real time data base based on data and delays of disclosure stipulated by the user.
#' @param mts A \code{mts} with the series data.
#' @param delay A numeric vector with the delay in days the information is available after the reference month. Each element corresponds to the series in the respective column in \code{mts}. 
#' @param vintage The day when the data is supposed to be collected.
#' @return A \code{mts} with the series transformed.
#' @examples 
#' # Pseudo Real Time Data Base from data base BRGDP
#' PRTDB(mts = BRGDP, delay = c(1,30,60,90,20,10,30,60), vintage = "2017-10-01")
#' @import zoo lubridate
#' @export

PRTDB<-function(mts, delay, vintage = Sys.Date()){
  mts_new <- mts
  
  # define the last day of the month
  month_end<-as.Date(mts)+months(1)-days(1)
  
  # create list with release date
  release<-lapply(1:length(delay),function(x) month_end+days(delay)[x])
  
  # Eliminate information not available until the day
      for (i in 1:length(delay)){
      mts_new[release[[i]]>vintage,i] <- NA
      }
  
  mts_new <- ts(mts_new[!as.Date(mts_new)>vintage,],start=start(mts_new),frequency = frequency(mts_new))
  
  return(mts_new)
}