#' @title Monthly to quarterly transformation
#' @description It transforms a monthly time series in a quarterly one, selecting the last month of the quarter to represent the value of the quarter.
#' @param x a \code{ts} or \code{mts} in monthly frequency
#' @param reference_month a vector to define the reference month that will represent the quarter. Default is 3. The options are 1, 2, 3 or 'mean'.
#' @return The correspondent quarterly transformation.
#' @examples 
#' # Selecting only last month of matrix time series BRGDP:
#' month2qtr(BRGDP)
#' 
#' # Vehicle production in the quarter from vehicle production in the month
#' month2qtr(stats::filter(BRGDP[,3],c(1,1,1),sides=1))
#' 
#' @import zoo
#' @importFrom lubridate year quarter
#' @importFrom xts xts apply.quarterly
#' @export

month2qtr <-  function (x, reference_month = 3){
  
  if(!reference_month %in% c(1,2,3,"mean")){
    stop("The reference_month should be either 1,2 or 3 or \"mean\"")
  }
  
  data <- zoo::as.Date(x)
  ano_inicial <- as.numeric(substr(data[1], 1, 4))
  meses <- substr(data, 6, 7)
  
  if(reference_month == "mean"){
    data2 <- xts::xts(x, order.by = data)
    data3 <- xts::apply.quarterly(data2, mean)
    x.tri <- ts(data3, 
                start = c(lubridate::year(index(data3))[1], lubridate::quarter(index(data3))[1]),
                frequency = 4)
  }else{
    
    months <- c(reference_month, 3 + reference_month, 6 + reference_month, 9 + reference_month)
    months <- paste0(0,months)
    months <- substr(months,nchar(months)-1,nchar(months))
  
    ultimo_tri <- meses %in% months
    tri <- which(months == meses[which(ultimo_tri)[1]])
    x.tri <- ts(data.frame(x)[ultimo_tri, ], 
              start = c(ano_inicial, tri), frequency = 4)
  }
  
  # output
  return(x.tri)
  
}