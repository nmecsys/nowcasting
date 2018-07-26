#' @title Monthly to quarterly transformation
#' @description It transforms a monthly time series in a quarterly one, selecting the last month of the quarter to represent the value of the quarter.
#' @param x a \code{ts} or \code{mts} in monthly frequency
#' @param reference_month an integer to define the reference month that will represent the quarter. Default is 3. The options are 1, 2 or 3.
#' @return The correspondent quarterly transformation.
#' @examples 
#' # Selecting only last month of matrix time series BRGDP:
#' month2qtr(BRGDP)
#' 
#' # Vehicle production in the quarter from vehicle production in the month
#' month2qtr(stats::filter(BRGDP[,3],c(1,1,1),sides=1))
#' 
#' @import zoo
#' @export

month2qtr <-  function (x, reference_month = 3){
  
  if(!reference_month %in% c(1,2,3)){
    stop("The reference_month should be 1,2 or 3")
  }
  
  data <- zoo::as.Date(x)
  ano_inicial <- as.numeric(substr(data[1], 1, 4))
  meses <- substr(data, 6, 7)
  
  months <- c(reference_month,3+reference_month,6+reference_month,9+reference_month)
  list_months <- vector()
  for(i in 1:length(months)){
    if(months[i]<10)
      list_months[i] <- paste("0",toString(months[i]), sep = "")
    else
      list_months[i] <- toString(months[i])
  } 
  
  ultimo_tri <- meses %in% list_months
  tri <- which(list_months == meses[which(ultimo_tri)[1]])
  x.tri <- ts(data.frame(x)[ultimo_tri, ], start = c(ano_inicial, 
                                                     tri), frequency = 4)
  x.tri

}