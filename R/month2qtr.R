#' @title Monthly to quarterly transformation
#' @description It transforms a monthly time series in a quarterly one, selecting the last month of the quarter to represent the value of the quarter.
#' @param x a \code{ts} or \code{mts} in monthly frequency
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

month2qtr <- function(x){
  data <- zoo::as.Date(x)
  ano_inicial <- as.numeric(substr(data[1],1,4))
  meses <- substr(data,6,7)
  ultimo_tri <- meses %in% c("03","06","09","12")
  tri <- which(c("03","06","09","12") == meses[which(ultimo_tri)[1]])
  x.tri <- ts(data.frame(x)[ultimo_tri,], start = c(ano_inicial,tri), frequency = 4)
  x.tri
}