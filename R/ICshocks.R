#' @title Information criterion for determining the number of shocks in a factor model
#' @description The function gives the number of shocks that minimizes the information criterion.
#' @param x a dataset;
#' @param delta a real number within the range (0,1/2) for the sensitivity of the tolerance level to the size of the dataset;
#' @param m a finite positive real number defining the tolerance level;
#' @param r a positive integer corresponding to the number of factors;
#' @param p a positive integer corresponding to the number of lags to be considered within the model.
#' @return A \code{list} containing two elements:
#' 
#' \item{q_star}{The number of shocks minimizing the information criterion;}
#' \item{p}{The number of lags used.}
#' 
#' @references Bai, J., Ng, S. (2007). Determining the Number of Primitive Shocks in Factor Models. Journal of Business & Economic Statistics, 25(1), 52-60. <https://doi.org/10.1198/073500106000000413> 
#' @importFrom vars VARselect VAR
#' @import stats
#' @export

ICshocks<- function(x, delta = 0.1, m = 1, r = NULL, p = NULL){
  
  # discarting rows with missing values
  x <- na.omit(x)
  
  # Normalization of the database
  x <- as.matrix(x)
  Mx <- colMeans(x)
  Wx <- apply(x, MARGIN = 2, FUN = sd)
  for(i in 1:ncol(x)){
    x[,i] <- (x[,i] - Mx[i])/Wx[i]
  }
  
  # Other parameters: size of the database
  TT <- nrow(x)
  N <- ncol(x)
  
  # Checking parameters delta and m (if not specified by the user we have used the values 
  # from the Bai and Ng 2007 paper)
  if(delta <= 0 || delta >= 1/2){stop("Delta needs to be within the (0,1/2) interval")}
  if(m <= 0 || is.infinite(m) ){stop("m needs to be be within the (0, Inf) interval")}
  
  # if not specified we will use the ICP2 criterium from Bai and Ng (2002) as
  # done in Bai and Ng 2007 to estimate the number of factors
  if(is.null(r)){
    
    eigen <- eigen(cov(x))
    result <- c(1:20)*0
    
    for(i in 1:20){
      eigenvectors <- eigen$vectors[,1:i]
      factors <- x %*% eigenvectors
      V <- sum(diag(t(x - factors %*% t(eigenvectors)) %*% (x - factors %*% t(eigenvectors)))/(nrow(x)*ncol(x)))
      result[i] <- log(V)+i*(N+TT)/(N*TT)*log(min(N,TT))
    }
    
    r <- which.min(result)
  }
  
  # estimate the static factors using principal components
  eigen <- eigen(cov(x))
  eigenvectors <- eigen$vectors[,1:r]
  Factor <- x %*% eigenvectors
  
  # number of lags for the VAR
  if(is.null(p)){
    select <- vars::VARselect(y = Factor, lag.max = 12)
    p <- as.numeric(names(sort(table(select$selection), decreasing = TRUE)[1]))
  }
  
  # residuals from the VAR in F
  VAR_F <- suppressWarnings(vars::VAR(y = Factor, p = p, type = "const")) # suppressing warning on column names
  u <- NULL
  for(i in 1:length(VAR_F$varresult)){
   u <- cbind(u, VAR_F$varresult[[i]]$residuals)
  }
  uu <- stats::var(u)
  
  # eigenvalues of uu
  eigen_uu <- eigen(uu)
  
  # calculating the sum of the squared eigenvalues
  sum = 0
  for(i in 1:ncol(uu)){
    sum = sum + eigen_uu$values[i]^2
  }
  
  # calculating vector D1
  D <- c(1:r)*0
  for(i in 1:r-1){
    D[i] = (eigen_uu$values[i+1]^2/sum)^(1/2)
  }
  
  # calculating vector K3
  exponent <- 1/(2-delta)
  K3 <- D < (m/(min(N^exponent,TT^exponent)))
  
  # Output
  list(q_star = which(K3 == TRUE)[1], p = p)
  
}