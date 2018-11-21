#' @title Information criteria for determining the number of factors in a factors model
#' @description Minimizes the selected information criterion to determine the number of factors to be used in an approximate factor model.
#' @param x a dataset;
#' @param rmax a positive integer corresponding to the maximum number of factors for which the information criterion should be tested;
#' @param type a positive integer corresponding to the chosen information criterion (1,2,3). The default is 2.
#' @return A \code{list} containing two elements:
#' 
#' \item{r_star}{The number of factors minimizing the information criterion;}
#' \item{IC}{A vector of values of the information criterion for the number of factors within the selected range.}
#' 
#' @references Bai, J., Ng, S. (2002). Determining the Number of Factors in Approximate Factor Models. Econometrica, 70(1), 191-221. <doi:10.1111/1468-0262.00273>
#' @export

ICfactors <- function(x, rmax = NULL, type = 2){
  
  # discarting rows with missing values
  x <- na.omit(x)
  
  # defining rmax and checking if it is a positive integer
  if(is.null(rmax)){rmax = min(dim(x)[2],20)}
  else if(rmax < 1 || rmax != as.integer(rmax)){stop("rmax needs to be a positive integer")}
  else if(rmax > dim(x)[2]){rmax = dim(x)[2]}
  
  # checking if the type is correctly specified
  if((type %in% c(1,2,3)) == F){stop("The information criterium type must be either 1, 2, or 3")}
  
  # Normalizing the database
  x <- as.matrix(x)
  Mx <- colMeans(x)
  Wx <- apply(x, MARGIN = 2, FUN = sd)
  
  for(i in 1:ncol(x)){
    x[,i] <- (x[,i] - Mx[i])/Wx[i]
  }
  
  # Determining the size of the database
  TT <- nrow(x)
  N <- ncol(x)
  
  if(type == 1){ 
    
    # Calculating the IC
    eigen <- eigen(cov(x))
    
    result <- c(1:rmax)*0
    
    for(r in 1:rmax){
      
      v <- eigen$vectors[,1:r]
      
      # Common factors
      factors <- x %*% v
      
      # Sum squared errors
      V <- sum(diag(t(x - factors %*% t(v)) %*% (x - factors %*% t(v)))/(nrow(x)*ncol(x)))
      
      result[r] <- log(V)+r*(N+TT)/(N*TT)*log(N*TT/(N+TT))
    }
    
    # plot
    graphics::plot(result, main = "ICR1", xlab = "Number of fators", ylab = "Index")
    graphics::points(which.min(result), result[which.min(result)], pch = 19, col ="red")
    
    # output  
    list(r_star = which.min(result)[1], IC = result)
    
  }else if(type == 2){
    
    # Calculating the IC
    eigen <- eigen(cov(x))
    
    result <- c(1:rmax)*0
    
    for(r in 1:rmax){
      
      v <- eigen$vectors[,1:r]
      
      # Common factors
      factors <- x %*% v
      
      # Sum squared errors
      V <- sum(diag(t(x - factors %*% t(v)) %*% (x - factors %*% t(v)))/(nrow(x)*ncol(x)))
      
      result[r] <- log(V)+r*(N+TT)/(N*TT)*log(min(N,TT))
    }
    
    #plot
    graphics::plot(result, main = "ICR2", xlab = "Number of fators", ylab = "Index")
    graphics::points(which.min(result), result[which.min(result)], pch = 19, col ="red")
    
    # output  
    list(r_star = which.min(result)[1], IC = result)
    
  }else if(type == 3){
   
    # Calculating the IC
    eigen <- eigen(cov(x))
    
    result <- c(1:rmax)*0
    
    for(r in 1:rmax){
      
      v <- eigen$vectors[,1:r]
      
      # Common factors
      factors <- x %*% v
      
      # Sum squared errors
      V <- sum(diag(t(x - factors %*% t(v)) %*% (x - factors %*% t(v)))/(nrow(x)*ncol(x)))
      
      result[r] <- log(V)+r*(log(min(N,TT))/min(N,TT))
    }
    
    #plot
    graphics::plot(result, main = "ICR3", xlab = "Number of fators", ylab = "Index")
    graphics::points(which.min(result), result[which.min(result)], pch = 19, col ="red")
    
    # output  
    list(r_star = which.min(result)[1], IC = result) 
    
  }
  
}