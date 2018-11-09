# #' @title Impact
# #' @description Estimate the change between a previous forecasting and a new one.
# #' @param out.old an output from nowcast function (old).
# #' @param out.new another output from nowcast function (newer).
# #' @param Y.old an \code{numeric} forecasting of y variable using informations in \code{out.old}.
# #' @param Y.new another \code{numeric} forecasting of y variable (newer) using informations in \code{out.new}.
# #' @param period a \code{character} vector reporting the period to evaluate the impact. The vector must have one or two positions to indicate the period range: "yyyy-mm" or c("yyyy-mm","yyyy-mm").
# #' @return A \code{list} containing two elements:
# #' \item{impact}{the impact of each variable in the Y.new - Y.old change.}
# #' \item{change}{the difference between Y.new and Y.old.}
# #' @references COLOCAR REFERENCIAS
# #' @import lubridate zoo
# #' @export

impact <- function(out.old = NULL, out.new = NULL, Y.old = NULL, Y.new = NULL, period = NULL){
  
  
  if(is.null(out.old) |  is.null(out.new) | is.null(Y.old) |  is.null(Y.new) |  is.null(period)){
    stop("The arguments out.old, out.new, Y.old, Y.new, and period cannot be NULL.")
  }

  for(i in 1:length(period)){
    if(is.Date(as.Date(as.yearmon(period[i])) != TRUE)){
      stop("The argument period should have elements of the form 'yyyy-mm'.")
    }
  }
  
  if(length(period) == 1){
    begin <- as.Date(as.yearmon(period))
    end <- begin
  }else if(length(period) == 2){
    begin <- as.Date(as.yearmon(period[1]))
    end <- as.Date(as.yearmon(period[2]))
  } else{
    stop("The argument period should have 1 or 2 dates in the form 'yyyy-mm'.")
  }
  
  if(length(Y.old)!=1 | !is.numeric(Y.old)){
    stop("The argument Y.old should be a numeric of length one")
  }
  
  if(length(Y.new)!=1 | !is.numeric(Y.new)){
    stop("The argument Y.new should be a numeric of length one")
  }
  
  X_old <- out.old$xfcst
  X_new <- out.new$xfcst
  coefficients_reg <- out.old$reg$coefficients[-1]
  factor_loading <- tryCatch(out.old$factors$Lambda[,1:ncol(out.old$factors$dynamic_factors)],
                             error = function(e) out.old$factors$Lambda)
  
  rownames(factor_loading) <- colnames(X_old)
  
  # consider only common variables from X_new and X_old
  common_names <- c(colnames(X_old), colnames(X_new))[duplicated(c(colnames(X_old), colnames(X_new)))]
  
  # normalization of X
  s <- apply(X_old, MARGIN = 2, FUN = sd, na.rm = T)
  m <- apply(X_old, MARGIN = 2, FUN = mean, na.rm = T)
  
  for(i in 1:dim(X_new)[2]){
    X_new[,i] <- (X_new[,i]-m[i])/s[i] 
  }
  
  for(i in 1:dim(X_old)[2]){
    X_old[,i] <- (X_old[,i]-m[i])/s[i] 
  }
  
  # redefining the matrixes to compare
  X_new <- X_new[,common_names]
  X_old <- X_old[,common_names]
  factor_loading <- factor_loading[common_names,]
  
  # defining the period
  times <- as.Date(stats::time(X_old))
  begin_pos <- which(times==begin)
  end_pos <- which(times==end) 
  
  # defining what changed within the period
  X_new <- X_new[begin_pos:end_pos,] 
  X_old <- X_old[begin_pos:end_pos,] 
  X_diff <- X_new - X_old
  if(length(period) == 1){
    X_diff <- matrix(X_diff, nrow = 1)
  }
  colnames(X_diff) <- common_names  
  # identify variables that changed
  cd <- which(colSums(X_diff) != 0)
  
  # calculating the impact
  impact <- NULL
  for(i in cd){
    X_diff1 <- X_diff
    X_diff1[,-i] <- 0
    F_diff <- solve(t(factor_loading) %*% factor_loading) %*% (t(factor_loading) %*% t(X_diff1))
    impact <- c(impact, sum(coefficients_reg %*% F_diff))
  }
  
  total_impact <- sum(impact)
  impact <- impact * (Y.new - Y.old) / total_impact
  
  # output
  return(list(impact = data.frame(variable_name = common_names[cd], impact = impact), 
              change = Y.new - Y.old))
  
}

