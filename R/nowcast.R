#' @title Nowcasting of a quarterly time series using a dynamic factor model.
#' @description Estimate nowcasting and forecasting models for quarterly or monthly time series. For more details read the Vignettes.
#' @param formula An object of class "formula": a symbolic description of the model to be fitted.
#' @param data A monthly time series matrix (\code{mts}) of stationary variables. 
#' @param r number of commom factors.
#' @param q Dynamic rank. Number of error terms.
#' @param p AR order of factor model.
#' @param method There are three options: \code{"2s"} (two stages without factors aggregation as in Giannone et al. 2008); \code{"2s_agg"} (two stages with factors aggregation); \code{"EM"} (Expected Maximization as in Bańbura et al. 2011).
#' @param blocks a matrix that defines the variables loaded into the factors.
#' @param frequency A vector of integers indicating the frequency of the variables: 4 for quarterly, 12 for monthly.
#' @return A \code{list} containing two elements:

#' \item{yfcst}{the original \code{y} series and its in-sample and out-of-sample estimations.}
#' \item{reg}{regression model between \code{y} and the estimated factors. Not available for EM method.}
#' \item{factors}{the estimated factors and DFM model coefficients.}
#' \item{xfcst}{the original regressors and their out-of-sample estimations.}
#' 
#' @references Giannone, D., Reichlin, L., & Small, D. (2008). Nowcasting: The real-time informational content of macroeconomic data. Journal of Monetary Economics, 55(4), 665-676.<doi:10.1016/j.jmoneco.2008.05.010>
#' 
#' Bańbura, M., & Rünstler, G. (2011). A look into the factor model black box: publication lags and the role of hard and soft data in forecasting GDP. International Journal of Forecasting, 27(2), 333-346. <doi:10.1016/j.ijforecast.2010.01.011>
#' 
#' Bańbura M., Giannone, D. & Reichlin, L. (2011). Nowcasting, in Michael P. Clements and David F. Hendry, editors, Oxford Handbook on Economic Forecasting, pages 193-224, January 2011. <doi:10.1093/oxfordhb/9780195398649.001.0001>
#' 
#' @examples
#' \dontrun{
#' ### Method 2s (Using the Mariano and Murasawa aggregation method on the variables)
#' data(USGDP)
#' gdp_position <- which(colnames(USGDP$base) == "RGDPGR")
#' base <- Bpanel(base = USGDP$base[,-gdp_position],
#'                trans = USGDP$legend$Transformation[-gdp_position],
#'                aggregate = TRUE)
#' data <- cbind(USGDP$base[,"RGDPGR"], base)
#' colnames(data) <- c("RGDPGR", colnames(base))
#' frequency <- c(4, rep(12, ncol(data) -1))
#' now2s <- nowcast(formula = RGDPGR ~ ., data = data, r = 2, p = 2, q = 2,
#'                  method = '2s', frequency = frequency)
#'
#'
#' ### Method 2s_agg (Using the Mariano and Murasawa aggregation method on the factors)
#' data <- Bpanel(base = USGDP$base,
#'                trans = USGDP$legend$Transformation,
#'                aggregate = FALSE)
#' frequency <- c(rep(12, ncol(data) -1), 4)
#' now2s_agg <- nowcast(formula = RGDPGR ~ ., data = data, r = 2, p = 2, q = 2, 
#'                      method = '2s_agg', frequency = frequency)
#' 
#'
#' ### Method EM
#' # Replication of the NY FED nowcast
#' data(NYFED)
#' base <- NYFED$base
#' blocks <- NYFED$blocks$blocks
#' trans <- NYFED$legend$Transformation
#' frequency <- NYFED$legend$Frequency
#' data <- Bpanel(base = base, trans = trans, NA.replace = F, na.prop = 1)
#' nowEM <- nowcast(formula = GDPC1 ~ ., data = data, r = 1, p = 1, 
#'                  method = "EM", blocks = blocks, frequency = frequency)

#' 
#' }
#' @seealso \code{\link[nowcasting]{base_extraction}}
#' @export

nowcast <- function(formula, data, r = NULL, q = NULL, p = NULL, method = 'EM', blocks = NULL, frequency = NULL){
  
  # Checking user inputs
    
    # check formula
    if(is.character(formula)){
      formula <- as.formula(formula)
    }
    
    # the number of factors, shocks, and lags
    if(is.null(q) | is.null(r) | is.null(p)){
      warnings('Parameters q, r and p must be specified.')
    }
    
    # the frequencies of the variables
    if(length(frequency) != ncol(data)){
      stop("the length of the frequency vector must be the same as the number of variables in the data object")
    }
  
    if(sum(!frequency%in% c(12,4))!=0){
      stop("The frequencies should be a vector of numerics taking values 4 (quarterly) or 12 (monthly)")
    }
  
  # preparing the data
  k <- model.frame(formula, data, na.action = NULL)
  x <- ts(k[,-1], start = start(data), frequency = 12)
  
  y_position <- which(colnames(data) == colnames(k)[1])
  freq_y <- frequency[y_position]
  
  if(freq_y == 4){
    y <- month2qtr(ts(k[,1], start = start(data), frequency = 12))
  }else{
    y <- ts(k[,1], start = start(data), frequency = 12)
  }  
  
  
  # selecting the method
  if(method == '2s'){
    factors <- FactorExtraction(x, q = q, r = r, p = p)
    fatores <- factors$dynamic_factors
    prev <- bridge(y,fatores,freq_y)
    
    # undo normalization 
    fit <- as.matrix(factors$dynamic_factors) %*% t(factors$eigen$vectors[,1:r])
    colnames(fit) <- colnames(x)
    s <- apply(x, MARGIN = 2, FUN = sd, na.rm = T)
    M <- apply(x, MARGIN = 2, FUN = mean, na.rm = T)
    
    x1 <- fit
    fore_x <- x[,colnames(x) %in% colnames(fit)]
    for(i in colnames(fit)){
      x1[,i] <- s[i] * fit[,i] + M[i]
      fore_x[is.na(fore_x[,i]), i] <- x1[is.na(fore_x[,i]), i]
    }
    
    names(factors) <- c("dynamic_factors", "A", "Lambda","BB","Psi","initx","initV","eigen","std","mean")
    res <- list(yfcst = prev$main, reg = prev$reg, factors = factors, xfcst = fore_x)

  }else if(method == '2s_agg'){
    factors <- FactorExtraction(x, q = q, r = r, p = p)
    
    fatores <- stats::filter(factors$dynamic_factors, c(1,2,3,2,1), sides = 1)
    prev <- bridge(y,fatores,freq_y)
    
    aux_fator_month <- cbind(rep(1/9, length(zoo::as.Date(factors$dynamic_factors))),
                             factors$dynamic_factors)
    
    # undo normalization 
    fit <- as.matrix(factors$dynamic_factors) %*% t(factors$eigen$vectors[,1:r])
    colnames(fit) <- colnames(x)
    s <- apply(x, MARGIN = 2, FUN = sd,na.rm=T)
    M <- apply(x, MARGIN = 2, FUN = mean,na.rm=T)
    
    x1 <- fit
    fore_x <- x[,colnames(x) %in% colnames(fit)]
    for(i in colnames(fit)){
      x1[,i] <- s[i] * fit[,i] + M[i]
      fore_x[is.na(fore_x[,i]), i] <- x1[is.na(fore_x[,i]), i]
    }
    
    names(factors) <- c("dynamic_factors", "A", "Lambda","BB","Psi","initx","initV","eigen","std","mean")
    res <- list(yfcst = prev$main, reg = prev$reg, factors = factors, xfcst = fore_x)
    
  }else if(method == 'EM'){
    
    # checking validity of inputs
    if(p > 5){
      stop('Parameter p must be less or equal to 5.')
    }
    
    if(is.null(blocks)){
      stop("The block structure determining which variables load into which factors should be specified.")
    }
    
    if(!is.null(q)){
      message("Obs: for this estimation method the number of common shocks is assumed to be equal to the number of factors, i.e. q = r.")
      }

    # rewrite blocks as matrix
    if(!is.matrix(blocks)){blocks <- as.matrix(blocks)}
    
    # determine the number of blocks
    n_blocks <- dim(blocks)[2]
    
    x <- ts(model.frame(formula, data, na.action = NULL), start = start(data), frequency = frequency(data))
    
    # new frequency
    new_frequency <- c(frequency[y_position], frequency[-y_position])
    blocks <- rbind(blocks[y_position,],blocks[-y_position,])
    # determine the number of quarterly series
    nQ <- sum(new_frequency==4)
    
    # preparing X
    
      # 1) all quarterly series should be positioned at the last columns
      
        # reshuffle vector
        idx <- cumsum(rep(1,dim(x)[2]))
        V_Q <- which(new_frequency==4)
        if(is.integer(V_Q) && length(V_Q) == 0){
          idx_new <- idx
        }else{
          idx_M <- idx[-V_Q]
          idx_new <- c(idx_M,V_Q)
        }
        
        # adapting data base
        x <- x[,idx_new]
        
        # adapting blocks
        blocks <- as.matrix(blocks[idx_new,])
        
      # 2) position of the target variable
      y_pos <- which(colnames(x)==colnames(k)[1])
    
    #   
    Par <- list(r = rep(r,n_blocks), # Number of common factors
                p = p, # Number of lags in autoregressive of factor (same for all factors)
                max_iter = 500, # max number of itereations for the EM loop
                i_idio = c(rep(T,dim(x)[2]-nQ), rep(F,nQ)),
                Rconstr = matrix(c(
                  c(2,3,2,1),
                  c(-1,0,0,0),
                  c(0,-1,0,0),
                  c(0,0,-1,0),
                  c(0,0,0,-1))
                  ,4,5), 
                q = matrix(rep(0,4),4,1), 
                nQ = nQ, # Number of quarterly series
                blocks = blocks # Block loadings
    )
    
    Res <- EM_DFM_SS_block_idioQARMA_restrMQ(x,Par)
    
    # recovering the factors
    idx_factor <- c(1:r)
    if(dim(blocks)[2]>1){
      idx_factor_aux <- lapply(X = seq(1,dim(blocks)[2]-1), FUN = function(x){(x*r*5 + 1):(x*r*5 + r)})
      for(j in 1:length(idx_factor_aux)){idx_factor <- append(idx_factor, idx_factor_aux[[j]])}
    }
    
    # Factors and estimated parameters
    factors <- list(dynamic_factors = ts(Res$FF[,idx_factor], start = start(x), frequency = 12))
    colnames(factors$dynamic_factors) <- as.vector(sapply(X = 1:dim(blocks)[2],FUN = function(X){paste0("Block",X,"_factor",1:r)}))
    
    fore_x <- ts(Res$X_sm, start = start(x), frequency = 12)
    colnames(fore_x) <- colnames(x)
    
    # y monthly
    if(new_frequency[idx_new[y_pos]]==12){
      yprev <- ts(Res$X_sm[,y_pos], start = start(x), frequency = 12)
      y <- x[,y_pos]
    }
    
    # y quarterly
    if(new_frequency[idx_new[y_pos]]==4){
      yprev <- month2qtr(ts(Res$X_sm[,y_pos], start = start(x), frequency = 12))
      y <- month2qtr(x[,y_pos])
    }
    
    # Observed and forecast y
    Y <- cbind(y,yprev,yprev)
    Y[is.na(Y[,1]),2] <- NA
    Y[!is.na(Y[,1]),3] <- NA
    colnames(Y) <- c('y','in','out')
    
    # forecast X
    fore_x <- fore_x[,-y_pos]
    colnames(fore_x) <- colnames(x)[-y_pos]
    
    res <- list(yfcst = Y, 
                factors = factors, 
                xfcst = fore_x,
                Res = Res
    )
    
  }
  
  # output
  return(res)
  
}
