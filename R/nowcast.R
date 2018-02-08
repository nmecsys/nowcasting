#' @title Nowcasting of a quarterly time serie using a dynamic factor.
#' @description Estimate nowcasting and foreacasting for a quarterly time serie. For more details read the Vignettes.
#' @param y Stationary quarterly time-series 
#' @param x A time series matrix (\code{mts}) representing the regressors of interest. The series must be stationary.
#' @param q Dynamic rank. Number of error terms.
#' @param r Static rank or number of factors (r>=q, for methods 2sq and 2sm).
#' @param p AR order of factors.
#' @param method 2sq: Two stages quarterly as in Giannone et al. 2008; 2sm: Two stages monthly as in Bańbura and Runstler 2011; EM: Expected Maximization as in Bańbura et al. 2011
#' @param blocks only for EM method. Select which factors impact the variables (global, nominal or real).
#' @return A \code{list} containing two elements:
#' 
#' A \code{mts} named \code{main} contains the original serie, the estimation in the sample, the estimation out of the sample;
#' 
#' A \code{list} named \code{factors} contains the estimated factors and coeffients.
#' 
#' A \code{mts} named \code{fore_x} contains the output of all regressors.
#' 
#' A \code{mts} named \code{month_y} contains the a monthly measure for GDP. 
#' 
#' @references Giannone, D., Reichlin, L., & Small, D. (2008). Nowcasting: The real-time informational content of macroeconomic data. Journal of Monetary Economics, 55(4), 665-676.<doi:10.1016/j.jmoneco.2008.05.010>
#' 
#' Bańbura, M., & Rünstler, G. (2011). A look into the factor model black box: publication lags and the role of hard and soft data in forecasting GDP. International Journal of Forecasting, 27(2), 333-346. <doi:10.1016/j.ijforecast.2010.01.011>
#' 
#' Bańbura M., Giannone, D. & Reichlin, L. (2011). Nowcasting, in Michael P. Clements and David F. Hendry, editors, Oxford Handbook on Economic Forecasting, pages 193-224, January 2011. <doi:10.1093/oxfordhb/9780195398649.001.0001>
#' 
#' @examples
#' \dontrun{
#' # nowcast function examples:
#' ### Method 2sq
#' pib<-BRGDP[,8]
#' y<-month2qtr(diff(diff(pib,3),12))
#' x<-Bpanel(BRGDP[,-8],rep(4,dim(BRGDP)[2]),aggregate = T)
#' q<-1
#' r<-2
#' p<-1
#' now_2sq<-nowcast(y,x,q,r,p,method = '2sq')
#'
#' ### Method 2sm
#' pib<-BRGDP[,8]
#' y<-month2qtr(diff(diff(pib,3),12))
#' x<-Bpanel(BRGDP[,-8],rep(4,dim(BRGDP)[2]),aggregate = F)
#' now_2sm<-nowcast(y,x,q,r,p,method = '2sm')
#'
#' ### Method EM
#' y<-month2qtr(diff(diff(pib,3),12))
#' x<-Bpanel(BRGDP[,-8],rep(4,dim(BRGDP)[2]),aggregate = F)
#' now_em<-nowcast(y,x,q,r,p,'EM')
#' }
#' @seealso \code{\link[nowcasting]{base_extraction}}
#' @export

nowcast <- function(y, x, q = NULL, r = NULL, p = NULL,method='2sq',blocks = NULL){

  if(is.null(q) | is.null(r) | is.null(p)){
    warnings('Parameters q, r and p must be specified.')
  }
  
  if(method=='2sq'){
    factors <- FactorExtraction(x, q = q, r = r, p = p)
    fatores <- factors$dynamic_factors
    prev <- bridge(y,fatores)

    # voltar da padronização
    fit <- as.matrix(factors$dynamic_factors) %*% t(factors$eigen$vectors[,1:r])
    colnames(fit) <- colnames(x)
    s <- apply(x, MARGIN = 2, FUN = sd,na.rm=T)
    M <- apply(x, MARGIN = 2, FUN = mean,na.rm=T)
    # x <- x
    # z <- x
    # for(i in 1:dim(x)[2]){
    #   z[,i] <- (x[,i] - M[i])/s[i]
    # }
    x1<-fit
    fore_x<-x[,colnames(x) %in% colnames(fit)]
    for(i in colnames(fit)){
      x1[,i]<-s[i]*fit[,i]+M[i]
      fore_x[is.na(fore_x[,i]),i] <- x1[is.na(fore_x[,i]),i]
    }
 
    names(factors) <- c("dynamic_factors", "A", "Lambda","BB","Psi","initx","initV","eigen")
    res <-list(yfcst = prev$main, reg = prev$reg, factors = factors, xfcst = fore_x)
    
  }else if(method=='2sm'){
    factors <- FactorExtraction(x, q = q, r = r, p = p)
    fatores <- stats::filter(factors$dynamic_factors, c(1,2,3,2,1), sides = 1)
    prev <- bridge(y,fatores)
    
    # aux_month<-prev$reg$coefficients*cbind(rep(1,length(zoo::as.Date(factors$dynamic_factors))),factors$dynamic_factors)
    # month_y<-ts(rowSums(aux_month),start=start(factors$dynamic_factors),freq=12)
    aux_fator_month<-cbind(rep(1/9,length(zoo::as.Date(factors$dynamic_factors))),factors$dynamic_factors)
    month_y<-ts(aux_fator_month%*%prev$reg$coefficients,start=start(factors$dynamic_factors),frequency=12)
    
    # voltar da padronização
    fit <- as.matrix(factors$dynamic_factors) %*% t(factors$eigen$vectors[,1:r])
    colnames(fit) <- colnames(x)
    s <- apply(x, MARGIN = 2, FUN = sd,na.rm=T)
    M <- apply(x, MARGIN = 2, FUN = mean,na.rm=T)
    # x <- x
    # z <- x
    # for(i in 1:dim(x)[2]){
    #   z[,i] <- (x[,i] - M[i])/s[i]
    # }
    x1<-fit
    fore_x<-x[,colnames(x) %in% colnames(fit)]
    for(i in colnames(fit)){
      x1[,i]<-s[i]*fit[,i]+M[i]
      fore_x[is.na(fore_x[,i]),i] <- x1[is.na(fore_x[,i]),i]
    }
    
    names(factors) <- c("dynamic_factors", "A", "Lambda","BB","Psi","initx","initV","eigen")
    res <- list(yfcst = prev$main, reg = prev$reg, factors = factors, xfcst = fore_x, month_y = month_y)
    
    
  }else if(method=='EM'){
    
    if(p > 5){
      stop('Parameter p must be less than 5.')
    }
    
    # y1<-qtr2month(y)
    # y1[rep(which(!is.na(y1)),each=2)-c(2,1)]<-rep(y1[!is.na(y1)],each=2)
    # X<-cbind(x,y1)
    X<-cbind(x,qtr2month(y))
    if(is.null(blocks)){
      blocks<-matrix(rep(1,dim(X)[2]*3),dim(X)[2],3)
    }
    Par<-list(r=rep(r,3),p=p,max_iter=50,i_idio=c(rep(T,dim(x)[2]),F),
              Rconstr = matrix(c(
                c(2,3,2,1),
                c(-1,0,0,0),
                c(0,-1,0,0),
                c(0,0,-1,0),
                c(0,0,0,-1))
                ,4,5),
              q = matrix(rep(0,4),4,1),nQ = 1,
              blocks = blocks)
    Res<-EM_DFM_SS_block_idioQARMA_restrMQ(X,Par)
    
    factors<-list(dynamic_factors = ts(Res$FF[,c(1:r, (r*5 + 1):(r*5 + r), (2*r*5 + 1):(2*r*5 + r))], start = start(x), freq = 12),
                  A = Res$A, C = Res$C, Q = Res$Q, R = Res$R, initx = Res$Z_0,
                  initV = Res$V_0)
    colnames(factors$dynamic_factors) <- c(paste0("globalFactor",1:r),paste0("nominalFactor",1:r),paste0("realFactor",1:r))
    
    if(is.null(blocks)){
      factors$dynamic_factors <- factors$dynamic_factors[,1:r]
    }
    
    # fore_x<-ts(Res$X_sm[,-dim(Res$X_sm)[2]],start=start(X),frequency = 12)
    fore_x<-ts(Res$X_sm,start=start(X),frequency = 12)
    yprev<-month2qtr(ts(Res$X_sm[,dim(Res$X_sm)[2]],start=start(X),frequency = 12))
    Y<-cbind(y,yprev,yprev)
    Y[is.na(Y[,1]),2]<-NA
    Y[!is.na(Y[,1]),3]<-NA
    colnames(Y)<-c('y','in','out')
    
    ind<-c(1:r,1:r+r*5,1:r+r*5*2,dim(Res$C)[2]-4)
    month_y<-ts(Res$Mx[length(Res$Mx)]/9+Res$FF[,ind]%*%Res$C[nrow(Res$C),ind]*Res$Wx[length(Res$Wx)],start=start(X),frequency = 12)

    # Essa é uma medida trimestral do PIB acumulado nos últimos três meses
    # month_y<-ts(Res$X_sm[,dim(Res$X_sm)[2]],start=start(X),frequency = 12)
    fore_x = fore_x[,-dim(fore_x)[2]]
    colnames(fore_x)<-colnames(x)
    
    names(factors) <- c("dynamic_factors", "T", "Z","Q","R","initx","initV")
    res <- list(yfcst = Y,factors = factors, xfcst = fore_x, month_y = month_y)
    
  }


  return(res)
  
}
