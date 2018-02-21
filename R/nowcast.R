#' @title Nowcasting of a quarterly time series using a dynamic factor model.
#' @description Estimate nowcasting and forecasting for a quarterly series. For more details read the Vignettes.
#' @param y Stationary quarterly time series. 
#' @param x A monthly time series matrix (\code{mts}) representing regressors variables. The series must be stationary.
#' @param q Dynamic rank. Number of error terms.
#' @param r number of commom factors.
#' @param p AR order of factor model.
#' @param method There are three options: \code{"2sq"}: "Two stages: quarterly factors" as in Giannone et al. 2008; \code{"2sm"}: "Two stages: monthly factors" as in Bańbura and Runstler 2011; \code{"EM"}: Expected Maximization as in Bańbura et al. 2011.
#' @param blocks a binary matrix Nx3 that characterizes the regressors variables in global (1st column), nominal (2nd column) and real (3rd column). If \code{NULL}, the matrix assume 1 for all cells.
#' @return A \code{list} containing two elements:
#' 
#' \item{yfcst}{the original \code{y} series and its in-sample and out-of-sample estimations.}
#' \item{reg}{regression model between \code{y} and the estimated factors. Not available for EM method.}
#' \item{factors}{the estimated factors and DFM model coefficients.}
#' \item{xfcst}{the original regressors and their out-of-sample estimations.}
#' \item{month_y}{the monthly measure for quarterly \code{y} variable. Only available for EM method.}
#' 
#' @references Giannone, D., Reichlin, L., & Small, D. (2008). Nowcasting: The real-time informational content of macroeconomic data. Journal of Monetary Economics, 55(4), 665-676.<doi:10.1016/j.jmoneco.2008.05.010>
#' 
#' Bańbura, M., & Rünstler, G. (2011). A look into the factor model black box: publication lags and the role of hard and soft data in forecasting GDP. International Journal of Forecasting, 27(2), 333-346. <doi:10.1016/j.ijforecast.2010.01.011>
#' 
#' Bańbura M., Giannone, D. & Reichlin, L. (2011). Nowcasting, in Michael P. Clements and David F. Hendry, editors, Oxford Handbook on Economic Forecasting, pages 193-224, January 2011. <doi:10.1093/oxfordhb/9780195398649.001.0001>
#' 
#' @examples
#' \dontrun{
#' ### Method 2sq (two stages: quarterly factors)
#' gdp <- month2qtr(x = USGDP$base[,"RGDPGR"])
#' gdp_position <- which(colnames(USGDP$base) == "RGDPGR")
#' base <- Bpanel(base = USGDP$base[,-gdp_position],
#'                trans = USGDP$legend$Transformation[-gdp_position],
#'                aggregate = TRUE)
#' now2sq <- nowcast(y = gdp, x = base, r = 2, p = 2, q = 2, method = '2sq')
#'
#' ### Method 2sm (two stages: monthly factors)
#' base <- Bpanel(base = USGDP$base[,-gdp_position],
#'                trans = USGDP$legend$Transformation[-gdp_position],
#'                aggregate = F)
#' now2sm <- nowcast(y = gdp, x = base, r = 2, p = 2, q = 2, method = '2sm')
#'
#' ### Method EM
#' # selecting and transforming y  
#' gdp <- month2qtr(x = USGDPshort$base[,"GDPUS"])
#' gdp <- ts(c(gdp,NA,NA,NA,NA), start = start(gdp), frequency = 4)
#' gdp_stationary <- gdp/lag(gdp, k = -1) -1
#' gdp_position <- which(colnames(USGDPshort$base) == "GDPUS")
#' 
#' # selecting and transforming x 
#' base <- USGDPshort$base[,-gdp_position]
#' trans <- USGDPshort$legend[-gdp_position,"transformation"]
#' stationaryBase <- cbind(base[,trans == 1]/lag(base[,trans == 1], k = -1) - 1, diff(base[,trans == 2]))
#' colnames(stationaryBase) <- colnames(base)[c(which(trans == 1),which(trans == 2)) ]
#' stationaryBase <- stationaryBase[,colnames(base)]
#' 
#' # DFM estimation via EM
#' blocks <- matrix(c(1,0,1,1,0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,
#'                    0,1,1,0,1,1,0,1,0,1,1,1,0,1,1,0,1,0,1,1,1,0,1,1,0,1,
#'                    1,0,1,0,1,1,0,1,1,0,1,1,0,1,1,1,0,1,0,1,1,0,1,1,1,0), byrow = T, ncol = 3)
#' nowEM <- nowcast(y = gdp_stationary, x = stationaryBase, r = 1, p = 1, q = 1, method = 'EM', blocks = blocks)
#' }
#' @seealso \code{\link[nowcasting]{base_extraction}}
#' @export

nowcast <- function(y, x, q = NULL, r = NULL, p = NULL, method='2sq', blocks = NULL){

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
    res$reg$call$data <- NULL

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
    res$reg$call$data <- NULL
    
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
    
    factors<-list(
      dynamic_factors = ts(Res$FF[,c(1:r, (r*5 + 1):(r*5 + r), (2*r*5 + 1):(2*r*5 + r))], start = start(x), frequency = 12),
      T = Res$A,
      Z = Res$C,
      #xBart = ts(Res$X_sm, start = start(x), freq = 12),
      #xt = ts(Res$X_sm[,-ncol(Res$X_sm)], start = start(x), freq = 12),
      muBar = Res$Mx,
      mu = Res$Mx[-length(Res$Mx)],
      sigma = Res$Wx,
      Q = Res$Q, R = Res$R, initx = Res$Z_0, initV = Res$V_0)
    
    colnames(factors$dynamic_factors) <- c(paste0("globalFactor",1:r),paste0("nominalFactor",1:r),paste0("realFactor",1:r))
    
    
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
    
    #names(factors) <- c("dynamic_factors", "T", "Z","Q","R","initx","initV")
    res <- list(yfcst = Y, factors = factors, xfcst = fore_x, month_y = month_y)
    
  }


  return(res)
  
}
