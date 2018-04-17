# #' @import R.matlab matlab zoo Matrix
# #' @export 


# Auxiliar functions ------------------------------------------------------

para_const <- function(X = NULL, P = NULL, lag = NULL){
  
  Z_0 <- P$Z_0
  V_0 <- P$V_0
  A <- P$A
  C <- P$C
  Q <- P$Q
  R <- P$R
  Mx <- P$Mx
  Wx <- P$Wx
  
  # %--------------------------------------------------------------------------
  # % Preparation of the data
  # %--------------------------------------------------------------------------
  TT <- nrow(X)
  N <- ncol(X)
  nomes<-names(X)
  
  #% Standardise x
  # xNaN <- (X-repmat(Mx,TT,1))/repmat(Wx,TT,1)
  xNaN<-(X-kronecker(Mx,rep(1,TT)))/kronecker(Wx,rep(1,TT))
  
  y <- t(xNaN)
  
  #%final run of the Kalman filter
  out <- runKF_lag(y, A, C, Q, R, Z_0, V_0, lag)
  # Zsmooth <- out$Zsmooth
  # P <- out$P
  Zsmooth <- out$xsmooth
  P <- out$Vsmooth
  
  Zsmooth <- t(Zsmooth)
  x_sm <- Zsmooth[2:nrow(Zsmooth),] %*% t(C)
  # X_sm <- matlab::repmat(Wx,TT,1) %*% x_sm + matlab::repmat(Mx,TT,1)
  X_sm <- kronecker(Wx,rep(1,TT)) * x_sm + kronecker(Mx,rep(1,TT))
  
  # %--------------------------------------------------------------------------
  # %   Loading the structure with the results
  # %--------------------------------------------------------------------------
  
  X_sm<-as.data.frame(X_sm)
  names(X_sm)<-nomes
  Res<-list(X_sm=X_sm,P=P)
  
  # Res$P <- P
  # Res$X_sm <- X_sm
  
  # output
  return(Res)
}

runKF_lag <- function(y = NULL, A = NULL, C = NULL, Q = NULL, R = NULL, x_0 = NULL, Sig_0 = NULL, k = NULL){
  
  n <- nrow(C)
  r <- ncol(C)
  
  if(k > 0){
    C <- cbind(C, zeros(n,k*r))
    # A <- bdiag(A,zeros(k*r))
    A <- magic::adiag(A,zeros(k*r))
    A[(r+1):nrow(A), 1:(k*r)] <- eye(k*r)
    # Q <- bdiag(Q,zeros(k*r))
    Q <- magic::adiag(Q,zeros(k*r))
    x <- zeros(k*r,1)
    # colnames(x) <- colnames(x_0)
    x_0 <- rbind(x_0,x)
    # Sig_0 <- bdiag(Sig_0,zeros(k*r,k*r))
    Sig_0 <- magic::adiag(Sig_0,zeros(k*r,k*r))
  }
  
  S <- SKF_lag(y,C,R,A,Q, x_0, Sig_0)
  S <- FIS(y,C,R,A,Q,S)
  
  xsmooth <- S$AmT[1:r,]
  Vsmooth <- S$PmT
  
  # output
  return(list(xsmooth = xsmooth, Vsmooth = Vsmooth))
  
}

SKF_lag <-function(Y,Z,R,TT,Q,A_0,P_0){
  # %______________________________________________________________________
  # % Kalman filter for stationary systems with time-varying system matrices
  # % and missing data.
  # %
  # % The model is        y_t   = Z * a_t + eps_t       
  # %                     a_t+1 = TT * a_t + u_t       
  # %
  # %______________________________________________________________________
  # % INPUT  
  # %        Y         Data                                 (nobs x n)  
  # % OUTPUT 
  # %        S.Am       Predicted state vector  A_t|t-1      (nobs x m)  
  # %        S.AmU      Filtered  state vector  A_t|t        (nobs+1 x m)  
  # %        S.Pm       Predicted covariance of A_t|t-1      (nobs x m x m)  
  # %        S.PmU      Filtered  covariance of A_t|t        (nobs+1 x m x m)  
  # %        S.loglik   Value of likelihood function
  # 
  # % Output structure & dimensions
  
  n <- dim(Z)[1]
  m <- dim(Z)[2]
  nobs  <- size(Y,2)
  
  S<-list()
  S$Am <- array(NA,c(m,nobs))
  S$Pm <- array(NA,c(m,m,nobs))
  S$AmU <- array(NA,c(m,nobs+1))
  S$PmU <- array(NA,c(m,m,nobs+1))
  S$loglik <- 0
  
  
  # ______________________________________________________________________
  Au <- A_0;  # A_0|0;
  Pu <- P_0;  # P_0|0
  
  S$AmU[,1]    = Au;
  S$PmU[,,1]  = Pu;
  
  
  
  for(t in 1:nobs){
    #       t
    # A = A_t|t-1   & P = P_t|t-1
    
    A <- TT%*%Au;
    P <- TT%*%Pu%*%t(TT) + Q;
    P <-  0.5*(P+t(P))
    
    # handling the missing data
    res_MissData <- MissData(Y[,t],Z,R)
    
    y_t <- res_MissData$y
    Z_t <- res_MissData$C
    R_t <- res_MissData$R
    L_t <- res_MissData$L
    
    # if(is.null(y_t)){
    if(sum(is.na(y_t))==length(y_t)){
      Au <- A
      Pu <- P
    } else {
      
      
      if(!is.matrix(Z_t)){
        Z_t<-t(as.matrix(Z_t))
      }
      
      
      PZ  <- P%*%t(Z_t)
      # iF  <- ginv(Z_t%*%PZ + R_t)
      iF  <- solve(Z_t%*%PZ + R_t)
      PZF <- PZ%*%iF
      
      V <- y_t - Z_t%*%A
      Au <- A  + PZF%*%V
      Pu <- P  - PZF%*%t(PZ)
      Pu <-  0.5*(Pu+t(Pu))
      
    }
    
    S$Am[,t] <- A
    S$Pm[,,t] <- P
    
    # Au = A_t|t   & Pu = P_t|t
    
    S$AmU[,t+1] <- Au
    S$PmU[,,t+1] <- Pu
    S$loglik <- S$loglik + 0.5*(log(det(iF))  - t(V)%*%iF%*%V)
  } # t
  
  if(sum(is.na(y_t))==length(y_t)){
    S$KZ <- zeros(m,m)
  }else{
    S$KZ <- PZF%*%Z_t
  }
  
  return(S)
  
}

# Main function -----------------------------------------------------------

# library(R.matlab)
# library(matlab)
# library(zoo)
# library(Matrix)
# R_new <- readMat("C:/Users/guilherme.branco/Desktop/EM-transcription/arquivos pra fç EMstep/R_new.mat")
# R_new <- R_new$R.new[,,1]
# names(R_new)[1]<-'X_sm'
# names(R_new)[9]<-'Z_0'
# names(R_new)[10]<-'V_0'
# X_old <- data.frame(read.csv("C:/Users/guilherme.branco/Desktop/EM-transcription/arquivos pra fç EMstep/X_old.csv", header = F))
# X_new <- data.frame(read.csv("C:/Users/guilherme.branco/Desktop/EM-transcription/arquivos pra fç EMstep/X_new.csv", header = F))
# # 
# Q = R_new
# t_fcst = 187   # Qual a data do nowcast?
# # v_news = NULL
# v_news<-"V25" # Qual a série de interesse? Pode ser um vetor?


News_DFM_ML <- function(X_old = NULL, X_new = NULL, Q = NULL, t_fcst = NULL, v_news = NULL){
  
  nomes<-names(X_new)
  r <- ncol(Q$C)
  TT <- nrow(X_new)
  N <- ncol(X_new)
  gList <- unique(unlist(Q$Groups))
  groupnews <- zeros(1,length(gList))
  singlenews <- data.frame(zeros(1,N))
  names(singlenews)<-names(X_new)
  gain <- NULL
  gainSer <- NULL
  actual <- NULL
  fore <- NULL
  filt <- NULL
  
  ### Na nova vintage saiu a informação da variável que estou fazendo forecast na data de interesse?
  if(sum(as.numeric(!is.na(X_new[t_fcst,v_news])))>0){    # Caso permita um vetor como série de interesse
    Res_old <- para_const(X_old, Q, 0)
    temp <- X_new[t_fcst,v_news] - Res_old$X_sm[t_fcst,v_news]
    singlenews[v_news] <- temp
    # groupnews[, gList %in% Q$Groups[v_news]] <- temp
    y_old <- Res_old$X_sm[t_fcst,v_news]
    y_new <- X_new[t_fcst,v_news]
  }else{
    
    Mx = Q$Mx
    Wx = Q$Wx
    
    miss_old <- is.na(X_old)
    miss_new <- is.na(X_new)
    temp <- miss_old - miss_new
    v_miss <- as.vector(which(colSums(temp == 1)==1))
    # x <- t_fcst*v_miss
    t_miss <- as.vector(unlist(apply(temp, MARGIN = 2, FUN = function(x) find(x == 1))))
    lag <- t_fcst - t_miss
    k <- max(c(abs(lag), max(lag)-min(lag)))
    
    C <- Q$C
    R <- t(Q$R)
    
    n_news <- length(lag)
    
    Res_old <- para_const(X_old, Q, k)
    Res_new <- para_const(X_new, Q, 0)
    
    y_old <- Res_old$X_sm[t_fcst,v_news]
    y_new <- Res_new$X_sm[t_fcst,v_news]
    
    # if(!is.empty(t_miss)){
    if(!is.null(t_miss)){
      P <- Res_old$P[,,2:dim(Res_old$P)[3]]
      P1 <- NULL
      for(i in 1:length(lag)){
        h <- abs(t_fcst-t_miss[i])
        m <- max(t_miss[i], t_fcst)
        if(t_miss[i] > t_fcst){
          Pp <- t(P[1:r,(h*r+1):(h*r+r),m])
        }else{
          Pp <- P[1:r,(h*r+1):(h*r+r),m]
        }
        P1 <- cbind(P1, Pp %*% C[v_miss[i],1:r])
      }
      
      innov <- NULL
      for(i in 1:length(t_miss)){
        X_new_norm <- (X_new[t_miss[i],v_miss[i]] - Mx[,v_miss[i]])/Wx[,v_miss[i]]
        X_sm_norm <- (Res_old$X_sm[t_miss[i],v_miss[i]] - Mx[,v_miss[i]])/Wx[,v_miss[i]]
        innov[i] <- X_new_norm-X_sm_norm
      }
      
      ins <- ncol(innov)
      P2 <- NULL
      p2 <- NULL
      WW <- matrix(NA,dim(R)[1],dim(R)[2])
      for(i in 1:length(lag)){
        for(j in 1:length(lag)){
          h <- abs(lag[i]-lag[j])
          m <- max(t_miss[i],t_miss[j])
          if(t_miss[j] > t_miss[i]){
            Pp <- t(P[1:r,(h*r+1):((h+1)*r),m])
          }else{
            Pp <- P[1:r,(h*r+1):((h+1)*r),m]
          }
          if(v_miss[i] == v_miss[j] & t_miss[i] != t_miss[j]){
            WW[v_miss[i],v_miss[j]] <- 0
          }else{
            WW[v_miss[i],v_miss[j]] <- R[v_miss[i],v_miss[j]]
          }
          p2 <- cbind(p2,C[v_miss[i],1:r] %*% Pp %*% C[v_miss[j],1:r] + WW[v_miss[i],v_miss[j]])
        }
        # colnames(p2) <- colnames(P2)
        P2 <- rbind(P2,p2)
        p2 <- NULL
      }
      
      totnews <- Wx[nomes==v_news] %*% C[nomes==v_news,1:r] %*% P1 %*% solve(P2) %*% innov
      temp <- Wx[nomes==v_news] %*% C[nomes==v_news,1:r] %*% P1 %*% solve(P2) * innov
      gain <- Wx[nomes==v_news] %*% C[nomes==v_news,1:r] %*% P1 %*% solve(P2)
      
      # Ainda preciso introduzir esse parâmetro na lista Q
      # gainSer <- Q$Series[v_miss]
      
      singlenews <- zeros(max(t_miss)-min(t_miss)+1,N)
      actual <- zeros(N,1)
      fore <- zeros(N,1)
      filt <- zeros(N,1)
      
      for(i in 1:length(innov)){
        singlenews[t_miss[i]-min(t_miss)+1,v_miss[i]] <- temp[i]
        actual[v_miss[i],1] <- X_new[t_miss[i],v_miss[i]];
        fore[v_miss[i],1] <- Res_old$X_sm[t_miss[i],v_miss[i]]
        filt[v_miss[i],] <- gain[i]/Wx[v_miss[i]]
      }
      
      singlenews <- sum(singlenews)

      # Falta implementar por grupo de notícias      
      # for(i in 1:length(gList)){
      #   groupnews[i] <- gain[, Q$Groups[v_miss] %in% gList[i]] %*% t(innov[Q$Groups[v_miss] %in% gList[i]])
      # }
      # 
      v_miss0 = unique(v_miss)
      idx <- which(v_miss0 %in% v_miss)
      j <- idx
      v_miss <- v_miss0
      gain <- gain[idx]
      # gainSer <- gainSer[idx]
      
    }
  }
  
  # output
  return(list(OldFcst = y_old,
              NewFcst = y_new,
              GroupNews = groupnews,
              SerNews = singlenews,
              gainT = gain,
              serGainT = gainSer,
              Actual = actual,
              Fcst = fore,
              Filt = filt))
  
}


# t_fcst = 188   # Qual a data do nowcast?
# v_news<-c("V25") # Qual a série de interesse? Pode ser um vetor?
# 
# news<-News_DFM_ML(X_old,X_new,Q,t_fcst,v_news)
# 
# # O que é cada output?
# 
# news$OldFcst      # forecast feita com a vintage antiga
# news$NewFcst      # gorecast feita com a vintage nova
# news$GroupNews    # ??
# news$SerNews      # Peso de cada série na resivão final
# news$gainT
# news$serGainT
# news$Actual
# news$Fcst
# news$Filt
