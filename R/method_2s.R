#' @importFrom stats cov end fitted lm median na.omit predict quantile sd start ts tsp window as.ts frequency
#' @import zoo matlab corpcor

bridge <- function(y,x){
  
  # y: ts (trimestral)
  # x: fatores (mensais - output da função FactorExtraction)
  
  # tranformar fatores mensais em trimestrais, selecionando o último fator
  # fatoresTS <- stats::ts(x[,-1], end = as.numeric(c(substr(x[nrow(x),1],1,4),
  #                                                   substr(x[nrow(x),1],6,7))), frequency = 12)
  fatoresTS <- x
  fatoresTRI <- month2qtr(fatoresTS)
  
  
  # estimação do modelo de regressão
  dados <- cbind(y, fatoresTRI)
  colnames(dados) <- c("Y", paste0("Factor",1:ncol(data.frame(fatoresTRI))))
  reg <- stats::lm(Y ~ ., data = na.omit(data.frame(dados)))
  fit <- stats::ts(fitted(reg), end = end(na.omit(dados)), frequency = 4)
  
  Qmax <- max(which(!is.na(dados[,1])))
  edge<-zoo::as.Date(dados)[Qmax]
  
  # previsão
  # newbase <- data.frame(dados[-(1:(Qmax-1)),-1])
  newbase <- data.frame(dados[-(1:Qmax),-1])
  colnames(newbase) <- paste0("Factor",1:ncol(data.frame(fatoresTRI)))
  
  ## função auxiliar
  # tail.ts <- function(data,n) {
  #   data <- as.ts(data)
  #   window(data,start=tsp(data)[2]-(n-1)/frequency(data))
  # }
  
  # ano<-as.numeric(substr(edge,1,4))
  # tri<-as.numeric(substr(quarters(edge),2,2))
  ano<-as.numeric(substr(edge+months(3),1,4))
  tri<-as.numeric(substr(quarters(edge+months(3)),2,2))
  
  prev <- stats::ts(predict(object = reg, newdata = newbase),
                    start = c(ano,tri),
                    frequency = 4) 
  
  dados_pib<-cbind(y,fit,prev)
  
  colnames(dados_pib) <- c("y", "in","out")
  
  
  # RETORNAR PREVISÃO DENTRO E FORA DA AMOSTRA
  list(main = dados_pib,reg = reg)
}

FactorExtraction <- function(x = NULL,q = NULL,r = NULL,p = NULL, 
                             A = NULL,C = NULL,Q = NULL,R = NULL,
                             initx = NULL, initV = NULL,
                             ss = NULL, MM = NULL, n.prevs = NULL){
  
  # The model
  # x_t = C F_t + \xi_t
  # F_t = AF_{t-1} + B u_t
  # R = E(\xi_t \xi_t')
  # Q = BB'
  # u_t ~ WN(0,I_q)
  # initx = F_0
  # initV = E(F_0 F_0')
  # ss: std(x) 
  # MM: mean(x) 
  
  # q: dynamic rank
  # r: static rank (r>=q)
  # p: ar order of the state vector (default p=1)
  
  # Ft: estimated factors
  # VF: estimation variance for the common factors
  
  
  # x = series1
  # q = 2; r = 2; p = 1;
  # A = NULL; C = NULL; Q = NULL;R = NULL;
  # initx = NULL; initV = NULL; ss = NULL; MM = NULL
  
  #finalx <- x
  #x <- finalx[,1:(ncol(finalx)-1)]
  
  # coluna1 <- x[,1]
  # x <- x[,-1]
  
  datas<-zoo::as.Date(x)
  # x<-data.frame(x)
  
  # Base dimension
  TT <- nrow(x)
  N <- ncol(x)
  
  # Number of missings - Here the missings are only new information
  n.missings <- colSums(is.na(x))
  m <- max(n.missings)
  
  # Number of arguments inputed
  n.arg <- sum(c(!is.null(x),!is.null(q),!is.null(r),!is.null(p),
                 !is.null(A),!is.null(C),!is.null(Q),!is.null(R),
                 !is.null(initx), !is.null(initV), 
                 !is.null(ss), !is.null(MM)))
  
  if(n.arg < 5){ # Estimate parameters if they are not inputed
    
    z <- x[1:(TT - m),]      # ONLY complete information for PCA
    s <- apply(z, MARGIN = 2, FUN = sd)
    M <- apply(z, MARGIN = 2, FUN = mean)
    
    for(i in 1:N){
      x[,i] <- (x[,i] - M[i])/s[i]
    }
    z <- x[1:(TT - m),]
    
    parametros <- pcatodfm(z,q,r,p)
    
    A <- parametros$A
    C <- parametros$C
    Q <- parametros$Q
    R <- parametros$R
    initx <- parametros$initx
    initV <- parametros$initV
    a <- parametros$eigen
    #print(A)
  }else{
    
    # If parameters are inputed only need to standardize
    
    for(i in 1:N){
      x[,i] <- (x[,i] - MM[i])/ss[i]
    }
    z <- x[1:(TT - m),]
    
  }
  
  # Parameters for state space model
  
  AA <- array(A, dim = c(nrow(A), ncol(A), TT))
  QQ <- array(Q, dim = c(nrow(Q), ncol(Q), TT))
  CC <- array(C, dim = c(nrow(C), ncol(C), TT))
  miss <- is.na(x)
  RR <- array(NA, dim = c(nrow(R), ncol(R), TT))
  
  for(i in 1:TT){
    Rtemp <- diag(R)
    Rtemp[miss[i,]] <- 1e32
    RR[,,i] <- diag(Rtemp)
  }
  
  xx <- x
  xx[is.na(x)] <- 0 # Arbitrary value for missing
  
  # KALMAN SMOOTHER DIAG
  
  resul <- kalman_smoother_diag(t(xx), AA, CC, QQ, RR, initx, initV, list('model',1:TT))
  
  xsmooth <- resul$xsmooth
  Vsmooth <- resul$Vsmooth
  VVsmooth <- resul$VVsmooth
  loglik <- resul$loglik
  
  VF <- Vsmooth
  ind <- size(VF,3)
  fatores <-  t(xsmooth)
  
  nomes_colunas <- c("data", paste0("Factor",1:ncol(fatores)))
  fator_final <- data.frame(datas, fatores)
  colnames(fator_final) <- nomes_colunas
  
  x<-fator_final
  fatoresTS <- stats::ts(x[,-1], end=as.numeric(c(substr(datas[length(datas)],1,4),
                                                  substr(datas[length(datas)],6,7)))
                         ,frequency = 12)
  
  if(p > 1){
    fatoresTS <- fatoresTS[,1:r]
  }
  list(dynamic_factors = fatoresTS,A = A,C = C,Q = Q,R =  R,initx =  initx,initV =  initV,eigen = a)
}


kalman_filter_diag <- function(y, A, C, Q, R, init_x, init_V, varagin){
  
  # y = t(xx);
  # A = AA;
  # C = CC;
  # Q = QQ;
  # R = RR;
  # init_x = initx;
  # init_V = initV;
  # varagin = list("model", model, "u", u, "B", B)
  
  os <- nrow(y)
  TT <- ncol(y)
  ss <- dim(A)[1] # size of state space
  
  # parâmetros default
  model <- ones(1,TT)
  u <- NULL
  B <- NULL
  ndx <- NULL
  
  args <- varagin
  nargs <- length(args)
  
  for(i in seq(1, nargs, by = nargs-1)){
    if(args[i] == "model"){ model <- args[[i+1]]
    }else if(args[i] == "u"){ u <- args[[i+1]]
    }else if(args[i] == "B"){ B <- args[[i+1]]
    }else if(args[i] == "ndx"){ ndx <- args[[i+1]]
    }
  }
  
  x <- zeros(ss, TT)
  V <- zeros(ss, ss, TT)
  VV <- zeros(ss, ss, TT)
  
  loglik <- 0
  LL <- NULL
  for(t in 1:TT){
    m <- model[t]
    if(t == 1){
      prevx <- init_x
      prevV <- init_V
      initial <- 1
    }else{
      prevx <- t(matrix(x[,t-1]))
      prevV <- V[,,t-1]
      initial <- 0
    }
    
    if(is.null(u)){
      
      #print(paste0("entrei aqui no primeiro", t))
      resul <- kalman_update_diag(A[,,m], C[,,m], Q[,,m], R[,,m], y[,t], prevx, prevV, list('initial', initial))
      x[,t] <- resul$xnew
      V[,,t] <- resul$Vnew
      LL[t] <- resul$loglik
      VV[,,t] <- resul$VVnew
      
    }else if(is.null(ndx)){
      
      #print(paste0("entrei aqui em null ndx", t))
      resul <- kalman_update_diag(A[,,m], C[,,m], Q[,,m], R[,,m], y[,t], prevx, prevV, list('initial', initial, 'u', u[,t], 'B', B[,,m]))
      x[,t] <- resul$xnew
      V[,,t] <- resul$Vnew
      LL[t] <- resul$loglik
      VV[,,t] <- resul$VVnew
      
    }else{
      
      #print(paste0("entrei aqui em else", t))
      i <- ndx[t];
      # copy over all elements; only some will get updated
      x[,t] <- prevx;
      prevP <- solve(prevV);
      prevPsmall <- prevP[i,i]
      prevVsmall <- solve(prevPsmall)
      
      resul <- kalman_update_diag(A[i,i,m], C[,i,m], Q[i,i,m], R[,,m], y[,t], prevx[i], prevVsmall, list('initial', initial, 'u', u[,t], 'B', B[i,,m]))
      x[i,t] <- resul[[1]]
      smallV <- resul[[5]]
      LL[t] <- resul[[3]]
      VV[i,i,t] <- resul[[4]]
      
      smallP <- solve(smallV)
      prevP[i,i] <- smallP
      V[,,t] <- solve(prevP)
    }
    
    loglik <- loglik + LL[t]
  }
  
  list(x = x, V = V, VV = VV, loglik = loglik)
  
  
}

kalman_smoother_diag <- function(y, A, C, Q, R, init_x, init_V, varagin){
  
  # y = t(xx);
  # A = AA;
  # C = CC;
  # Q = QQ;
  # R = RR;
  # init_x = initx;
  # init_V = initV;
  # varagin = list("model", 1:TT)
  
  os <- nrow(y)
  TT <- ncol(y)
  ss <- dim(A)[1]
  
  # parâmetros default
  model <- ones(1,TT)
  u <- NULL
  B <- NULL
  
  args = varagin
  nargs = length(args)
  
  for(i in seq(1, nargs, by = nargs-1)){
    if(args[i] == "model"){ model <- args[[i+1]]
    }else if(args[i] == "u"){ u <- args[[i+1]]
    }else if(args[i] == "B"){ B <- args[[i+1]]
    }
  }
  
  xsmooth <- zeros(ss, TT)
  Vsmooth <- zeros(ss, ss, TT)
  VVsmooth <- zeros(ss, ss, TT)
  
  # Forward pass
  resul <- kalman_filter_diag(y, A, C, Q, R, init_x, init_V, list("model", model, "u", u, "B", B))
  xfilt <- resul$x
  Vfilt <- resul$V
  VVfilt <- resul$VV
  loglik <- resul$loglik
  
  # Backward pass
  xsmooth[,TT] <- xfilt[,TT]
  Vsmooth[,,TT] <- Vfilt[,,TT]
  
  for(t in (TT - 1):1){
    m <- model[t+1]
    if(is.null(B)){
      
      resul <- smooth_update(xsmooth[,t+1], Vsmooth[,,t+1], xfilt[,t], Vfilt[,,t], 
                             Vfilt[,,t+1], VVfilt[,,t+1], A[,,m], Q[,,m], NULL, NULL)
      
      xsmooth[,t] <- resul$xsmooth
      Vsmooth[,,t] <- resul$Vsmooth
      VVsmooth[,,t+1] <- resul$VVsmooth
      
    }else{
      
      resul <- smooth_update(xsmooth[,t+1], Vsmooth[,,t+1], xfilt[,t], Vfilt[,,t], 
                             Vfilt[,,t+1], VVfilt[,,t+1], A[,,m], Q[,,m],  B[,,m], u[,t+1])
      
      xsmooth[,t] <- resul$xsmooth
      Vsmooth[,,t] <- resul$Vsmooth
      VVsmooth[,,t+1] <- resul$VVsmooth
      
    }
  }
  
  VVsmooth[,,1] <- zeros(ss,ss)
  
  list(xsmooth = xsmooth, Vsmooth = Vsmooth, VVsmooth = VVsmooth, loglik = loglik)
}

kalman_update_diag <- function(A, C, Q, R, y, x, V, varagin){
  
  # A = A[,,m]
  # C = C[,,m]
  # Q = Q[,,m]
  # R = R[,,m]
  # y = y[,t]
  # x = prevx
  # V = prevV
  # varagin <- list('initial', initial)
  
  # parâmetros default
  u <- NULL
  B <- NULL
  initial <- 0
  
  args <- varagin
  nargs <- length(args)
  
  for(i in seq(1, nargs, by = nargs-1)){
    if(args[i] == "u"){ u <- args[[i+1]]
    }else if(args[i] == "B"){ B <- args[[i+1]]
    }else if(args[i] == "initial"){ initial <- args[[i+1]]
    }
  }
  
  #  xpred(:) = E[X_t+1 | y(:, 1:t)]
  #  Vpred(:,:) = Cov[X_t+1 | y(:, 1:t)]
  
  if(initial != 0){
    if(is.null(u)){ 
      xpred <- x
    }else{ 
      xpred <- x + B %*% u
    }
    Vpred <- V
  }else{
    if(is.null(u)){
      xpred <- t(A %*% t(x))
    }else{
      xpred <- t(A %*% t(x) + B %*% u)
    }
    Vpred <- A %*% V %*% t(A) + Q
  }
  
  e <- y - C %*% t(xpred) # error (innovation)
  n <- length(e)
  ss <- nrow(A)
  if(is.null(ss)){ ss <- length(A) }
  
  d <- size(e,1)
  
  S <- C %*% Vpred %*% t(C) + R
  GG <- t(C) %*% diag(1/diag(R)) %*% C
  Sinv <- diag(1/diag(R)) - diag(1/diag(R)) %*% C %*% corpcor::pseudoinverse(eye(ss) + Vpred %*% GG) %*% Vpred %*% t(C) %*% diag(1/diag(R)) # works only with R diagonal
  
  # Sinv = inv(S)
  
  ################################################
  
  detS <- prod(diag(R)) %*% det(eye(ss) + Vpred %*% GG)
  denom <- (2*pi)^(d/2)*sqrt(abs(detS))
  mahal <- rowSums(t(e) %*% Sinv %*% e)
  loglik <- -0.5*mahal - log(denom)
  
  ################################################
  
  K <- Vpred %*% t(C) %*% Sinv # Kalman gain matrix
  
  # If there is no observation vector, set K = zeros(ss).
  xnew <- t(xpred) + (K %*% e)              #csi_est(t\t) formula 13.6. 5    
  Vnew <- (eye(ss) - K %*% C) %*% Vpred    #P(t\t) formula 13.2.16 hamilton
  VVnew <- (eye(ss) - K %*% C) %*% A %*% V
  
  list(xnew = xnew, Vnew = Vnew, VVnew = VVnew, loglik = loglik)          
  
}

outliers_correction <- function(x, k_ma = 3){
  # x é um série temporal
  
  # encontrar missings
  missing <- is.na(x)
  
  # Função criada no R para superar os problemas da função do matlab
  # outlier são as obs que ultrapassam 4* distância interquartilica
  outlier <- abs(x - median(x, na.rm = T)) > 4*abs(quantile(x, probs = 0.25, na.rm = T) -  quantile(x, probs = 0.75, na.rm = T)) & !missing
  
  ### Problema 2: Usa NA para calcular o tamanho do vetor
  # TT <- length(x)   # problem 1
  # TT <- sum(!is.na(x))   # Solução: utilizar apenas as observações completas
  
  ### Problema 3: utilizar a função round do matlab gera viés
  # round2 = function(x, n) {
  #   posneg = sign(x)
  #   z = abs(x)*10^n
  #   z = z + 0.5
  #   z = trunc(z)
  #   z = z/10^n
  #   z*posneg
  # }
  # outlier <- abs(x - sort(x)[round2(TT*1/2,0)]) > 4*abs(sort(x)[round2(TT*1/4,0)] -  sort(x)[round2(TT*3/4,0)]) & !missing
  # possível solução: utilizar a função round do R que controla pelo viés.
  
  Z <- x
  
  # substituir outliers e missings pela mediana
  Z[outlier] <- median(x, na.rm = T)
  Z[missing] <- median(x, na.rm = T)
  # Z[outlier] <- sort(x)[round2(TT*1/2,0)]
  # Z[missing] <- sort(x)[round2(TT*1/2,0)]
  
  # Média móvel de ordem K
  xpad <- c(Z[1]*ones(k_ma,1), Z, Z[length(Z)]*ones(k_ma,1))
  x_ma <- xpad*NA
  for(j in (k_ma + 1):(length(xpad) - k_ma)){
    x_ma[j - k_ma] = mean(xpad[(j - k_ma):(j + k_ma)])
  }
  x_ma <- x_ma[1:length(x)]
  
  Z[outlier] <- x_ma[outlier]
  Z[missing] <- x_ma[missing]
  
  Z
}

pcatodfm <- function(x, q, r, p){
  
  # x é a base de dados em formato ts
  # q é o número de choques nos fatores
  # r é a quantidade de fatores
  # p é o grau do polinômio autorregressivo
  # x = x_padronizado
  # r = 2
  # q = 2
  # p = 1
  
  # finalx2 = z
  # x = z
  x <- as.matrix(x)
  Mx <- colMeans(x)
  Wx <- apply(x, MARGIN = 2, FUN = sd)
  
  for(i in 1:ncol(x)){
    x[,i] <- (x[,i] - Mx[i])/Wx[i]
  }
  
  # tamanho da base
  TT <- nrow(x)
  N <- ncol(x)
  
  # restrição
  if(r < q){ stop("q must be less than or equal to r") }
  
  # nlag 
  nlag <- p - 1
  
  # A temporária - matriz de zeros
  A_temp <- t(zeros(r,r*p))
  # matriz identidade
  I <- diag(r*p)
  if(p == 1){ 
    A <- A_temp 
  }else{
    A <- rbind(t(A_temp), I[1:(nrow(I) - r), 1:ncol(I)]) 
  }
  
  Q <- zeros(r*p,r*p)
  Q[1:r,1:r] <- diag(r)
  
  # autovalores e autovetores
  a <- eigen(cov(x))
  a1 <- a    # save eigen
  d <- a$values[1:r]
  v <- a$vectors[,1:r]
  
  # estimativa dos fatores comuns
  EF <- x %*% v
  # estimativa da matriz de covariância do termo de erro na equação dos fatores comuns
  R <- diag(diag(cov(x - x %*% v %*% t(v))))
  
  if(p > 0){
    # estimar o modelo autoregressivo para os fatores VAR: F(t) =  A1*F(t-1) + ... + Ap*F(t-p) + e(t)
    z <- EF
    Z <- NULL
    for(kk in 1:p){
      Z <- cbind(Z, z[(p - kk + 1):(nrow(z) - kk),])
    }
    z <- z[(p + 1):nrow(z),]
    
    # Estimador OLS para a matriz de transição do VAR
    A_temp <- solve(t(Z) %*% Z) %*% t(Z) %*% z
    A[1:r,1:(r*p)] <- t(A_temp)
    
    # Q
    e <- z - Z%*%A_temp # VAR residuals
    H <- cov(e) # VAR covariance matrix
    
    if(r == q){ # The covariance matrix of the VAR residuals is of full rank
      
      Q[1:r,1:r] <- H
      
    }else{ #The covariance matrix of the VAR residuals has reduced rank
      
      a <- eigen(H)
      P <- a$vectors[,1:q]
      
      if(is.matrix(P)){
        M <- diag(a$values[1:q])
        P <- P %*% diag(sign(P[1,]))
        Q[1:r,1:r] <- P %*% M %*% t(P)
      }else{
        M <- a$values[1:q]
        P <- P * sign(P[1])
        Q[1:r,1:r] <- P * M * t(P)
      }
      
      u_orth <- e %*% P %*% (M^(-.5)) # extracting the common shocks
      e_pc <- e %*% P %*% t(P)
      #Q[1:r,1:r] <- P %*% M %*% t(P)
      
    }
  }
  
  # Condições iniciais pro filtro de kalman
  
  if(p > 0){
    
    z <- EF
    Z <- NULL
    
    for(kk in 0:nlag){
      Z <- cbind(Z, z[(nlag - kk + 1):(nrow(z) - kk),]) # stacked regressors (lagged SPC)
    }
    
    initx <- t(Z[1,]) 
    initV <- matlab::reshape(corpcor::pseudoinverse(eye(size(kronecker(A,A),1))- kronecker(A,A)) %*% matrix(Q, ncol = 1), r*p, r*p)
    
  }else{
    
    initx <- NA
    initV <- NA
    
  }
  
  if(nlag != 0){
    C <- cbind(v, zeros(N,r*(nlag)))
  }else{
    C <- as.matrix(v)
  }
  
  list(A = A, C = C, Q = Q, R = R, initx = initx, initV = initV, eigen = a1)
}

smooth_update <- function(xsmooth_future, Vsmooth_future, xfilt, Vfilt,  Vfilt_future, VVfilt_future, A, Q, B, u){
  
  # xsmooth_future = xsmooth[,t+1]
  # Vsmooth_future = Vsmooth[,,t+1]
  # xfilt = xfilt[,t]
  # Vfilt_future = Vfilt[,,t+1]
  # Vfilt = Vfilt[,,t]
  # 
  # VVfilt_future = VVfilt[,,t+1]
  # A = A[,,m]
  # Q = Q[,,m]
  # B = NULL
  # u = NULL
  
  
  if(is.null(B)){
    xpred <- A %*% xfilt
  }else{
    xpred <- A %*% xfilt + B %*% u
  }
  
  Vpred <- A %*% Vfilt %*% t(A) + Q # Vpred = Cov[X(t+1) | t]
  J <- Vfilt %*% t(A) %*% corpcor::pseudoinverse(Vpred) # smoother gain matrix
  xsmooth <- xfilt + J %*% (xsmooth_future - xpred)
  Vsmooth <- Vfilt + J %*% (Vsmooth_future - Vpred) %*% t(J)
  
  if(is.matrix(Vfilt_future)){
    VVsmooth_future <- VVfilt_future + (Vsmooth_future - Vfilt_future) %*% corpcor::pseudoinverse(Vfilt_future) %*% VVfilt_future
  }else{
    VVsmooth_future <- VVfilt_future + (Vsmooth_future - Vfilt_future) %*% 1/Vfilt_future %*% VVfilt_future
    
  }
  list(xsmooth = xsmooth, Vsmooth = Vsmooth, VVsmooth_future = VVsmooth_future)
}
