#' @importFrom magic adiag


remNaNs_spline <-function(X,options){
  
  TT <- dim(X)[1]
  N <- dim(X)[2]
  k <- options$k
  indNaN <- is.na(X)
  
  if(options$method == 1){ # replace all the missing values (método Giannone et al. 2008)
    for (i in 1:N){  
      x = X[,i]
      x[indNaN[,i]] = median(x,na.rm = T);
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]]
      X[,i]=x;
    }
    
  }else if(options$method == 2){ # replace missing values after removing leading and closing zeros
    
    rem1 <- (rowSums(indNaN)>N*0.8)
    nanLead <- which(rem1)
    # nanEnd <- which(rem1[length(rem1):1])
    # nanLE <- c(nanEnd,nanLead)
    nanLE<-nanLead
    X<-X[-nanLE,]
    indNaN=is.na(X)
    for (i in 1:N){  
      x = X[,i]
      isnanx = is.na(x)
      t1 = min(which(!isnanx))
      t2 = max(which(!isnanx))
      
      x1<-stats::spline(x[t1:t2],xout = 1:(t2-t1+1))
      xx<-x1$y
      x[t1:t2]<-x1$y
      isnanx<-is.na(x)
      x[isnanx] <- median(x,na.rm = T)
      
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]]
      X[,i]=x;
    }
    
  }else if(options$method == 3){ # only remove rows with leading and closing zeros
    
    rem1 <- (rowSums(indNaN)==N)
    nanLead <- which(rem1)
    # nanEnd <- which(rem1[length(rem1):1])
    # nanLE <- c(nanEnd,nanLead)
    nanLE<-nanLead
    if(length(nanLE) != 0){
      X <- X[-nanLE,]
    }
    indNaN=is.na(X)
    
  }else if(options$method == 4){ # remove rows with leading and closing zeros & replace missing values
    
    rem1 <- (rowSums(indNaN)==N)
    nanLead <- which(rem1)
    # nanEnd <- which(rem1[length(rem1):1])
    # nanLE <- c(nanEnd,nanLead)
    nanLE<-nanLead
    X<-X[-nanLE,]
    indNaN=is.na(X)
    for (i in 1:N){  
      x = X[,i]
      isnanx = is.na(x)
      t1 = min(which(!isnanx))
      t2 = max(which(!isnanx))
      
      x1<-stats::spline(x[t1:t2],xout = 1:(t2-t1+1))
      xx<-x1$y
      x[t1:t2]<-x1$y
      isnanx<-is.na(x)
      x[isnanx] <- median(x,na.rm = T)
      
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]]
      X[,i]=x;
    }
    
  }else if(options$method == 5){
    indNaN=is.na(X)
    
    for (i in 1:N){  
      x = X[,i]
      isnanx = is.na(x)
      t1 = min(which(!isnanx))
      t2 = max(which(!isnanx))
      
      x1<-stats::spline(x[t1:t2],xout = 1:(t2-t1+1))
      xx<-x1$y
      x[t1:t2]<-x1$y
      isnanx<-is.na(x)
      x[isnanx] <- median(x,na.rm = T)
      
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]]
      X[,i]=x;
    }
  }
  
  return(list(X = X,indNaN=indNaN))
  
}

InitCond<-function(xNaN,r,p,blocks,optNaN,R_mat,q,nQ,i_idio){

  x<-xNaN
  Rcon<-R_mat
  
  # library(magic)
  
  pC = size(Rcon,2)
  ppC = max(p,pC)
  n_b = size(blocks,2)
  
  OPTS<-list()
  OPTS$disp=0
  
  res_remNaNs_spline <- remNaNs_spline(x,optNaN)
  xBal<- res_remNaNs_spline$X
  indNaN <- res_remNaNs_spline$indNaN
  TT <- dim(xBal)[1]
  N <- dim(xBal)[2]
  NM <- N-nQ
  
  xNaN = xBal
  
  for(i in 1:N){
    xNaN[indNaN[,i],i] <- NA
  }  
  
  
  C = {}
  A = {}
  Q = {}
  initV = {}
  
  res = xBal
  resNaN = xNaN
  indNaN[1:pC-1,] <- T
  
  
  for(i in 1:n_b){    # roda o loop em cada um dos blocos (geral, real e nominal)
    r_i<-r[i]
    
    ########################
    # Observation equation #
    ########################
    
    C_i = zeros(N,r_i*ppC)
    idx_i = find(blocks[,i])
    
    # Atenção aqui funciona pois a(s) última(S) variável(is) da base de dados é(são) trimestral(is)!
    idx_iM = idx_i[idx_i<NM+1];   # índice representando variável mesal
    idx_iQ = idx_i[idx_i>NM];     # índice representando variável trimestral
    
    eig<-eigen(cov(res[,idx_iM]))
    v<-eig$vectors[,1:r_i]
    d<-eig$values[1:r_i]
    
    C_i[idx_iM,1:r_i] = v
    f = as.matrix(res[,idx_iM])%*%as.matrix(v)
    for(kk in 0:(max(p+1,pC)-1)){
      if(kk == 0){
        FF<-f[(pC-kk):(dim(f)[1]-kk),]
      }else{
        FF <- cbind(FF,f[(pC-kk):(dim(f)[1]-kk),])
      }
    }
    
    Rcon_i = kronecker(Rcon,eye(r_i))
    q_i = kronecker(q,zeros(r_i,1));
    ff = FF[,1:(r_i*pC)]
    
    for(j in idx_iQ){     # Coeficiente "loadings" de Variáveis trimestrais
      xx_j = resNaN[pC:dim(resNaN)[1],j]
      if(sum(!is.na(xx_j)) < size(ff,2)+2){
        xx_j = res[pC:dim(res)[1],j]
      }
      ff_j = ff[!is.na(xx_j),]
      xx_j = xx_j[!is.na(xx_j)]
      iff_j = solve(t(ff_j)%*%ff_j)
      Cc = iff_j%*%t(ff_j)%*%xx_j
      Cc = Cc - iff_j%*%t(Rcon_i)%*%solve(Rcon_i%*%iff_j%*%t(Rcon_i))%*%(Rcon_i%*%Cc-q_i);
      C_i[j,1:(pC*r_i)] <- t(Cc)
    }
    
    ff = rbind(zeros(pC-1,pC*r_i),ff)
    res = res - ff%*%t(C_i)
    resNaN = res
    for(i_aux in 1:dim(indNaN)[2]){
      resNaN[indNaN[,i_aux],i_aux] <- NA
    }
    C <- cbind(C,C_i)
    
    #######################    
    # Transition Equation #
    #######################  
    z <- FF[,1:r_i]
    Z <- FF[,(r_i+1):(r_i*(p+1))]
    A_i = t(zeros(r_i*ppC,r_i*ppC))
    A_temp = solve(t(Z)%*%Z)%*%t(Z)%*%z
    A_i[1:r_i,1:(r_i*p)] <- t(A_temp)
    A_i[(r_i+1):dim(A_i)[1],1:(r_i*(ppC-1))] <- eye(r_i*(ppC-1))
    
    ##########################
    
    Q_i = zeros(ppC*r_i,ppC*r_i)
    e = z  - Z%*%A_temp         # VAR residuals
    Q_i[1:r_i,1:r_i] = cov(e);  # VAR covariance matrix
    
    initV_i = matlab::reshape(solve(eye((r_i*ppC)^2)-kronecker(A_i,A_i))%*%c(Q_i),r_i*ppC,r_i*ppC);
    
    if(is.null(A)){
      A<-A_i
    }else{
      A <- magic::adiag(A,A_i)  
    }
    
    if(is.null(Q)){
      Q<-Q_i
    }else{
      Q <- magic::adiag(Q,Q_i)  
    }
    
    if(is.null(initV)){
      initV<-initV_i
    }else{
      initV <- magic::adiag(initV,initV_i)  
    }
    
    # linha 401 do código em matlab
    
  }
  
  
  R = diag(apply(resNaN, 2, stats::var, na.rm = T))
  
  eyeN = eye(N)
  eyeN<-eyeN[,i_idio]
  
  # Initial conditions
  C=cbind(C,eyeN)
  
  ii_idio = find(i_idio)
  n_idio = length(ii_idio)
  B = zeros(n_idio)
  S = zeros(n_idio)
  
  BM<-zeros(n_idio)
  SM<-zeros(n_idio)
  
  for (i in 1:n_idio){
    R[ii_idio[i],ii_idio[i]] <- 1e-04
    
    res_i = resNaN[,ii_idio[i]]
    # number of leading zeros
    # leadZero = max( find( t(1:TT) == cumsum(is.na(res_i)) ) )
    # endZero = max( find( TT:1 == cumsum(is.na(res_i[length(res_i):1])) ) );
    # res_i<-res_i[(leadZero+1):(length(res_i)-endZero)]    
    
    res_i<-res_i[!is.na(res_i)]
    
    BM[i,i] = solve(t(res_i[1:(length(res_i)-1)])%*%res_i[1:(length(res_i)-1)])%*%t(res_i[1:(length(res_i)-1)])%*%res_i[2:length(res_i)] 
    SM[i,i] = stats::var(res_i[2:length(res_i)]-res_i[1:(length(res_i)-1)]*BM[i,i])
    # SM[i,i] = var(res_i[2:length(res_i)]-res_i[1:(length(res_i)-1)]*B[i,i])
    # ATENÇÃO: Aqui os autores usam B[i,i], porém esse valor é 0. Então eu uso BM[i,i]
  }
  
  initViM = diag(1/diag(eye(size(BM,1))-BM^2))%*%SM;
  
  
  C<-cbind(C,rbind(zeros(NM,5*nQ),t(kronecker(eye(nQ),c(1,2,3,2,1)))))
  Rdiag<-diag(R)
  sig_e <- Rdiag[(NM+1):N]/19
  Rdiag[(NM+1):N] <- 1e-04
  R = diag(Rdiag)
  
  rho0<-0.1
  
  BQ <- kronecker(eye(nQ),rbind(cbind(rho0,zeros(1,4)),cbind(eye(4),zeros(4,1))))
  temp = zeros(5)
  temp[1,1] = 1
  if(is.matrix(sig_e)){
    SQ = kronecker(diag((1-rho0^2)*sig_e),temp)
  }else{
    SQ = kronecker((1-rho0^2)*sig_e,temp)
  }
  
  initViQ = matlab::reshape(solve(eye((5*nQ)^2)-kronecker(BQ,BQ))%*%c(SQ),5*nQ,5*nQ)
  
  # BQ = kronecker(eye(nQ),rbind(zeros(1,5),cbind(eye(4),zeros(4,1))))
  # temp = zeros(5)
  # temp[1,1] = 1
  # if(is.matrix(sig_e)){
  #   SQ = kronecker(diag(sig_e),temp)
  # }else{
  #   SQ = kronecker(diag(as.matrix(sig_e)),temp)  
  # }
  # temp = matrix(c(19, 16, 10, 4, 1, 16, 19, 16, 10, 4, 10, 16, 19, 16, 10, 4, 10, 16, 19, 16, 1, 4, 10, 16, 19),
  #               5,5)
  # if(is.matrix(sig_e)){
  #   initViQ = kronecker(diag(sig_e),temp);
  # }else{
  #   initViQ = kronecker(sig_e,temp);
  # }
  
  A1<-magic::adiag(A,BM,BQ)
  Q1<-magic::adiag(Q, SM, SQ)
  
  A<-A1
  Q<-Q1
  
  # Initial conditions
  initZ = zeros(size(A,1),1); ##[randn(1,r*(nlag+1))]';
  initV = magic::adiag(initV, initViM, initViQ)
  
  return(list(A = A, C = C, Q = Q, R = R, initZ = initZ, initV = initV))
  
}

EM_DFM_SS_block_idioQARMA_restrMQ<-function(X,Par){
  
  # library(matlab)
  
  thresh = 1e-4
  r = Par$r
  p = Par$p
  max_iter = Par$max_iter
  i_idio = Par$i_idio
  R_mat = Par$Rconstr
  q = Par$q
  nQ = Par$nQ
  blocks = Par$blocks
  
  ### Prepara??o dos dados
  
  TT <- dim(X)[1]
  N <- dim(X)[2]
  
  ### Standardise X
  Mx = colMeans(X,na.rm=T)
  Wx = sapply(1:N,function(x) sd(X[,x],na.rm = T))
  # xNaN = (X-repmat(Mx,TT,1))/repmat(Wx,TT,1)
  xNaN <- (X-kronecker(t(Mx),rep(1,TT)))/kronecker(t(Wx),rep(1,TT))
  
  ### Initial conditions
  
  # Removing missing values (for initial estimator)
  optNaN<-list()
  optNaN$method = 2; # Remove leading and closing zeros
  optNaN$k = 3;
  
  res_InitCond<-InitCond(xNaN,r,p,blocks,optNaN,R_mat,q,nQ,i_idio);
  
  A<-res_InitCond$A
  C<-res_InitCond$C
  Q<-res_InitCond$Q
  R<-res_InitCond$R
  Z_0<-res_InitCond$initZ
  V_0<-res_InitCond$initV
  
  # some auxiliary variables for the iterations
  previous_loglik = -Inf
  num_iter = 0
  LL = -Inf
  converged = F
  
  # y for the estimation is WITH missing data
  y = t(xNaN)
  
  #--------------------------------------------------------------------------
  #THE EM LOOP
  #--------------------------------------------------------------------------
  
  #The model can be written as
  #y = C*Z + e;
  #Z = A*Z(-1) + v
  #where y is NxT, Z is (pr)xT, etc
  
  #remove the leading and ending nans for the estimation
  optNaN$method = 3
  y_est <- remNaNs_spline(xNaN,optNaN)
  y_est_indNaN<-t(y_est$indNaN)
  y_est<-t(y_est$X)
  
  
  while ((num_iter < max_iter) & !converged){
    
    # message(num_iter)

    res_EMstep = EMstep(y_est, A, C, Q, R, Z_0, V_0, r,p,R_mat,q,nQ,i_idio,blocks)
    # res_EMstep <- list(C_new, R_new, A_new, Q_new, Z_0, V_0, loglik)
    
    C = res_EMstep$C_new;
    R = res_EMstep$R_new;
    A = res_EMstep$A_new;
    Q = res_EMstep$Q_new;
    Z_0<-res_EMstep$Z_0
    V_0<-res_EMstep$V_0
    loglik<-res_EMstep$loglik
    
    # Checking convergence
    if (num_iter>2){
      res_em_converged = em_converged(loglik, previous_loglik, thresh,1)
      # res_em_converged<-list(converged,decrease[num_iter+1])
      
      converged<-res_em_converged$converged
      decreasse<-res_em_converged$decrease
    }
    
    LL <- cbind(LL, loglik)
    previous_loglik <- loglik
    num_iter <-  num_iter + 1
  }
  
  # final run of the Kalman filter
  res_runKF = runKF(y, A, C, Q, R, Z_0, V_0)
  # res_runKF = runKF(y_est, A, C, Q, R, Z_0, V_0)
  Zsmooth<-t(res_runKF$xsmooth)
  x_sm <- Zsmooth[2:dim(Zsmooth)[1],]%*%t(C)
  
  Res<-list()
  
  # Res$X_sm <- repmat(Wx,TT,1)*x_sm+repmat(Mx,TT,1)
  Res$X_sm <- kronecker(t(Wx),rep(1,TT))*x_sm + kronecker(t(Mx),rep(1,TT))
  Res$FF <- Zsmooth[2:dim(Zsmooth)[1],]
  
  #--------------------------------------------------------------------------
  #  Loading the structure with the results
  #--------------------------------------------------------------------------
  Res$C <- C;
  Res$R <- R;
  Res$A <- A;
  Res$Q <- Q;
  Res$Mx <- Mx;
  Res$Wx <- Wx;
  Res$Z_0 <- Z_0;
  Res$V_0 <- V_0;
  Res$r <- r;
  Res$p <- p;
  
  return(Res)
  
  
}

EMstep <- function(y = NULL, A = NULL, C = NULL, Q = NULL, R = NULL, Z_0 = NULL, V_0 = NULL, 
                   r = NULL, p = NULL, R_mat = NULL, q = NULL, nQ = NULL, i_idio = NULL, blocks = NULL){
  
  # y=y_est
  
  # message('EMstep antes dos parâmetros')
  
  n <- size(y,1)
  TT <- size(y,2)
  nM <- n - nQ
  pC <- size(R_mat,2)
  ppC <- max(p,pC)
  n_b <- size(blocks,2)
  
  # Compute the (expected) sufficient statistics for a single Kalman filter sequence.
  
  #Running the Kalman filter with the current estimates of the parameters
  res_runKF = runKF(y, A, C, Q, R, Z_0, V_0);
  
  
  Zsmooth<-res_runKF$xsmooth
  Vsmooth<-res_runKF$Vsmooth
  VVsmooth<-res_runKF$VVsmooth
  loglik<-res_runKF$loglik
  
  A_new <- A
  Q_new <- Q
  V_0_new <- V_0
  
  # message('EMstep antes loop 1:nb')
  
  
  for(i in 1:n_b){
    
    # message(i)
    
    r_i <- r[i]
    rp <- r_i*p
    if(i == 1){
      rp1 <- 0*ppC
    }else{
      rp1 <- sum(r[1:(i-1)])*ppC
    }
    A_i <- A[(rp1+1):(rp1+r_i*ppC), (rp1+1):(rp1+r_i*ppC)]
    Q_i <- Q[(rp1+1):(rp1+r_i*ppC), (rp1+1):(rp1+r_i*ppC)]
    
    if(rp==1){
    EZZ <- t(Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) %*% Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)] +
      sum(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),2:dim(Vsmooth)[3]])  # E(Z'Z)
    EZZ_BB <- t(Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)]) %*% Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)] + 
      sum(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),1:(dim(Vsmooth)[3]-1)]) #E(Z(-1)'Z_(-1))
    EZZ_FB <- t(Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) %*% Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)] + 
      sum(VVsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),]) #E(Z'Z_(-1))
  
    }else{
    EZZ <- (Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) %*% t(Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) +
      apply(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),2:dim(Vsmooth)[3]],c(1,2),sum)  # E(Z'Z)
    EZZ_BB <- (Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)]) %*% t(Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)]) + 
      apply(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),1:(dim(Vsmooth)[3]-1)],c(1,2),sum) #E(Z(-1)'Z_(-1))
    EZZ_FB <- (Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) %*% t(Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)]) + 
      apply(VVsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),],c(1,2),sum) #E(Z'Z_(-1))
    }
    # message('após EZZ-.')
    
    

    A_i[1:r_i,1:rp] <- EZZ_FB[1:r_i,1:rp] %*% solve(EZZ_BB[1:rp,1:rp])
    Q_i[1:r_i,1:r_i] <- (EZZ[1:r_i,1:r_i] - A_i[1:r_i,1:rp] %*% t(matrix(EZZ_FB[1:r_i,1:rp],r_i,rp))) / TT
    
    # message('depois de Q_i')
    
    A_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] <- A_i 
    Q_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] <- Q_i;
    V_0_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] <- Vsmooth[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC),1]
    
    # message('depois de V_0_new')
    
  }
  
  # message('EMstep depois loop 1:nb')
  
  rp1 <- sum(r)*ppC
  niM <- sum(i_idio[1:nM])
  
  # idiosyncratic
  EZZ <- diag(diag(Zsmooth[(rp1+1):nrow(Zsmooth),2:ncol(Zsmooth)] %*% t(Zsmooth[(rp1+1):nrow(Zsmooth),2:ncol(Zsmooth)]))) + diag(diag(apply(Vsmooth[(rp1+1):dim(Vsmooth)[1],(rp1+1):dim(Vsmooth)[2],2:dim(Vsmooth)[3]], MARGIN = 1:2, FUN = sum))) #E(Z'Z)
  EZZ_BB <- diag(diag(Zsmooth[(rp1+1):nrow(Zsmooth),1:(ncol(Zsmooth)-1)] %*% t(Zsmooth[(rp1+1):nrow(Zsmooth),1:(ncol(Zsmooth)-1)]))) + diag(diag(apply(Vsmooth[(rp1+1):dim(Vsmooth)[1],(rp1+1):dim(Vsmooth)[2],1:(dim(Vsmooth)[3]-1)], MARGIN = 1:2, FUN = sum)))  #E(Z(-1)'Z_(-1))
  EZZ_FB <- diag(diag(Zsmooth[(rp1+1):nrow(Zsmooth),2:ncol(Zsmooth)] %*% t(Zsmooth[(rp1+1):nrow(Zsmooth),1:(ncol(Zsmooth)-1)]))) + diag(diag(apply(VVsmooth[(rp1+1):dim(VVsmooth)[1],(rp1+1):dim(VVsmooth)[2],], MARGIN = 1:2, FUN = sum))) #E(Z'Z_(-1)) 
  
  A_i <- EZZ_FB %*% diag(1/diag(EZZ_BB))
  Q_i <- (EZZ - A_i %*% t(EZZ_FB)) / TT
  
  A_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] <- A_i[1:niM,1:niM] 
  Q_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] <- Q_i[1:niM,1:niM]
  V_0_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] <- diag(diag(Vsmooth[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM),1]))
  
  Z_0 <- Zsmooth[,1] #zeros(size(Zsmooth,1),1); #
  
  # nanY <- is.nan(y)
  nanY<-is.na(y)
  y[nanY] <- 0
  
  # LOADINGS
  C_new <- C
  
  # message('EMstep antes Blocks')
  
  # Blocks
  bl <- unique(blocks)
  n_bl <- size(bl,1)
  bl_idxM <- NULL
  bl_idxQ <- NULL
  R_con <- NULL
  q_con <- NULL
  
  # message('EMstep antes segundo loop 1:nb')
  
  for(i in 1:n_b){
    bl_idxQ <- cbind(bl_idxQ, repmat(bl[,i],1,r[i]*ppC))
    bl_idxM <- cbind(bl_idxM, repmat(bl[,i],1,r[i]), zeros(n_bl,r[i]*(ppC-1)))
    if(i == 1){
      R_con <- kronecker(R_mat,eye(r[i]))
    }else{
      R_con <- as.matrix(Matrix::bdiag(R_con, kronecker(R_mat,eye(r[i]))))
    }
    q_con <- rbind(q_con, zeros(r[i]*size(R_mat,1),1))
  }
  
  # message('EMstep depois segundo loop 1:nb')
  
  bl_idxM <- bl_idxM == 1
  bl_idxQ <- bl_idxQ == 1
  
  #idio
  i_idio_M <- i_idio[1:nM]
  n_idio_M <- length(find(i_idio_M))
  c_i_idio <- cumsum(i_idio)
  
  for(i in 1:n_bl){
    bl_i <- bl[i,]
    rs <- sum(r[bl_i == 1])
    
    idx_i <- NULL
    for(k in 1:nrow(blocks)){
      idx_i[k] <- sum(blocks[k,] == bl_i) == 3
    }
    idx_i <- find(idx_i)
    
    # MONTHLY
    idx_iM <- idx_i[idx_i < (c(nM) + 1)]
    n_i <- length(idx_iM)
    
    denom <- zeros(n_i*rs,n_i*rs)
    nom <- zeros(n_i,rs)
    
    i_idio_i <- i_idio_M[idx_iM] == 1
    i_idio_ii <- c_i_idio[idx_iM]
    i_idio_ii <- i_idio_ii[i_idio_i]
    
    for(t in 1:TT){
      nanYt <- diag(!nanY[idx_iM,t])
      nn <- sum(bl_idxM[i,])
      denom <- denom + kronecker(Zsmooth[bl_idxM[i,],t+1][1:nn] %*% t(Zsmooth[bl_idxM[i,],t+1][1:nn]) + Vsmooth[bl_idxM[i,],bl_idxM[i,],t+1][1:nn,1:nn], nanYt)
      nom <- nom + y[idx_iM,t] %*% t(Zsmooth[bl_idxM[i,],t+1][1:nn]) - nanYt[,i_idio_i] %*% (Zsmooth[rp1+i_idio_ii,t+1] %*% t(Zsmooth[bl_idxM[i,],t+1][1:nn]) + Vsmooth[rp1+i_idio_ii,bl_idxM[i,],t+1][,1:nn])
    }
    
    vec_C <- solve(denom) %*% c(nom)
    C_new[idx_iM,bl_idxM[i,]][,1:nn] <- matlab::reshape(vec_C,n_i,rs)
    
    # QUARTERLY
    idx_iQ <- idx_i[idx_i>c(nM)]
    rps <- rs*ppC
    
    R_con_i <- R_con[,bl_idxQ[i,]]
    q_con_i <- q_con
    no_c <- !(rowSums(R_con_i == 0) == ncol(R_con_i))
    R_con_i <- R_con_i[no_c,]
    q_con_i <- q_con_i[no_c,] 
    
    
    for(j in idx_iQ){
      denom <- zeros(rps,rps)
      nom <- zeros(1,rps)
      idx_jQ <- j-c(nM)
      
      if(i != 1){
        i_idio_jQ <- (rp1+n_idio_M+5*(idx_jQ-1)+1):(rp1+n_idio_M+5*idx_jQ)
        V_0_new[i_idio_jQ,i_idio_jQ] <- Vsmooth[i_idio_jQ,i_idio_jQ,1]
        A_new[i_idio_jQ[1],i_idio_jQ[1]] <- A_i[i_idio_jQ[1]-rp1,i_idio_jQ[1]-rp1]
        Q_new[i_idio_jQ[1],i_idio_jQ[1]] <- Q_i[i_idio_jQ[1]-rp1,i_idio_jQ[1]-rp1]
        
        for(t in 1:TT){
          nanYt <- as.vector(!nanY[j,t])*1
          nn2 <- sum(bl_idxQ[i,])
          denom <- denom + kronecker(Zsmooth[bl_idxQ[i,],t+1][1:nn2] %*% t(Zsmooth[bl_idxQ[i,],t+1][1:nn2]) + Vsmooth[bl_idxQ[i,],bl_idxQ[i,],t+1][1:nn2,1:nn2],nanYt)
          nom <- nom + y[j,t] %*% t(Zsmooth[bl_idxQ[i,],t+1][1:nn2])
          nom <- nom - nanYt %*% (matrix(c(1,2,3,2,1), nrow = 1) %*% Zsmooth[i_idio_jQ,t+1] %*% t(Zsmooth[bl_idxQ[i,],t+1][1:nn2]) +
                                    matrix(c(1,2,3,2,1), nrow = 1) %*% Vsmooth[i_idio_jQ,bl_idxQ[i,],t+1][,1:nn2])
        }
        C_i <- solve(denom) %*% t(nom)
        C_i_constr <- C_i - solve(denom) %*% t(R_con_i) %*% solve(R_con_i %*% solve(denom) %*% t(R_con_i)) %*% (R_con_i %*% C_i - q_con_i)
        nn3 <- sum(bl_idxQ[i,])
        C_new[j,bl_idxQ[i,]][1:nn3] <- C_i_constr
      }
    }
  }
  
  R_new <- zeros(n,n)
  for(t in 1:TT){
    nanYt <- diag(!nanY[,t])*1 == 1
    R_new <- R_new + (y[,t] - nanYt %*% C_new %*% Zsmooth[,t+1]) %*% t(y[,t] - nanYt %*% C_new %*% Zsmooth[,t+1]) + nanYt %*% C_new %*% Vsmooth[,,t+1] %*% t(C_new) %*% nanYt + (eye(n)-nanYt) %*% R %*% (eye(n)-nanYt)
  }
  
  R_new <- R_new/TT
  RR <- diag(R_new) #RR(RR<1e-2) = 1e-2;
  RR[i_idio_M] <- 1e-04
  RR[(nM+1):length(RR)] <- 1e-04
  R_new <- diag(RR)
  
  if(!is.matrix(Z_0)){
    Z_0<-matrix(Z_0,length(Z_0),1)
  }
  
  # output
  return(list(C_new = C_new, R_new = R_new, A_new = A_new, Q_new = Q_new, Z_0 = Z_0, V_0 = V_0, loglik = loglik))
  
}

FIS <- function(Y,Z,R,TT,Q,S){
  
  # library(corpcor)
  
  # %______________________________________________________________________
  # % Fixed intervall smoother (see Harvey, 1989, p. 154)
  # % FIS returns the smoothed state vector AmT and its covar matrix PmT             
  # % Use this in conjunction with function SKF
  # %______________________________________________________________________
  # % INPUT  
  # %        Y         Data                                 (nobs x n)  
  # %        S Estimates from Kalman filter SKF                                                          
  # %          S.Am   : Estimates     a_t|t-1                  (nobs x m) 
  # %          S.Pm   : P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
  # %          S.AmU  : Estimates     a_t|t                    (nobs x m) 
  # %          S.PmU  : P_t|t   = Cov(a_t|t)               (nobs x m x m)       
  # % OUTPUT 
  # %        S Smoothed estimates added to above
  # %          S.AmT  : Estimates     a_t|T                    (nobs x m) 
  # %          S.PmT :  P_t|T   = Cov(a_t|T)               (nobs x m x m)
  # %          S.PmT_1 : Cov(a_ta_t-1|T)
  # %        where m is the dim of state vector and t = 1 ...T is time
  
  m<-dim(S$Am)[1]
  nobs<-dim(S$Am)[2]
  
  S$AmT           = zeros(m,nobs+1)
  S$PmT           = array(0,c(m,m,nobs+1))
  S$AmT[,nobs+1] <- S$AmU[,nobs+1]
  S$PmT[,,nobs+1] <- S$PmU[,,nobs+1]
  
  S$PmT_1<-array(0,c(m,m,nobs))
  S$PmT_1[,,nobs] <- (eye(m)-S$KZ)%*%TT%*%S$PmU[,,nobs]
  
  pinv<-corpcor::pseudoinverse(S$Pm[,,nobs])
  
  J_2 <- S$PmU[,,nobs]%*%t(TT)%*%pinv
  
  for (t in nobs:1){ 
    PmU <- S$PmU[,,t]
    Pm1 <- S$Pm[,,t]
    P_T <- S$PmT[,,t+1]
    P_T1 <- S$PmT_1[,,t]
    
    J_1 <- J_2
    
    S$AmT[,t] <- S$AmU[,t] + J_1%*%(S$AmT[,t+1] - TT%*%S$AmU[,t])
    S$PmT[,,t] <- PmU + J_1%*%(P_T - Pm1)%*%t(J_1) 
    
    if(t>1){
      pinv<-corpcor::pseudoinverse(S$Pm[,,t-1])
      J_2 <- S$PmU[,,t-1]%*%t(TT)%*%pinv
      S$PmT_1[,,t-1] = PmU%*%t(J_2)+J_1%*%(P_T1-TT%*%PmU)%*%t(J_2)
    }
  }
  
  return(S)
  
}

SKF <-function(Y,Z,R,TT,Q,A_0,P_0){
  #Y = y; Z = C; TT = A; A_0 = x_0; P_0 = Sig_0
  
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
      iF  <- solve(Z_t%*%PZ + R_t)
      PZF <- PZ%*%iF
      
      V <- y_t - Z_t%*%A
      Au <- A  + PZF%*%V
      Pu <- P  - PZF%*%t(PZ)
      Pu <-  0.5*(Pu+t(Pu))
      S$loglik <- S$loglik + 0.5*(log(det(iF))  - t(V)%*%iF%*%V)
    }
    
    S$Am[,t] <- A
    S$Pm[,,t] <- P
    
    # Au = A_t|t   & Pu = P_t|t
    
    S$AmU[,t+1] <- Au
    S$PmU[,,t+1] <- Pu
  } # t
  
  if(sum(is.na(y_t))==length(y_t)){
    S$KZ <- zeros(m,m)
  }else{
    S$KZ <- PZF%*%Z_t
  }
  
  return(S)
  
}

MissData <- function(y,C,R){
  
  # library(matlab)
  
  # ______________________________________________________________________
  #  PROC missdata                                                        
  #  PURPOSE: eliminates the rows in y & matrices Z, G that correspond to     
  #           missing data (NaN) in y                                                                                  
  #  INPUT    y             vector of observations at time t  (n x 1 )    
  #           S             KF system matrices             (structure)
  #                       must contain Z & G
  #  OUTPUT   y             vector of observations (reduced)   (# x 1)     
  #           Z G           KF system matrices     (reduced)   (# x ?)     
  #           L             To restore standard dimensions     (n x #)     
  #                         where # is the nr of available data in y
  # ______________________________________________________________________
  
  # [y,C,R,L] 
  
  ix <- !is.na(y)
  e  <- eye(length(y))
  L  <- e[,ix]
  
  y <-    y[ix]
  C  <-  C[ix,]  
  R  <-  R[ix,ix]
  
  return(list(y=y,C=C,R=R,L=L))
  
}

# %%%  Replication files for:
# %%%  ""Nowcasting", 2010, (by Marta Banbura, Domenico Giannone and Lucrezia Reichlin), 
# %%% in Michael P. Clements and David F. Hendry, editors, Oxford Handbook on Economic Forecasting.
# %%%
# %%% The software can be freely used in applications. 
# %%% Users are kindly requested to add acknowledgements to published work and 
# %%% to cite the above reference in any resulting publications
# %--------------------------------------------------------------------------
# % KALMAN FILTER
# %--------------------------------------------------------------------------

runKF <- function(y, A, C, Q, R, x_0, Sig_0){
  
  # x_0 = Z_0; Sig_0 = V_0
  S <- SKF(y,C,R,A,Q, x_0, Sig_0);
  S <- FIS(y,C,R,A,Q,S);
  
  xsmooth <- S$AmT;
  Vsmooth <- S$PmT;
  VVsmooth <- S$PmT_1;
  loglik <- S$loglik;
  
  return(list(xsmooth = xsmooth,Vsmooth = Vsmooth,VVsmooth = VVsmooth,loglik = loglik))  
  
}

em_converged <- function(loglik = NULL, previous_loglik = NULL, threshold = NULL, check_increased = NULL){
  # EM_CONVERGED Has EM converged?
  # [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
  #
  # We have converged if the slope of the log-likelihood function falls below 'threshold',
  # i.e., |f(t) - f(t-1)| / avg < threshold,
  # where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
  # 'threshold' defaults to 1e-4.
  #
  # This stopping criterion is from Numerical Recipes in C p423
  #
  # If we are doing MAP estimation (using priors), the likelihood can decrase,
  # even though the mode of the posterior is increasing.
  
  nargin <- 4 - sum(c(is.null(loglik), is.null(previous_loglik), is.null(threshold), is.null(check_increased)))
  
  if(nargin < 3){threshold <- 1e-4}
  if(nargin < 4){check_increased <- 1}
  
  converged <- 0
  decrease <- 0
  
  if(!is.null(check_increased)){
    if(loglik - previous_loglik < -1e-3){ # allow for a little imprecision
      print(paste(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik))
      decrease <- 1
    }
  }
  
  delta_loglik <- abs(loglik - previous_loglik)
  avg_loglik <- (abs(loglik) + abs(previous_loglik) + 2.2204e-16)/2
  
  # message((delta_loglik / avg_loglik))
  
  if((delta_loglik / avg_loglik) < threshold){converged <- 1}
  
  # output
  list(converged = converged, decrease = decrease)
  
}

