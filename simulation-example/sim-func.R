#####################################################################################
# Codes for simulation studies: 
# The functions for generating simulated data & evaluating performances of competitors,
# which are used to support numerical simulation studies.
#####################################################################################

#####################################################################################
# Functions for generation of simulated data
#####################################################################################
generate.data = function(N, para.list){
  Mu0.list = para.list$Mu0.list
  Theta0.list = para.list$Theta0.list
  Sigma0.list = para.list$Sigma0.list
  
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: generate.data
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the simulated data.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R packages: mvtnorm
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ N: M0 * 1 vector, the sample sizes of subgroups.
  ## @ Mu0.list: a list including M0 mean vectors (p * 1).
  ## @ Theta0.list: a list including M0 precision matrices (p * p).
  ## @ Sigma0.list: a list including M0 correlation matrices (p * p).
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list "whole.data" including:
  ## @ L0: n * 1 vector, the subgroup labels to which each sample belongs.
  ## @ Mu0: M0 * p matrix, M0 mean vectors.
  ## @ Theta0: M0 * p * p array, M0 precision matrices.
  ## @ data: n * p matrix, the design matrix.
  ## @ n_all: int, the total sample size.
  ## @ M0: int, the true number of subgroups.
  ## ---------------------------------------------------------------------------------------------------------------
  
  M0 = length(Mu0.list)
  p = length(Mu0.list[[1]])
  Mu0=matrix(0,M0,p);L0=NULL;Theta0=array(0, dim = c(p, p, M0));data=NULL
  for (k in 1:M0) {
    Mu0[k,] <- Mu0.list[[k]]
    L0 <- c(L0,rep(k,N[k]))
  }
  for (k in 1:M0) {
    Theta0[,,k] <- as.matrix(Theta0.list[[k]])
  }
  for (k in 1:M0) {
    data <- rbind(data,mvrnorm(N[k],Mu0[k,],Sigma0.list[[k]]))
  }
  n_all = dim(data)[1]
  whole.data=list()
  whole.data$L0=L0
  whole.data$Mu0=Mu0
  whole.data$Theta0=Theta0
  whole.data$data=data
  whole.data$n_all=n_all
  whole.data$M0=M0
  return(whole.data)
}

generate.aux.data = function(N0=600, para.list.aux){
  K = length(para.list.aux)
  aux.whole.data = list()
  A.data = list()
  for (k in 1:K) {
    M0k = length(para.list.aux[[k]]$Mu0.list)
    N.k = rep(N0, M0k)
    whole.data.k = generate.data(N.k, para.list.aux[[k]])
    aux.whole.data[[k]] = whole.data.k
    A.data[[k]] = whole.data.k$data
  }
  
  return(list(A.data=A.data,aux.whole.data=aux.whole.data))
  
}

gen.target.para = function(p=100, M0=3, type="power.law", mue=2, nonnum=p/10, 
                           sameSigma=F, seed=NULL){
  
  if(p/nonnum <= M0){
    print("M0 is too large")
    break
  }
  
  Mu0.list = list()
  Sigma0.list = list()
  Theta0.list = list()
  for (k in 1:M0) {
    mu0 = rep(0, p)
    mu0[nonnum*(k-1)+1:nonnum] = mue
    Mu0.list[[k]] = mu0
  }
  
  if(!sameSigma){
    for (k in 1:M0) {
      if(!is.null(seed)){set.seed(p*seed+k)}
      
      if(type=="power.law"){
        Theta0.list[[k]] = Power.law.network0(p,s=10,umin=0.4,umax=0.5)
        Sigma0.list[[k]] = solve(Theta0.list[[k]])
      }
      
      if(type=="tridiag"){
        Theta0.list[[k]] = tridiag.cor(p, s=1, rho=min(k,4.5)/10)
        Sigma0.list[[k]] = solve(Theta0.list[[k]])
      }
      
    }
  } else {
    if(!is.null(seed)){set.seed(p*seed)}
    
    if(type=="power.law"){
      Theta0 = Power.law.network0(p,s=10,umin=0.4,umax=0.5)
      Sigma0 = solve(Theta0)
    }
    
    if(type=="tridiag"){
      Theta0 = tridiag.cor(p, s=1, rho=0.3)
      Sigma0 = solve(Theta0)
    }
    
    for (k in 1:M0) {
        Theta0.list[[k]] = Theta0
        Sigma0.list[[k]] = Sigma0
    }
    
  }
  
  
  return(list(Mu0.list=Mu0.list, Theta0.list=Theta0.list, Sigma0.list=Sigma0.list, M0=M0))
  
  
}

gen.aux.para = function(n=200, p=100, para.list.target, K=6, K.A=3, M0.A.vec=c(2,3,4,3,3,3), 
                        c0=0.1, prob0=0.1, mue=2, nonnum=p/10 ){
  if(length(M0.A.vec) < K){ print("Error: the length of M0.A.vec < K"); break }
  para.aux = list()
  Mu0.list = para.list.target$Mu0.list
  Sigma0.list = para.list.target$Sigma0.list
  M0 = length(Sigma0.list)
  for (k in 1:K) {
    if(k <= K.A){
      M0k = M0.A.vec[k]
      Mu0.list.k = list()
      Sigma0.list.k = list()
      Theta0.list.k = list()
      for (m in 1:M0k) {
        # mu0 = rep(0, p)
        # mu0[nonnum*(m-1)+1:nonnum] = mue
        # Mu0.list.k[[m]] = mu0
        if(m==k){
          uni.up = c0 * sqrt( log(p) / n )
          Mu0.list.k[[m]] = Mu0.list[[m]] + rbinom(p,size=1,prob=prob0)*runif(p,-uni.up,uni.up)
          # Theta0.list.k[[m]] = Theta0.list[[m]]
          Sigma0 = Sigma0.list[[m]]
          delta.k = matrix(rbinom(p^2,size=1,prob=prob0)*runif(p^2,-uni.up,uni.up), ncol=p)
          Sig.k = Sigma0 %*% (delta.k + diag(1,p))
          Sig.k = (Sig.k+t(Sig.k))/2
          if(min(eigen(Sig.k)$values)<0.05){
            Sig.k = Sig.k + diag(0.1-min(eigen(Sig.k)$values),p)
          }
          Sigma0.list.k[[m]] = Sig.k
          Theta0.list.k[[m]] = solve(Sig.k)
          
        } else {
          mu0 = rep(0, p)
          mu0[nonnum*(m-1)+1:nonnum] = 2*mue
          Mu0.list.k[[m]] = mu0
          
          # Theta0.list.k[[m]] = diag(p)
          # Sigma0.list.k[[m]] = diag(p)
          Theta0.list.k[[m]] = Power.law.network0(p,s=5,umin=0.4,umax=0.5)
          Sigma0.list.k[[m]] = solve(Theta0.list.k[[m]])
        }
      }
      para.k = list(Mu0.list=Mu0.list.k, Theta0.list=Theta0.list.k, Sigma0.list=Sigma0.list.k)
      para.aux[[k]] = para.k
    } else{
      KK = sum(M0.A.vec[1:K.A])
      diffe.t = matrix(10,KK,M0)
      for (m in 1:M0) {
        v = 1
        for (kk in 1:K.A) {
          for (mk in 1:M0.A.vec[kk]) {
            diffe.t[v,m] = sum((para.list.target$Theta0.list[[m]] %*% para.aux[[kk]]$Sigma0.list[[mk]] - diag(p))^2)
            v = v+1
          }
        }
      }
      
      v = which.max(apply(diffe.t, 1, mean) / apply(diffe.t, 1, sd))
      M0.A.veckk = c(0,M0.A.vec)
      for (kk in 1:K.A) {
        if((sum(M0.A.veckk[1:kk]) < v) & (sum(M0.A.veckk[1:(kk+1)]) >= v)){
          k0 = kk
          mk0 = v - sum(M0.A.veckk[1:kk])
        }
      }
      
      M0k = M0.A.vec[k]
      Mu0.list.k = list()
      Sigma0.list.k = list()
      Theta0.list.k = list()
      for (m in 1:M0k) {
        mu0 = rep(0, p)
        mu0[nonnum*(m-1)+1:nonnum] = 2*mue
        Mu0.list.k[[m]] = mu0
        Theta0.list.k[[m]] = para.aux[[k0]]$Theta0.list[[mk0]]
        Sigma0.list.k[[m]] = para.aux[[k0]]$Sigma0.list[[mk0]]
      }
      # M0k = M0.A.vec[k]
      # Mu0.list.k = list()
      # Sigma0.list.k = list()
      # Theta0.list.k = list()
      # for (m in 1:M0k) {
      #   mu0 = rep(0, p)
      #   mu0[nonnum*(m-1)+1:nonnum] = mue
      #   Mu0.list.k[[m]] = mu0
      #   # Theta0.list.k[[m]] = diag(p)
      #   # Sigma0.list.k[[m]] = diag(p)
      #   Theta0.list.k[[m]] = Power.law.network0(p,s=1,umin=0.5,umax=0.75)
      #   Sigma0.list.k[[m]] = solve(Theta0.list.k[[m]])
      # }
      para.k = list(Mu0.list=Mu0.list.k, Theta0.list=Theta0.list.k, Sigma0.list=Sigma0.list.k)
      para.aux[[k]] = para.k
    }
    
  }
  
  return(para.aux)
}

gen.aux.para.overall = function(n=200, p=100, para.list.target, K=5, K.A=3, M0.A.vec=rep(para.list.target$M0,K), 
                        c0=0.1, prob0=0.1, mue=2, nonnum=p/10 ){
  # auxilary domains: overall similarity or non-informative completely
  if(length(M0.A.vec) < K){ print("Error: the length of M0.A.vec < K"); break }
  para.aux = list()
  Mu0.list = para.list.target$Mu0.list
  Sigma0.list = para.list.target$Sigma0.list
  M0 = length(Sigma0.list)
  for (k in 1:K) {
    if(k <= K.A){
      M0k = M0.A.vec[k]
      Mu0.list.k = list()
      Sigma0.list.k = list()
      Theta0.list.k = list()
      
      for (m in 1:M0k) {
        uni.up = c0 * sqrt( log(p) / n )
        Mu0.list.k[[m]] = Mu0.list[[min(m,M0)]] + rbinom(p,size=1,prob=prob0)*runif(p,-uni.up,uni.up)
        
        Sigma0 = Sigma0.list[[min(m,M0)]]
        set.seed(2025)
        delta.k = matrix(rbinom(p^2,size=1,prob=prob0)*runif(p^2,-uni.up,uni.up), ncol=p)
        Sig.k = Sigma0 %*% (delta.k + diag(1,p))
        Sig.k = (Sig.k+t(Sig.k))/2
        if(min(eigen(Sig.k)$values)<0.05){
          Sig.k = Sig.k + diag(0.1-min(eigen(Sig.k)$values),p)
        }
    
        Sigma0.list.k[[m]] = Sig.k
        Theta0.list.k[[m]] = solve(Sig.k)
      }
      
      para.k = list(Mu0.list=Mu0.list.k, Theta0.list=Theta0.list.k, Sigma0.list=Sigma0.list.k)
      para.aux[[k]] = para.k
    } else{
     
      M0k = M0.A.vec[k]
      Mu0.list.k = list()
      Sigma0.list.k = list()
      Theta0.list.k = list()
      
      for (m in 1:M0k) {
        mu0 = rep(0, p)
        mu0[nonnum*(m-1)+1:nonnum] = 2*mue
        Mu0.list.k[[m]] = mu0
        
        Sigma0 = Sigma0.list[[min(m,M0)]]
        uni.up = 100 * c0 * sqrt( log(p) / n )
        set.seed(2025)
        delta.k = matrix(rbinom(p^2,size=1,prob=prob0)*runif(p^2,-uni.up,uni.up), ncol=p)
        Sig.k = Sigma0 %*% (delta.k + diag(1,p))
        Sig.k = (Sig.k+t(Sig.k))/2
        if(min(eigen(Sig.k)$values)<0.05){
          Sig.k = Sig.k + diag(0.1-min(eigen(Sig.k)$values),p)
        }
        
        Sigma0.list.k[[m]] = Sig.k
        Theta0.list.k[[m]] = solve(Sig.k)
      }
      

      para.k = list(Mu0.list=Mu0.list.k, Theta0.list=Theta0.list.k, Sigma0.list=Sigma0.list.k)
      para.aux[[k]] = para.k
    }
    
  }
  
  return(para.aux)
}

Power.law.network0 = function(p,s=10,umin=0.4,umax=0.5){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: Power.law.network
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the s-block power-law precision matrices.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding data: No
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ p: Dimensions of the precision matrix.
  ## @ s: The number of sub-networks.
  ## @ umin: The lower bound of non-zero elements on non-diagonal elements.
  ## @ umax: The upper bound of non-zero elements on non-diagonal elements.
  ## @ I2: The replacement blocks for the precision matrix of the second subgroup.
  ## @ I3: The replacement blocks for the precision matrix of the third subgroup.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ A list including The precision matrices of three subgroups.
  ## ---------------------------------------------------------------------------------------------------------------
  
  pp=p/s
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  submatrix=list()
  for (ss in 1:s) {
    g = sample_pa(pp, m=2,directed = F)
    Eg = as.data.frame(as_edgelist(g))
    subi = diag(1,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }
    for (i in 1:pp) {
      subi[i,i] = sum(abs(subi[i,setdiff(1:pp,i)]))+0.1
    }
    submatrix[[ss]]=subi
  }
  
  A=submatrix[[1]]
    if(s > 1){
      for (ss in 2:s) {
      A=bdiag(A,submatrix[[ss]])
    }
  }
  A = as.matrix(A)
  
  return(A)
}

tridiag.cor = function(p, s=1, rho=0.3){
  
  ## ---------------------------------------------------------------------------------------------------------------
  ## The name of the function: tridiag.cor
  ## ---------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Generating the tri-diagonal precision matrices.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Required preceding data: No
  ## ---------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ p: Dimensions of the precision matrix.
  ## @ s: The number of sub-networks.
  ## @ rho: The value of non-zero elements on non-diagonal elements.
  ## ---------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ A: The precision matrix.
  ## ---------------------------------------------------------------------------------------------------------------
  
  m=p/s
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  sig.G = matrix(0,m,m)
  for (j in 1:m) {
    ncol0 = c(j-1,j,j+1)
    sig = c(rho, 1, rho)
    ncol = ncol0[which(ncol0>0 & ncol0<=m)]
    sig = sig[which(ncol0>0 & ncol0<=m)]
    sig.G[j,ncol] = sig
  }
  submatrix=list()
  for (ss in 1:s) {
    submatrix[[ss]] = sig.G
  }
  A=submatrix[[1]]
  if(s > 1){
    for (ss in 2:s) {
      A=bdiag(A,submatrix[[ss]])
    }
  }
  return(A)
}


################################################################################
# Functions for evaluating performances of proposed methods and alternatives:
################################################################################
Evaluation.GGMM = function(data, mu_hat, Theta_hat, Mu0, Theta0, M0, L.mat, L0, prob){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: Esti.error
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Gauging performance of the proposed and alternative approaches
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: 
  ##            R functions: f.den.vec() 
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: ( (q+1)*(p+1)-1 ) * s matrix, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ mu_hat: M0_hat * p matrix, the estimated mean vectors of M0_hat subgroups.
  ## @ Theta_hat: p * p * M0_hat array, the estimated precision matrices of M0_hat subgroups.
  ## @ other input parameters: Similar to the previous definition.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## The vector including:
  ## @ K: The estimated number of subgroups.
  ## @ CE: The sub-grouping error
  ## @ CME: The mean squared error (MSE) for the mean vectors.
  ## @ PME: The mean squared error (MSE) for the precision matrices.
  ## @ TPR/FPR: The true and false positive rates for the off-diagonal elements of the precision matrices.
  ## -----------------------------------------------------------------------------------------------------------------
  
  p = dim(mu_hat)[2]
  K_hat = dim(mu_hat)[1]
  n_all = dim(data)[1]
  if(K_hat == M0){
    num = rep(0,M0)
    numk = NULL
    for (k in 1:M0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    mu_hat = mu_hat[num,]
    Theta_hat = Theta_hat[,,num]
    #  CME & PME
    CME = sqrt(sum((mu_hat - Mu0)^2))/M0
    PME = sqrt(sum((Theta_hat - Theta0)^2))/M0
    
    #  TPR & FPR
    TPk = (apply((Theta_hat!=0) + (Theta0!=0) == 2,3,sum) - p) / (apply(Theta0!=0,3,sum) - p)
    TPk[(is.na(TPk))] = 1
    TPR = sum(TPk) / M0
    FPk = apply((Theta_hat!=0) + (Theta0==0) == 2,3,sum) / apply(Theta0==0,3,sum)
    FPk[(is.na(FPk))] = 0
    FPR = sum(FPk[(!is.na(FPk))]) / M0
    
    #  CE
    f.mat = matrix(0,n_all,M0)
    L.mat = matrix(0,n_all,M0)
    for(k.ind in 1:K_hat) {
      f.mat[,k.ind]=f.den.vec( data, as.numeric(mu_hat[k.ind,]), Theta_hat[,,k.ind] )                  
    }
    for(k.ind in 1:M0) {
      for(i in 1:n_all) {
        L.mat[i,k.ind] = prob[k.ind] * f.mat[i,k.ind] / prob %*% f.mat[i,]
      }
    }
    member = apply(L.mat,1,which.max)
    
    aa = L0
    cap_matrix0 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = member
    cap_matrix1 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)
  } else{
    num = rep(0,K_hat)
    for (k in 1:K_hat) {
      mu_hatk = mu_hat[k,]
      errork.mu = apply((t(Mu0) - mu_hatk)^2,2,sum)
      Theta_hatk = Theta_hat[,,k]
      errork.Theta = apply((Theta0 - rep(Theta_hatk,M0))^2,3,sum)
      errork = errork.mu + errork.Theta
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    Mu0.re = Mu0[num,]
    Theta0.re = Theta0[,,num]
    if(K_hat == 1){
      Mu0.re = as.matrix(t(Mu0.re))
      Theta0.re = as.array(Theta0.re)
      dim(Theta0.re) <- c(p,p,K_hat)
    }
    
    #  CME & PME
    CME = sqrt(sum((mu_hat - Mu0.re)^2))/K_hat
    PME = sqrt(sum((Theta_hat - Theta0.re)^2))/K_hat
    
    #  TPR & FPR
    TPk = (apply((Theta_hat!=0) + (Theta0.re!=0) == 2,3,sum) - p) / (apply(Theta0.re!=0,3,sum) - p)
    TPk[(is.na(TPk))] = 1
    TPR = sum(TPk) / K_hat
    FPk = apply((Theta_hat!=0) + (Theta0.re==0) == 2,3,sum) / apply(Theta0.re==0,3,sum)
    FPk[(is.na(FPk))] = 0
    FPR = sum(FPk[(!is.na(FPk))]) / K_hat
    
    #  CE
    num = rep(0,M0)
    numk = NULL
    L.hat = rep(0,n_all)
    for (k in 1:M0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
      L.hat[which(L0 == k)] = numk
    }
    if(K_hat < M0){L.hat=L0}
    
    f.mat = matrix(0,n_all,K_hat)
    L.mat = matrix(0,n_all,K_hat)
    for(k.ind in 1:K_hat) {
      f.mat[,k.ind]=f.den.vec( data, as.numeric(mu_hat[k.ind,]), Theta_hat[,,k.ind] )                  
    }
    for(k.ind in 1:K_hat) {
      for(i in 1:n_all) {
        L.mat[i,k.ind] = prob[k.ind] * f.mat[i,k.ind] / prob %*% f.mat[i,]
      }
    }
    member = apply(L.mat,1,which.max)
    
    aa = L.hat
    cap_matrix0 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = member
    cap_matrix1 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)
  }
  index = as.data.frame(t(c(as.numeric(K_hat == M0),K_hat,CE,CME,PME,TPR,FPR)))
  names(index) = c("Per","K","CE","CME","PME","TPR","FPR")
  
  return(index)
}

Esti.error.JGL = function(data, mu_hat, Theta_hat, Mu0, Theta0, M0, L0, member){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: Esti.error.JGL
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##      Gauging performance of the K-means + JGL 
  ##      (compared to the function Esti.error(), the input is slightly different due to the specificity of K-means).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: 
  ##            R functions: f.den.vec() 
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: ( (q+1)*(p+1)-1 ) * s matrix, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ mu_hat: M0_hat * p matrix, the estimated mean vectors of M0_hat subgroups.
  ## @ Theta_hat: p * p * M0_hat array, the estimated precision matrices of M0_hat subgroups.
  ## @ other input parameters: Similar to the previous definition.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## The vector including:
  ## @ K: The estimated number of subgroups.
  ## @ CE: The sub-grouping error
  ## @ CME: The mean squared error (MSE) for the mean vectors.
  ## @ PME: The mean squared error (MSE) for the precision matrices.
  ## @ TPR/FPR: The true and false positive rates for the off-diagonal elements of the precision matrices.
  ## -----------------------------------------------------------------------------------------------------------------
  
  p = dim(mu_hat)[2]
  K_hat = dim(mu_hat)[1]
  n_all = dim(data)[1]
  if(K_hat == M0){
    num = rep(0,M0)
    numk = NULL
    L.hat = rep(0,n_all)
    for (k in 1:M0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
      L.hat[which(L0 == k)] = numk
    }
    mu_hat = mu_hat[num,]
    Theta_hat = Theta_hat[,,num]
    #  CME & PME
    CME = sqrt(sum((mu_hat - Mu0)^2))/M0
    PME = sqrt(sum((Theta_hat - Theta0)^2))/M0
    
    #  TPR & FPR
    TPk = (apply((Theta_hat!=0) + (Theta0!=0) == 2,3,sum) - p) / (apply(Theta0!=0,3,sum) - p)
    TPk[(is.na(TPk))] = 1
    TPR = sum(TPk) / M0
    FPk = apply((Theta_hat!=0) + (Theta0==0) == 2,3,sum) / apply(Theta0==0,3,sum)
    FPk[(is.na(FPk))] = 0
    FPR = sum(FPk[(!is.na(FPk))]) / M0
    
    #  CE
    aa = L.hat
    cap_matrix0 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = member
    cap_matrix1 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)

    index = as.data.frame(t(c(K_hat,CE,CME,PME,TPR,FPR)))
    names(index) = c("K","CE","CME","PME","TPR","FPR")
    
    return(index)
  } else{
    num = rep(0,K_hat)
    for (k in 1:K_hat) {
      mu_hatk = mu_hat[k,]
      errork.mu = apply((t(Mu0) - mu_hatk)^2,2,sum)
      Theta_hatk = Theta_hat[,,k]
      errork.Theta = apply((Theta0 - rep(Theta_hatk,M0))^2,3,sum)
      errork = errork.mu + errork.Theta
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    Mu0.re = Mu0[num,]
    Theta0.re = Theta0[,,num]
    if(K_hat == 1){
      Mu0.re = as.matrix(t(Mu0.re))
      Theta0.re = as.array(Theta0.re)
      dim(Theta0.re) <- c(p,p,K_hat)
    }
    
    #  CME & PME
    CME = sqrt(sum((mu_hat - Mu0.re)^2))/K_hat
    PME = sqrt(sum((Theta_hat - Theta0.re)^2))/K_hat
    
    #  TPR & FPR
    TPk = (apply((Theta_hat!=0) + (Theta0.re!=0) == 2,3,sum) - p) / (apply(Theta0.re!=0,3,sum) - p)
    TPk[(is.na(TPk))] = 1
    TPR = sum(TPk) / K_hat
    FPk = apply((Theta_hat!=0) + (Theta0.re==0) == 2,3,sum) / apply(Theta0.re==0,3,sum)
    FPk[(is.na(FPk))] = 0
    FPR = sum(FPk[(!is.na(FPk))]) / K_hat
    
    #  CE
    num = rep(0,M0)
    numk = NULL
    L.hat = rep(0,n_all)
    for (k in 1:M0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
      L.hat[which(L0 == k)] = numk
    }
    if(K_hat < M0){L.hat=L0}
    if(K_hat == 2){
      m12=member[L0==1]
      if(sum(m12==1)<sum(m12==2)){
        n.L.hat = c(which(L.hat==1),which(L.hat==2))
        L.hat[n.L.hat]=c(L0[L0==2],L0[L0==1])
      }
    }
    aa = L.hat
    cap_matrix0 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = member
    cap_matrix1 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)
    
    index = as.data.frame(t(c(as.numeric(K_hat == M0),K_hat,CE,CME,PME,TPR,FPR)))
    names(index) = c("Per","K","CE","CME","PME","TPR","FPR")
    
    return(index)
  }
}

Evaluation.mtlgmm = function(data, mu_hat, Theta_hat, Mu0, Theta0, M0, L.mat, L0){
  
  p = dim(mu_hat)[2]
  K_hat = dim(mu_hat)[1]
  n_all = dim(data)[1]
  
  if(K_hat == M0){
    num = rep(0,M0)
    numk = NULL
    for (k in 1:M0) {
      errork = apply((t(mu_hat) - Mu0[k,])^2,2,sum)
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    mu_hat = mu_hat[num,]
    Theta_hat = Theta_hat[,,num]
    #  CME & PME
    CME = sqrt(sum((mu_hat - Mu0)^2))/M0
    PME = sqrt(sum((Theta_hat - Theta0)^2))/M0
    
    #  TPR & FPR
    TPk = (apply((Theta_hat!=0) + (Theta0!=0) == 2,3,sum) - p) / (apply(Theta0!=0,3,sum) - p)
    TPk[(is.na(TPk))] = 1
    TPR = sum(TPk) / M0
    FPk = apply((Theta_hat!=0) + (Theta0==0) == 2,3,sum) / apply(Theta0==0,3,sum)
    FPk[(is.na(FPk))] = 0
    FPR = sum(FPk[(!is.na(FPk))]) / M0
    
    #  CE
    member = L.mat
    
    aa = L0
    cap_matrix0 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = member
    cap_matrix1 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)
  } else{
    num = rep(0,K_hat)
    for (k in 1:K_hat) {
      mu_hatk = mu_hat[k,]
      errork.mu = apply((t(Mu0) - mu_hatk)^2,2,sum)
      Theta_hatk = Theta_hat[,,k]
      errork.Theta = apply((Theta0 - rep(Theta_hatk,M0))^2,3,sum)
      errork = errork.mu + errork.Theta
      numk = which(errork == min(errork))[1]
      num[k] = numk
    }
    Mu0.re = Mu0[num,]
    Theta0.re = Theta0[,,num]
    if(K_hat == 1){
      Mu0.re = as.matrix(t(Mu0.re))
      Theta0.re = as.array(Theta0.re)
      dim(Theta0.re) <- c(p,p,K_hat)
    }
    
    #  CME & PME
    CME = sqrt(sum((mu_hat - Mu0.re)^2))/K_hat
    PME = sqrt(sum((Theta_hat - Theta0.re)^2))/K_hat
    
    #  TPR & FPR
    TPk = (apply((Theta_hat!=0) + (Theta0.re!=0) == 2,3,sum) - p) / (apply(Theta0.re!=0,3,sum) - p)
    TPk[(is.na(TPk))] = 1
    TPR = sum(TPk) / K_hat
    FPk = apply((Theta_hat!=0) + (Theta0.re==0) == 2,3,sum) / apply(Theta0.re==0,3,sum)
    FPk[(is.na(FPk))] = 0
    FPR = sum(FPk[(!is.na(FPk))]) / K_hat
    
    #  CE
    member = L.mat
    
    aa = L0
    cap_matrix0 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix0[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    aa = member
    cap_matrix1 = matrix(0,n_all,n_all)
    for(i in 1:(n_all-1)){
      for (j in (i+1):(n_all)) {
        cap_matrix1[i,j] <- as.numeric(aa[i] == aa[j])
      }
    }
    CE = sum(abs(cap_matrix1-cap_matrix0)) / (n_all*(n_all-1)/2)
  }
  index = as.data.frame(t(c(as.numeric(K_hat == M0),K_hat,CE,CME,PME,TPR,FPR)))
  names(index) = c("Per","K","CE","CME","PME","TPR","FPR")
  
  return(index)
}


##################### Summary resulds of alternative methods  #####################
cha.res = function(all_store, only.mean=NULL){
  mean.vec = as.character(format(round(apply(all_store, 2, mean), digits = 4), nsmall = 4))
  sd.vec = as.character(format(round(apply(all_store, 2, sd), digits = 4), nsmall = 4))
  res = paste0(mean.vec,"(",sd.vec,")")
  if(!is.null(only.mean)){
    res[only.mean] = mean.vec[only.mean]
  }
  return(res)
}

