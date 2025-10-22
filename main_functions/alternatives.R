# -------------- The SCAN method (JMLR, 2018)----------------
SCAN.hao = function(data, K, lambda1 = 0.1, lambda2 = 0.1, lambda3 = 0.1,
                    eps = 1e-2, niter = 20, initial = T, initialize, concave=F){

  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: SCAN.hao
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Implementing the SCAN method (Hao et al., 2018), and the main contributor to this function is Botao Hao.
  ##            Hao, B., Sun, W. W., Liu, Y., & Cheng, G. (2018). Simultaneous clustering and estimation of heterogeneous graphical models.
  ##            The Journal of Machine Learning Research, 18(1), 7981-8038.
  ##            https://www.jmlr.org/papers/v18/17-019.html
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: f.den.vec()    Symmetrize()
  ##            R packages: JGL
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ K: int, a selected K.
  ## @ lambda1: a float value, the tuning parameter controlling the sparse of the mean parameter.
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the similarity of subgroups.
  ## @ eps: a float value, algorithm termination threshold.
  ## @ niter: int, Maximum number of cycles of the algorithm.
  ## @ initialization: the logical variable, whether to calculate the initial value, the default setting is T.
  ##                                     if initialization = F, the initial value uses initialize.
  ## @ initialize: A given initial value used if initialization = F.
  ## @ concave: a logical variable, whether to use the concave penalty on the mean parameters, the default setting is F.
  ##                                If concave = T, it can bring a big boost to the speed of the algorithm.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including:
  ## @ mu: K * p matrix, the estimated mean vectors of K subgroups.
  ## @ Omega: p * p * K array, the estimated precision matrices of K subgroups.
  ## @ prob: K * 1 vector, the estimated mixture probabilities of subgroups.
  ## @ L.mat: n * K matrix, the estimated probability that each sample belongs to each subgroup.
  ## @ bic: a float value, the BIC value corresponding the choice of given tuning parameters.
  ## @ fiterror: a float value, the value of the loss function (without penalty function) corresponding the choice of given tuning parameters.
  ## @ df: a float value, the penalty value for non-zero parameters corresponding the choice of given tuning parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------

  n <- dim(data)[1]
  p <- dim(data)[2]

  if(initial){
    out.initial = myinitial(data,K)
    prob = out.initial$prob
    mu = out.initial$Mu
    S = out.initial$S
    memb = out.initial$memb
    L.mat = matrix(0,n,K)
    for(jj in 1:n) L.mat[jj, memb[jj]]=1
  } else {
    mu = initialize$Mu
    S = initialize$S
    prob = initialize$prob
    memb = initialize$memb
    L.mat = matrix(0,n,K)
    for(jj in 1:n) L.mat[jj, memb[jj]]=1
  }

  f.mat = matrix(0,n,K)
  L.mat.old = matrix(0,n,K)
  mu.old = matrix(0, K, p)
  S.old = array(0, dim = c(p, p, K))
  Omega.old = array(0, dim = c(p, p, K))
  Omega = array(0, dim = c(p, p, K))
  member.old = rep(1,n)

  for(k.ind in 1:K){
    cond.no = function(delta){
      tempOmega = S[,,k.ind] +diag(delta,p)
      sing.values = svd(tempOmega)$d
      M = max(sing.values)
      m = min(sing.values)
      return(M/m-p)
    }

    lowbd = -1
    ttt = 0
    while(cond.no(lowbd) * cond.no(10) > 0 && ttt <=10){
      lowbd = lowbd/2
      ttt = ttt + 1
    }
    delta = uniroot(cond.no,interval=c(lowbd,10),tol = 0.00000001)$root
    S[,,k.ind] = S[,,k.ind] + diag(delta,p)
    Omega[,,k.ind] = solve(S[,,k.ind])
  }
  # avoid corner case in the following while loop.
  if(is.na(f.density(as.numeric(data[1,]), mu[1,], Omega[,,1])) == TRUE || f.density(as.numeric(data[1,]), mu[1,], Omega[,,1]) == Inf){
    for(k.ind in 1:K){
      Omega[,,k.ind] = diag(p)
    }
  }
  t = 0
  while((norm(mu.old-mu,type="2")/(norm(mu,type="2")+0.00001) >= eps ||
         norm(Omega.old-Omega,type="2")/(norm(Omega,type="2")+0.00001) >= eps) &&
        t < niter)
  {

    t = t + 1
    prob.old = prob
    mu.old = mu
    Omega.old = Omega
    L.mat.old = L.mat

    for(k.ind in 1:K) {
      f.mat[,k.ind]=f.den.vec( data, as.numeric(mu.old[k.ind,]), Omega.old[,,k.ind] )
    }

    for(k.ind in 1:K) {
      for(i in 1:n) {
        L.mat[i,k.ind] = prob.old[k.ind] * f.mat[i,k.ind] / prob.old %*% f.mat[i,]
      }
      prob[k.ind] = mean(L.mat[,k.ind])
    }

    if(!concave){
      for(j in 1:p){
        for(k.ind in 1:K) {
          tmp = rep(0,n)
          for(i in 1:n){
            tmp[i] = t(data[i,] - mu.old[k.ind,]) %*% Omega.old[,j,k.ind] + mu.old[k.ind,j] * Omega.old[j,j,k.ind]
          }
          if(n*lambda1 >= abs(t(L.mat[,k.ind])%*%tmp)){
            mu[k.ind,j] = 0
          }else{
            tmp2 = rep(0,n)
            for(i in 1:n){
              tmp2[i] = t(data[i,]) %*% Omega.old[,j,k.ind]
            }
            mu[k.ind,j] = (Omega.old[j,j,k.ind]*sum(L.mat[,k.ind]))^(-1) * as.numeric(  t(L.mat[,k.ind]) %*% tmp2 - sum(L.mat[,k.ind])*( t(mu.old[k.ind,]) %*% Omega.old[,j,k.ind] - mu.old[k.ind,j]*Omega.old[j,j,k.ind]) - n*lambda1*sign(mu.old[k.ind,j])  )
          }
        }
      }
    } else {
      nK = apply(L.mat,2,sum)
      for(j in 1:p){
        for(k.ind in 1:K) {
          tmp = t(t(data) - as.numeric(mu[k.ind,])) %*% Omega.old[,j,k.ind] + mu[k.ind,j] * Omega.old[j,j,k.ind]
          a=3
          hj = t(L.mat[,k.ind])%*%tmp
          tau_k = sqrt(apply((t(mu[k.ind,] - t(mu[-k.ind,])))^2,1,sum))
          v_k = sum(mcp_d(tau_k, lambda3, a) / (tau_k+0.00001) * mu[-k.ind,j])
          mcp_lambda1 = mcp_d(mu[k.ind,j], lambda1, a)
          v_k_hat = sum(mcp_d(tau_k, lambda3, a) / (tau_k+0.00001))

          mu[k.ind,j] = (hj + n*v_k) / (nK[k.ind] * Omega.old[j,j,k.ind] + n*v_k_hat + n*mcp_lambda1/(abs(mu[k.ind,j])+0.00001))
        }
      }
      mu[abs(mu) < 1e-3] <- 0
    }

    cluster.data = vector("list", length = K)
    member = apply(L.mat,1,which.max)
    if(length(unique(member))<K){
      # print("Warning: combined clusters!")
      break
    }
    breakflag = 0
    for(k in 1:K) {
      cluster.data[[k]] = data[member==k,]
      if(is.null(dim(cluster.data[[k]])) == TRUE){
        breakflag = 1
        break
      }
    }
    if(breakflag == 0){
      out.jgl <- JGL(Y=cluster.data, penalty="group",lambda2,lambda3,return.whole.theta=TRUE)
      for(k in 1:K) {
        Omega[,,k] = out.jgl$theta[[k]]
      }
    }else{
      break
    }
  }
  for(k in 1:K) {
    Omega[,,k] = Symmetrize(Omega[,,k])
  }
  #compute BIC
  if(length(unique(member))<K || norm(mu,type="2")==0){
    out.BIC = list()
    out.BIC$bic = 1e10
    out.BIC$fiterror = 1e10
    out.BIC$df = 1e10
  }else{
    out.BIC = BIC(data, mu, Omega, L.mat)
  }
  P <- list()
  P$member <-  member
  P$mu <-  mu
  P$Omega <- Omega
  P$prob <- prob
  P$L.mat <- L.mat
  P$bic <- out.BIC$bic
  P$fiterror <- out.BIC$fit.error
  P$df <- out.BIC$df
  return(P)
}

tuning.lambda.SCAN = function(lambda, data, K, initial.selection="K-means",
                              initialize, trace = F, concave=F){
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: tuning.lambda.FGGM
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Searching and selecting the optional tuning parameters under the adaptive BIC-type criterion
  ##            using the SCAN method.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: SCAN.hao()
  ##            R packages:  glasso  mnormt  JGL MASS  caret Matrix  SoDA
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ lambda: a list, the sequences of the tuning parameters (lambda1, lambda2, and lambda3).
  ## @ data: n * p matrix, the design matrix.
  ## @ K: int, a selected upper bound of K_0.
  ## @ initial.selection: the different initial values from two clustering methods, which can be selected from c("K-means","dbscan").
  ## @ initialize: A given initial values, which should be given when initial.selection is not in c("K-means","dbscan").
  ## @ trace: the logical variable, whether or not to output the number of identified subgroups during the search for parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list "result" including:
  ## @ Opt_lambda: the selected optional tuning parameters.
  ## @ Mu_hat.list: the estimated mean vectors of K0_hat subgroups corresponding all choices of given tuning parameters.
  ## @ Theta_hat.list: the estimated precision matrices of K0_hat subgroups corresponding all choices of given tuning parameters.
  ## @ prob.list: the estimated mixture probabilities of subgroups corresponding all choices of given tuning parameters.
  ## @ member.list: subgroup labels to which each sample belongs corresponding all choices of given tuning parameters.
  ## @ L.mat.list: the estimated probability that each sample belongs to each subgroup corresponding all choices of given tuning parameters.
  ## @ Opt_aBIC: the optional BIC value.
  ## @ Opt_num: the position of the optimal parameter for all given parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  t0 = proc.time()
  
  lambda1 = lambda$lambda1
  lambda2 = lambda$lambda2
  lambda3 = lambda$lambda3
  L1 = length(lambda1)
  L2 = length(lambda2)
  L3 = length(lambda3)
  L = L1+L2+L3

  aBIC = rep(0,L)
  n_all = dim(data)[1]
  # initialize
  if(initial.selection=="K-means"){
    out.initial = initialize_fuc(data,K)
  } else {out.initial = initialize}

  if(L == 3){
    l=1
    index = as.data.frame(matrix(0,l,length(nameindex)))
    names(index) = nameindex
    aBIC = rep(0,l)
    Mu_hat.list = list()
    Theta_hat.list = list()
    prob.list = list()
    L.mat.list = list()
    member.list = list()
    lam1 = lambda1;lam2 = lambda2;lam3 = lambda3;
    SCAN = SCAN.hao(data, K, lam1, lam2, lam3, initial = F, initialize=out.initial, concave=concave)
    mu_hat <- SCAN$mu; Theta_hat <- SCAN$Omega; prob=SCAN$prob; L.mat=SCAN$L.mat
    aBIC[l] = SCAN$bic; member = SCAN$member
    Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; prob.list[[l]]=prob;L.mat.list[[l]] = L.mat;member.list[[l]]=member
  } else {
    Mu_hat.list = list()
    Theta_hat.list = list()
    prob.list = list()
    L.mat.list = list()
    member.list = list()
    lam1 = median(lambda1)
    # search lam3
    lam1 = median(lambda1);lam2 = median(lambda2)
    for (l in 1:L3) {
      lam3 = lambda3[l]
      SCAN = SCAN.hao(data, K, lam1, lam2, lam3, initial = F, initialize=out.initial, concave=concave)
      mu_hat <- SCAN$mu; Theta_hat <- SCAN$Omega; prob=SCAN$prob; L.mat=SCAN$L.mat
      aBIC[l] = SCAN$bic; member = SCAN$member
      Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; prob.list[[l]]=prob;L.mat.list[[l]] = L.mat;member.list[[l]]=member
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 ="),c(round(lam1,2),round(lam2,2),round(lam3,2)),":"),paste("K =",as.numeric(dim(Theta_hat)[3]))))
      }
    }
    aBIC[is.na(aBIC)] = 10^10
    aBIC[aBIC==0] = abs(aBIC[aBIC!=0][1])*10
    n_lam3 = which(aBIC[1:L3] == min(aBIC[1:L3]))[1];lam3 = lambda3[n_lam3]

    # search lam2
    for (l2 in 1:L2) {
      lam2 = lambda2[l2];l = L3+l2
      SCAN = SCAN.hao(data, K, lam1, lam2, lam3, initial = F, initialize=out.initial, concave=concave)
      mu_hat <- SCAN$mu; Theta_hat <- SCAN$Omega; prob=SCAN$prob; L.mat=SCAN$L.mat; member = SCAN$member
      aBIC[l] = SCAN$bic
      Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; prob.list[[l]]=prob;L.mat.list[[l]] = L.mat; member.list[[l]]=member
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 ="),c(round(lam1,2),round(lam2,2),round(lam3,2)),":"),paste("K =",as.numeric(dim(Theta_hat)[3]))))
      }
    }
    aBIC[is.na(aBIC)] = 10^10
    aBIC[aBIC==0] = abs(aBIC[aBIC!=0][1])*10
    n_lam2 = which(aBIC[(L3+1):(L3+L2)] == min(aBIC[(L3+1):(L3+L2)]))[1];lam2 = lambda2[n_lam2]

    # search lam1
    for (l1 in 1:L1) {
      lam1 = lambda1[l1];l = L3+L2+l1
      SCAN = SCAN.hao(data, K, lam1, lam2, lam3, initial = F, initialize=out.initial, concave=concave)
      mu_hat <- SCAN$mu; Theta_hat <- SCAN$Omega; prob=SCAN$prob; L.mat=SCAN$L.mat; member = SCAN$member
      aBIC[l] = SCAN$bic
      Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; prob.list[[l]]=prob;L.mat.list[[l]] = L.mat; member.list[[l]]=member
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 ="),c(round(lam1,2),round(lam2,2),round(lam3,2)),":"),paste("K =",as.numeric(dim(Theta_hat)[3]))))
      }
    }
    aBIC[is.na(aBIC)] = 10^10
    aBIC[aBIC==0] = abs(aBIC[aBIC!=0][1])*10
    n_lam1 = which(aBIC[(L3+L2+1):(L3+L2+L1)] == min(aBIC[(L3+L2+1):(L3+L2+L1)]))[1];lam1 = lambda1[n_lam1]
  }
  n_lam = which(aBIC == min(aBIC))[1]
  Opt_aBIC = min(aBIC)

  opt_num = n_lam
  opt_Theta_hat = Theta_hat.list[[opt_num]]
  opt_Mu_hat = Mu_hat.list[[opt_num]]
  opt_prob = prob.list[[opt_num]]
  opt_L.mat = L.mat.list[[opt_num]]
  opt_member = member.list[[opt_num]]
  K_hat = dim(opt_Theta_hat)[3]

  time = as.numeric((proc.time() - t0)[3])
  
  result = list(Mu_hat.list=Mu_hat.list,
                Theta_hat.list=Theta_hat.list,prob.list=prob.list,L.mat.list=L.mat.list,
                Opt_aBIC=Opt_aBIC,BIC=aBIC,Opt_num=n_lam,member.list=member.list,
                opt_Theta_hat=opt_Theta_hat, opt_Mu_hat=opt_Mu_hat, opt_L.mat=opt_L.mat,
                opt_prob=opt_prob, opt_member=opt_member, K_hat=K_hat, time=time)
  return(result)
}

f.density = function(x.vec, mu1, Omega1){
  p = length(mu1)
  fdensity = as.numeric((((2*pi)^(-p/2)))*(det(Omega1))^(1/2) * exp((-1/2*t(x.vec-mu1)%*% Omega1 %*%(x.vec-mu1))))
  return(fdensity)
}

initialize_fuc = function(data, K, n.start = 100){

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: initialize_fuc
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Generating the initial values using K-means.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------

  n <- dim(data)[1]
  p <- dim(data)[2]
  Mu <- matrix(0, K, p)
  kmeans.clust <- kmeans(data, K, nstart = n.start)
  memb <- kmeans.clust$cluster
  prob <- kmeans.clust$size/n
  Theta <- array(0, dim = c(p, p, K))
  S <- array(0, dim = c(p, p, K))
  for(k in 1:K)
  {
    Mu[k,] <- t(colMeans(data[memb == k, , drop = FALSE]) )
    S[,,k]  <- cov(data[memb == k, , drop = FALSE])
    Theta[,,k] <- solve(S[,,k])
  }

  int.res <- list()
  int.res$prob <-  prob
  int.res$Mu <-  Mu
  int.res$Theta <- Theta
  int.res$S <- S
  int.res$memb <- memb
  return(int.res)
}

Symmetrize = function(X){

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: Symmetrize
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Symmetrize the precision matrices using the symmetrization strategy of Cai et al. (2016).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ X: p * p matrix.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ X: a symmetrical matrix.
  ## -----------------------------------------------------------------------------------------------------------------

  p = dim(X)[1]
  for(i in 1:p){
    for(j in i:p){
      if(X[i,j] < X[j, i]){
        X[j, i] = X[i, j]
      }else{
        X[i, j] = X[j, i]
      }
    }
  }
  return(X)
}

BIC = function(data, mu_hat, Theta_hat, L.mat){

  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: BIC
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Calculating the adaptive BIC-type criterion.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: f.den.vec()
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ mu_hat: K0_hat * p matrix, the estimated mean vectors of K0_hat subgroups.
  ## @ Theta_hat: p * p * K0_hat array, the estimated precision matrices of K0_hat subgroups.
  ## @ L.mat: n * K0_hat matrix, the estimated probability that each sample belongs to each subgroup.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list P including:
  ## @ fit.error: a float value, the value of the loss function (without penalty function).
  ## @ df: a float value, the penalty value for non-zero parameters corresponding the choice of given tuning parameters.
  ## @ bic: a float value, the BIC value corresponding the choice of given tuning parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------

  n = nrow(data)
  K = nrow(mu_hat)

  # fitting error
  pi_vec = apply(L.mat, 2, sum)/n
  fit.error_mat = matrix(0, n, K)
  for(k in 1:K) {
    fit.error_mat[,k] = pi_vec[k] * f.den.vec( data, as.numeric(mu_hat[k,]), Theta_hat[,,k] )
  }
  fit0 = apply(fit.error_mat, 1, sum)
  fit.error = sum(log( fit0 + min(fit0[fit0>0]) ))
  fit.error = - 2*fit.error

  # degrees of freedom
  for(i in 1:K){
    Theta_hat[upper.tri(Theta_hat[, , i], diag = T)] = 0
  }

  df =  log(n) * length(which(mu_hat != 0)) + 2 * length(which(Theta_hat != 0))
  bic = fit.error + df
  P = list()
  P$fit.error = fit.error
  P$df = df
  P$bic = bic
  return(P)
}

# -------------- The K-means + JGL method ----------------
kmeans.JGL = function(data, lambda2, lambda3, cluster.data, mu.hat, prob, member)
{

  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: kmeans.JGL
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Implementing the K-means + JGL method
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R packages: JGL
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ lambda2: a float value, the tuning parameter controlling the sparse of the precision matrix.
  ## @ lambda3: a float value, the tuning parameter controlling the similarity of subgroups.
  ## @ cluster.data: a list including data of all subgroups.
  ## @ mu.hat: The mean vectors estimated using K-means or hierarchical clustering.
  ## @ prob/member: Similar to the previous definition.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including:
  ## @ mu: K_hat * p matrix, the estimated mean vectors of K_hat subgroups.
  ## @ Theta: p * p * K_hat array, the estimated precision matrices of K_hat subgroups.
  ## @ prob: K_hat * 1 vector, the estimated mixture probabilities of subgroups.
  ## @ bic: a float value, the BIC value corresponding the choice of given tuning parameters.
  ## @ member: subgroup labels to which each sample belongs.
  ## ------------------------------------------------------------------------------------------------------------------------------------------

  p = dim(data)[2]
  K = length(unique(member))
  Theta = array(0, dim = c(p, p, K))
  out.jgl <- JGL(Y=cluster.data, penalty="group",lambda2,lambda3,return.whole.theta=TRUE)
  for(k in 1:K) {
    Theta[,,k] = out.jgl$theta[[k]]
  }
  bic = BIC.JGL(data, mu.hat, Theta, prob)
  JGL_res <- list();JGL_res$mu <- mu.hat;JGL_res$Theta <- Theta
  JGL_res$prob <- prob;JGL_res$bic <- bic$bic; JGL_res$member <- member
  return(JGL_res)
}

tuning.lambda.JGL = function(lambda, data, K_min=2,K_max=6,k.index= "gap", dbscan = F, trace = F, pca = F, pca.cut = 0.95)
{

  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: tuning.lambda.JGL
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Searching and selecting the optional tuning parameters under the adaptive BIC-type criterion
  ##            using the K-means + JGL method.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: NbClust()
  ##            R packages:  JGL
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ lambda: a list, the sequences of the tuning parameters (lambda1, lambda2, and lambda3).
  ## @ data: n * p matrix, the design matrix.
  ## @ K_min: int, a low bound of K_0.
  ## @ K_max: int, a upper bound of K_0.
  ## @ k.index: the selection criterion for K_hat using K-means, which can refer to R package NbClust.
  ## @ dbscan: a logical variable, whether to use hierarchical clustering.
  ## @ trace: the logical variable, whether or not to output the number of identified subgroups during the search for parameters.
  ## @ pca: the logical variable, whether to construct clustering to obtain the subgroups based on top principle components of data.
  ## @ pca.cut: When using pca, the selected cut-off value of cumulative proportion of variance.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list "result" including:
  ## @ Opt_lambda: the selected optional tuning parameters.
  ## @ Mu_hat.list: the estimated mean vectors of K0_hat subgroups corresponding all choices of given tuning parameters.
  ## @ Theta_hat.list: the estimated precision matrices of K0_hat subgroups corresponding all choices of given tuning parameters.
  ## @ prob.list: the estimated mixture probabilities of subgroups corresponding all choices of given tuning parameters.
  ## @ member.list: subgroup labels to which each sample belongs corresponding all choices of given tuning parameters.
  ## @ L.mat.list: the estimated probability that each sample belongs to each subgroup corresponding all choices of given tuning parameters.
  ## @ Opt_aBIC: the optional BIC value.
  ## @ Opt_num: the position of the optimal parameter for all given parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------

  t0 = proc.time()
  
  lambda1 = lambda$lambda1
  lambda2 = lambda$lambda2
  lambda3 = lambda$lambda3
  L1 = 0
  L2 = length(lambda2)
  L3 = length(lambda3)
  L = L1+L2+L3

  aBIC = rep(0,L)
  n_all = dim(data)[1]
  p = dim(data)[2]

  if(pca){
    cov = cov(data)
    pca = princomp(cov)
    pca = summary(pca)
    pca.var = as.numeric(pca$sdev)^2
    pca.cp = rep(0,p)
    for (jp in 1:p) {
      pca.cp[jp] = sum(as.numeric(pca.var/sum(pca.var))[1:jp])
    }
    r_p = which(pca.cp > pca.cut)[1]
    temp = eigen(cov)
    val = temp$values
    vec = temp$vectors
    pc = data %*% vec
    rp_X = as.matrix(pc[,1:r_p])

    if(!dbscan){
      nb_clust = NbClust(rp_X,distance = "euclidean",min.nc=K_min, max.nc=K_max, method = "kmeans",index = k.index,alphaBeale=0.1)
      member = nb_clust$Best.partition
    } else {
      out.initial = initialize_fuc.dbscan(rp_X,K)
      member = out.initial$memb
    }

  } else {
    if(!dbscan){
      nb_clust = NbClust(data,distance = "euclidean",min.nc=K_min, max.nc=K_max, method = "kmeans",index = k.index,alphaBeale=0.1)
      member = nb_clust$Best.partition
    } else {
      out.initial = initialize_fuc.dbscan(data,K)
      member = out.initial$memb
    }
  }





  K = length(unique(member))
  Theta = array(0, dim = c(p, p, K))
  mu.hat = matrix(0, K, p)
  prob = rep(0,K)
  cluster.data = vector("list", length = K)
  for(k in 1:K) {
    cluster.data[[k]] = data[member==k,]
    mu.hat[k,] = apply(cluster.data[[k]],2,mean)
    prob[k] = sum(member == k)/length(member)
  }

  Mu_hat.list = list()
  Theta_hat.list = list()
  prob.list = list()
  L.mat.list = list()
  member.list = list()

  # search lam3
  lam1 = median(lambda1);lam2 = median(lambda2)
  for (l in 1:L3) {
    lam3 = lambda3[l]
    SCAN = kmeans.JGL(data, lambda2=lam2, lambda3=lam3, cluster.data, mu.hat, prob, member)
    mu_hat <- SCAN$mu; Theta_hat <- SCAN$Theta; prob=SCAN$prob
    aBIC[l] = SCAN$bic; member <- SCAN$member
    Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat;
    prob.list[[l]]=prob; member.list[[l]] <- member
    Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; prob.list[[l]]=prob; member.list[[l]]=member
    if(trace){
      print(c(cat(paste(l,"lam1 lam2 lam3 ="),c(round(lam1,2),round(lam2,2),round(lam3,2)),":"),paste("K =",as.numeric(dim(Theta_hat)[3]))))
    }
  }
  aBIC[is.na(aBIC)] = 10^10
  aBIC[aBIC==0] = abs(aBIC[aBIC!=0][1])*10
  n_lam3 = which(aBIC[1:L3] == min(aBIC[1:L3]))[1];lam3 = lambda3[n_lam3]

  # search lam2
  for (l2 in 1:L2) {
    lam2 = lambda2[l2];l = L3+l2
    SCAN = kmeans.JGL(data, lambda2=lam2, lambda3=lam3, cluster.data, mu.hat, prob, member)
    mu_hat <- SCAN$mu; Theta_hat <- SCAN$Theta; prob=SCAN$prob
    aBIC[l] = SCAN$bic; member <- SCAN$member
    Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat;
    prob.list[[l]]=prob; member.list[[l]] <- member
    if(trace){
      print(c(cat(paste(l,"lam1 lam2 lam3 ="),c(round(lam1,2),round(lam2,2),round(lam3,2)),":"),paste("K =",as.numeric(dim(Theta_hat)[3]))))
    }
  }
  aBIC[is.na(aBIC)] = 10^10
  n_lam2 = which(aBIC[(L3+1):(L3+L2)] == min(aBIC[(L3+1):(L3+L2)]))[1];lam2 = lambda2[n_lam2]

  n_lam = which(aBIC == min(aBIC))[1]
  Opt_aBIC = min(aBIC)

  opt_num = n_lam
  opt_Theta_hat = Theta_hat.list[[opt_num]]
  opt_Mu_hat = Mu_hat.list[[opt_num]]
  opt_prob = prob.list[[opt_num]]
  opt_member = member.list[[opt_num]]
  K_hat = dim(opt_Theta_hat)[3]

  time = as.numeric((proc.time() - t0)[3])
  
  result = list(Mu_hat.list=Mu_hat.list, Theta_hat.list=Theta_hat.list, prob.list=prob.list,
                Opt_aBIC=Opt_aBIC, BIC=aBIC, opt_num=n_lam, member.list=member.list,
                opt_Theta_hat=opt_Theta_hat, opt_Mu_hat=opt_Mu_hat,
                opt_prob=opt_prob, opt_member=opt_member, K_hat=K_hat, time=time)
  return(result)
}

BIC.JGL = function(data, mu_hat, Theta_hat, prob){
  n = nrow(data)
  K = nrow(mu_hat)

  # fitting error
  pi_vec = prob
  fit.error_mat = matrix(0, n, K)
  fit.error = 0
  for(k in 1:K) {
    fit.error_mat[,k] = pi_vec[k] * f.den.vec( data, as.numeric(mu_hat[k,]), Theta_hat[,,k] )
  }
  fit0 = apply(fit.error_mat, 1, sum)
  fit.error = sum(log( fit0 + min(fit0[fit0>0]) ))
  fit.error = - 2*fit.error

  for(i in 1:K){
    Theta_hat[upper.tri(Theta_hat[, , i], diag = T)] = 0
  }

  df =  log(n) * length(which(mu_hat != 0)) + 2 * length(which(Theta_hat != 0))
  bic = fit.error + df
  P = list()
  P$fit.error = fit.error
  P$df = df
  P$bic = bic
  return(P)
}

f.den.vec = function(data, mu, Theta){

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: f.den.vec
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            calculate the density function values at each sample point.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ mu1: p * 1 vector, the mean vector.
  ## @ Omega1: p * p matrix, the precision matrix.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ fdensity: The density function values at each sample point.
  ## -----------------------------------------------------------------------------------------------------------------

  p = length(mu)
  fden = as.numeric( (2*pi)^(-p/2) * (det(Theta))^(1/2) * exp(-1/2*diag(t(t(data) - as.numeric(mu)) %*% Theta %*% (t(data) - as.numeric(mu)))) )
  return(fden)
}





# -------------- The mtlgmm method (Arxiv, 2023)----------------
tlgmm.summary = function(t.data, A.data, lambda_choice="fixed", initial_method = "kmeans"){
  t0 = proc.time()
  fit_mtl <- mtlgmm(x = A.data, initial_method = initial_method,
                    lambda_choice = lambda_choice, step_size = "lipschitz")
  time_mtl = as.numeric((proc.time() - t0)[3])
  fit_tl <- tlgmm(x = t.data, fitted_bar = fit_mtl, kappa0 = 1/3,
                  initial_method = initial_method, lambda_choice = lambda_choice, step_size = "lipschitz")
  time_tl = as.numeric((proc.time() - t0)[3])
  
  K_hat = 2
  n = dim(t.data)[1]
  p = dim(t.data)[2]
  
  Mu_hat = matrix(0, K_hat, p)
  Mu_hat[1, ] = fit_tl$mu1
  Mu_hat[2, ] = fit_tl$mu2
  
  Theta_hat = array(0, list(p,p,K_hat))
  Theta_hat[,,1] = solve(fit_tl$Sigma)
  Theta_hat[,,2] = Theta_hat[,,1]
  
  member = predict_gmm(w = fit_tl$w, mu1 = fit_tl$mu1, mu2 = fit_tl$mu2, beta = fit_tl$beta, newx = t.data)
  
  res = list(Mu_hat=Mu_hat, Theta_hat=Theta_hat, member=member,
             time_mtl=time_mtl, time_tl=time_tl)
}


