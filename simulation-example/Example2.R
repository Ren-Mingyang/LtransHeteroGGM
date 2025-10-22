################################################################################################
#Codes for conducting simulation studies
################################################################################################
rm(list = ls(all = TRUE))
ls()
###############################

library(HeteroGGM)
library(MASS)
library(NbClust)
library(JGL)
library(Matrix)
library(mtlgmm)
setwd("D:/SJTU/2mypaper/transfer_GGMM/2numerical_study/simulation-example")
R.files = list.files("../main_functions")
L=length(R.files)
for (l in 1:L) {
  source(paste0("../main_functions/",R.files[l]))
}
source("sim-func.R")


################ Simulation parameters ################
n = 200                    # The sample size of each subgroup
p = 100                    # The dimension of the precision matrix
M0 = 3                     # The true number of subgroups
n.vec = rep(n,M0)          # The sample sizes of M0 subgroups
M = 2*M0                   # The given upper bound of M0
K = 6                      # The number of auxiliary domains
K.A = 3                    # The number of auxiliary domains
M0.A.vec=c(2,3,4,3,3,3)    # The true number of subgroups in each auxiliary domains
N0 = 3*n                   # The sample size of each subgroup in auxiliary domains
mue0 = 2                   # The signal strength of non-zero mean elements

################ Generating true parameters ################
para.list.target = gen.target.para(p, M0, type="power.law", mue=mue0, seed=2025)
para.list.aux = gen.aux.para(n, p, para.list.target, K, K.A, M0.A.vec, mue=mue0)

################ Generating simulated data ################
target.whole.data = generate.data(n.vec, para.list.target)
aux.whole.data.all = generate.aux.data(N0, para.list.aux)
t.data = target.whole.data$data
A.data = aux.whole.data.all$A.data


t0 = proc.time()
# The proposed transfer GGMM
lambda.t = genelambda.obo(nlambda1=1,lambda1_max=0.25,lambda1_min=0.25,
                          nlambda2=1,lambda2_max=0.25,lambda2_min=0.25,
                          nlambda3=5,lambda3_max=5,lambda3_min=1)
lambda.t.refit = genelambda.obo(nlambda1=2,lambda1_max=0.5,lambda1_min=0.25,
                                nlambda2=10,lambda2_max=1,lambda2_min=0.05,
                                nlambda3=1,lambda3_max=1,lambda3_min=1)
lambda.t.refit$lambda3 = 0
lambda = genelambda.obo(nlambda1=2,lambda1_max=0.5,lambda1_min=0.25,
                        nlambda2=5,lambda2_max=1,lambda2_min=0.25,
                        nlambda3=5,lambda3_max=3,lambda3_min=1)
lambda$lambda3 = 0
lambda.A.list = list()
for (k in 1:K) {
  lambda.A.list[[k]] = lambda
}
M.A.vec = M0.A.vec
res.LocalTrans = trans_GGMM(t.data, lambda.t, M, A.data, lambda.A.list, M.A.vec,
                       cn.lam2=seq(1,0.1,length.out=10), trace=T, lambda.t.refit=lambda.t.refit)
# res.LocalTrans = trans_GGMM(t.data, lambda.t, M, A.data, lambda.A.list, M.A.vec,
#                        cn.lam2=seq(0.1,1.0,length.out=10),
#                        cov.method="weight", preselect.aux=1, sel.type="L2",
#                        trace=T, lambda.t.refit=lambda.t.refit)
res.target = res.LocalTrans$res.target
res.target0 = res.LocalTrans$res.target0
index.pro = Evaluation.GGMM(t.data, res.target$opt_Mu_hat, res.target$opt_Theta_hat, target.whole.data$Mu0, target.whole.data$Theta0, M0, res.target$opt_L.mat, target.whole.data$L0, res.target$opt_prob)
index.single.FGGM = Evaluation.GGMM(t.data, res.target0$opt_Mu_hat, res.target0$opt_Theta_hat, target.whole.data$Mu0, target.whole.data$Theta0, M0, res.target0$opt_L.mat, target.whole.data$L0, res.target0$opt_prob)

# tlgmm
res.tlgmm = tlgmm.summary(t.data, A.data)
index.tlgmm = Evaluation.mtlgmm(t.data, res.tlgmm$Mu_hat, res.tlgmm$Theta_hat, target.whole.data$Mu0, target.whole.data$Theta0, M0, res.tlgmm$member, target.whole.data$L0)

# SCAN
lambda.SCAN = genelambda.obo(nlambda1=5,lambda1_max=0.5,lambda1_min=0.1,
                             nlambda2=10,lambda2_max=2,lambda2_min=0.1,
                             nlambda3=10,lambda3_max=2,lambda3_min=0.1)
lambda.SCAN$lambda3 = 0
res.SCAN = tuning.lambda.SCAN(lambda.SCAN, t.data, M0)
index.SCAN = Evaluation.GGMM(t.data, res.SCAN$opt_Mu_hat, res.SCAN$opt_Theta_hat, target.whole.data$Mu0, target.whole.data$Theta0, M0, res.SCAN$opt_L.mat, target.whole.data$L0, res.SCAN$opt_prob)
res.SCAN.e = tuning.lambda.SCAN(lambda.SCAN, t.data, 2)
index.SCAN.e = Evaluation.GGMM(t.data, res.SCAN.e$opt_Mu_hat, res.SCAN.e$opt_Theta_hat, target.whole.data$Mu0, target.whole.data$Theta0, M0, res.SCAN.e$opt_L.mat, target.whole.data$L0, res.SCAN.e$opt_prob)

# K-means + JGL
lambda.JGL = genelambda.obo(nlambda2=20,lambda2_max=5,lambda2_min=0.01,
                             nlambda3=20,lambda3_max=5,lambda3_min=0.01)
res.JGL = tuning.lambda.JGL(lambda.JGL, t.data)
index.JGL = Esti.error.JGL(t.data, res.JGL$opt_Mu_hat, res.JGL$opt_Theta_hat, target.whole.data$Mu0, target.whole.data$Theta0, M0, target.whole.data$L0, res.JGL$opt_member)
proc.time() - t0



index.pro
index.tlgmm
index.single.FGGM
index.SCAN
index.SCAN.e
index.JGL

res.LocalTrans$time_trans
res.JGL$time
res.LocalTrans$time_init
res.SCAN$time
res.SCAN.e$time
res.JGL$time










