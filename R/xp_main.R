library(gtools)
library(dplyr)
library(RSpectra)
library(mvnfast)
library(Rcsdp)
library(MASS)
library(aricode)
library(matlib)
library(expm)
library(KRLS)
library(ggplot2)
library(igraph)
library(xtable)
library(reshape2)
library(ClusterR)

###################################################################################
# Unidentifiable community by considering each source of information individually
###################################################################################

# Parameters :
K = 3 
n = 1000
Z = pure_membership(n, K,rep(1/3,3))
convZ = convertZ(Z)
#Pi = 0.02*matrix(c(1.6,1.2,0.5,1.2,1.6,0.5,0.05,0.05,1.2),3,3)
Pi = 0.02*matrix(c(1.5,1.5,0.5,1.5,1.5,0.5,0.5,0.5,1.5),3,3)
P = Z%*% Pi %*% t(Z)

mu = rbind(c(0,0,1),c(-1,1,0),c(0,0,1)) # centroids of the GMM
sigma = diag(0.2,3) # covariance matrix
s = 0.2
# Generate random inputs (to be used for Matlab experiments with SDP-Comb)
rep =20
adj_mat_ls = lapply(1:rep,function(i) sample_IER(P))
covariates_ls = lapply(1:rep,function(i) t(sapply(convZ, function(j){return(mvrnorm(1,mu=mu[j,], Sigma =sigma))})))

# Compute the score for each clustering method
score2 = matrix(NA, nrow = 7, ncol = rep)
for (rep in 20:rep){
  L = laplacian_reg(adj_mat_ls[[rep]],0.01)
  clust_L = clust_on_mat_eig(L,3)
  
  Ker = gausskernel(X=covariates_ls[[rep]],sigma=sqrt(s))
  clust_K = clust_on_mat_eig(Ker,3)
  
  #lambda = grid_search_oracle(L, Ker, convZ, 3,stop=100,step=0.001)
  #clust_Lr = clust_on_mat_eig(L+ lambda*Ker,3)
  
  data_mix<-cbind(svds(adj_mat_ls[[rep]],3)$v, covariates_ls[[rep]])
  gmm<-GMM(data_mix,3)
  pr = predict_GMM(data_mix, gmm$centroids, gmm$covariance_matrices, gmm$weights)   
  clust_gmm = pr$cluster_labels
  
  Z_init =convertClust(clust_gmm)
  #Z_init= pure_membership(n, K,rep(1/3,3))
  #Z_init= convertClust(clust_L)
  
  #clust_map = iterCSBM_map(n,3,adj_mat_ls[[rep]],covariates_ls[[rep]],Z_init,15,s)
  clust_simp = iterCSBM_simp(n,3,adj_mat_ls[[rep]],covariates_ls[[rep]],Z_init,10,s)
  clust_ls = iterCSBM_snr(n,3,adj_mat_ls[[rep]],covariates_ls[[rep]],Z_init, 10,s)
  
  score2[1,rep] = NMI(clust_L,convZ)
  #score2[2,rep] = NMI(clust_Lr,convZ)
  score2[3,rep] = NMI(clust_gmm,convZ)
  score2[4,rep] = NMI(clust_ls,convZ)
  score2[5,rep] = NMI(clust_simp,convZ)
  #score2[6,rep] = NMI(clust_map,convZ)
  score2[7,rep] = NMI(clust_K,convZ)
}


# Plot
rownames(score2) =c("L-SC", "ORL-SC", "EM-Emb","IR-LS", "IR-sLS","IR-MAP","K-SC")
score2 = t(score2)
score2= score2%>%as.data.frame()
mScore2 = melt(score2)


p2 <- ggplot(mScore2, aes(factor(variable), value)) + 
  geom_boxplot()
p2

write.csv(score2,"score_csbm2_bis3.csv")
write.csv(score2,"score_csbm2.csv")

################################################################### 
# Comparison between iter_map and iter_csbm when the signal vary
###################################################################
 K=2
 mu = rbind(c(1,0),c(0,-1))
 n=1000
 p=4
 q=1
 Pi = (p-q)*diag(1,2)+q*matrix(1,2,2)
 c=1
 signal = c*log(n)*(sqrt(p)-sqrt(q))^2
 sigma= 1/(2*signal)
 

 xp3<-function(c){
   p=4*log(n)/n
   q=log(n)/n
   Pi = c*(p-q)*diag(1,2)+c*q*matrix(1,2,2)
   Z = pure_membership(n,2,rep(1/2,2))
   convZ = convertZ(Z)
   P = Z%*% Pi %*% t(Z)
   A = sample_IER(P)
   A=A%>%as.matrix()
   
   signal = c*log(n)
   sigma = 1/(2*signal)*diag(2)
   sigma2 = 1/(2*signal)
   Ker = sapply(convZ, function(j){return(mvrnorm(1,mu=mu[j,], Sigma =sigma))})%>%t()
   print(Ker[1,1])
   
   data_mix<-cbind(svds(A,2)$v, Ker)
   gmm<-GMM(data_mix,2)
   pr = predict_GMM(data_mix, gmm$centroids, gmm$covariance_matrices, gmm$weights)
   clust_gmm = pr$cluster_labels

   Z_init =convertClust(clust_gmm)
   # 
   # clust_map = iterCSBM_map(n,2,A,Ker,Z_init,15,sigma2)
   clust_simp = iterCSBM_simp(n,2,A,Ker,Z_init,15,sigma2)
   # 
   score2 = matrix(NA, nrow = 3, ncol = 1)
   # score2[1,1] = NMI(clust_map,convZ)
   score2[2,1] = NMI(clust_simp,convZ)
   score2[3,1] = NMI(clust_gmm,convZ)
   return(score2)
 }
 

 
 
 
 seq_c=seq(from=0.1,to=0.5,by=0.025)
 score=matrix(NA,nrow=3,ncol=4)
 for(c in seq_c){
   score[,]
 }
test= sapply(1:2,function(i){xp3(0.1)})
score[,1]=test%>%as.matrix()%>%rowMeans()

score2 = matrix(NA, nrow = 3, ncol = rep)
c=0.1
for(r in 1:rep){
  p=4*log(n)/n
  q=log(n)/n
  Pi = c*(p-q)*diag(1,2)+c*q*matrix(1,2,2)
  Z = pure_membership(n,2,rep(1/2,2))
  convZ = convertZ(Z)
  P = Z%*% Pi %*% t(Z)
  A = sample_IER(P)
  A=A%>%as.matrix()
  
  signal = c*log(n)
  sigma = 1/(2*signal)*diag(2)
  sigma2 = 1/(2*signal)
  Ker = sapply(convZ, function(j){return(mvrnorm(1,mu=mu[j,], Sigma =sigma))})%>%t()
  
  print(Ker[1,1])
  
  data_mix<-cbind(svds(A,2)$v, Ker)
  gmm<-GMM(data_mix,2)
  pr = predict_GMM(data_mix, gmm$centroids, gmm$covariance_matrices, gmm$weights)   
  clust_gmm = pr$cluster_labels
  
  Z_init =convertClust(clust_gmm)
  
  clust_map = iterCSBM_map(n,2,A,Ker,Z_init,15,sigma2)
  clust_simp = iterCSBM_simp(n,2,A,Ker,Z_init,15,sigma2)
  
  score2[1,r] = NMI(clust_map,convZ)
  score2[2,r] = NMI(clust_simp,convZ)
  score2[3,r] = NMI(clust_gmm,convZ)
  
}