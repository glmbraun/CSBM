# Experiences with random initialization

xp_test<-function(c){
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
  Ker = matrix(NA,nrow=n,ncol=1)
  for(i in 1:n){
    Ker[i,]=convZ[i]+sigma2*rnorm(1)
  }
  
  data_mix<-cbind(svds(A,2)$v, Ker)
  gmm<-GMM(data_mix,2)
  pr = predict_GMM(data_mix, gmm$centroids, gmm$covariance_matrices, gmm$weights)
  clust_gmm = pr$cluster_labels

  Z_init =pure_membership(n,2,rep(1/2,2))
  
  #Z_init=convertClust(clust_gmm)
  
  clust_map = iterCSBM_map(n,2,A,Ker,Z_init,15,sigma2)
  clust_simp = iterCSBM_simp(n,2,A,Ker,Z_init,15,sigma2)
  # 
   score2 = matrix(NA, nrow = 3, ncol = 1)
  score2[1,1] = NMI(clust_map,convZ)
  score2[2,1] = NMI(clust_simp,convZ)
  score2[3,1] = NMI(Z_init%>%convertZ(),convZ)
  return(score2)
}


library(matrixStats)

seq_c=seq(from=0.1,to=0.5,by=0.05)
score=matrix(NA,nrow=3,ncol=length(seq_c)+1)
score_var=matrix(NA,nrow=3,ncol=length(seq_c)+1)
for(i in 1:length(seq_c)){
  print(seq_c[i])
  test= sapply(1:20,function(t){xp_test(seq_c[i])})
  score[,i]=test%>%as.matrix()%>%rowMeans()
  score_var[,i]=test%>%as.matrix()%>%rowVars()
}

write.csv(score,"score_rd.csv",row.names = F)
write.csv(score_var,"score_var.csv",row.names = F)
