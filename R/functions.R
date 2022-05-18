###############################
##### Auxiliary functions #####
###############################

norm_vec2 <- function(x) sum(x^2)

#Convert a membership matrix to a vector 
convertZ<-function(Z){
  n<-dim(Z)[1]
  return(sapply(1:n, function(i) {match(1, Z[i,])}))
}

# Convert a vector of labels to a membership matrix
convertClust<-function(clust){
  n<-length(clust)
  k<-length(unique(clust))
  Z<-matrix(0,nrow=n,ncol=k)
  for(i in 1:n){Z[i,clust[i]]<-1}
  return(Z)
}


hard_threshold = function(M, delta){
  # delta : threshold value
  M[(M) <= delta] = 0
  return(M)
}


normalize_row<-function(M){
  return(diag(1/rowSums(abs(M)))%*%M)
}


laplacian <- function(A, tau = 0) {
  Atau = A+tau
  degree = apply(Atau,1,sum)
  diag(1/sqrt(degree))%*%Atau%*%diag(1/sqrt(degree))
}

laplacian_reg <- function(A, tau = 0) {
  degree = apply(A,1,sum) + tau
  diag(1/sqrt(degree))%*%A%*%diag(1/sqrt(degree))
}


clust_on_mat_eig<-function(A,k){
  return(kmeans(svds(A,k)$u,k)$cluster)
}

spectral.norm <- function( x )
{
  if ( !is.numeric( x ) ) {
    stop( "argument x is not numeric" )
  }
  if ( is.vector( x ) ) {
    return( sqrt( sum( x * x ) ) )
  }
  if ( !is.matrix( x ) ) {
    return( "argument x is not a matrix" )
  }
  A <- t(x) %*% x
  eigenA <- eigen( A )
  lambdaA <- eigenA$values
  maxLambdaA <- lambdaA[1]
  if ( maxLambdaA < 0 ) {
    stop( "t(x) %*% x is negative definite" )
  }
  return( sqrt( maxLambdaA ) )
}

# MAP
map_sbm <- function(A,k,Pi, z, i){
  j_values = 1:n
  j_values= j_values[j_values!= i]
  f<-function(j){A[i,j]*log(Pi[k, z[j]])+(1-A[i,j])*log(1-Pi[k, z[j]])}
  return(sapply(j_values,f)%>%sum())
}

map_csbm <- function(A,X, k, Pi, mu, z, i,sigma ){
  #sigma = 0.2 #to simplify in a first step
  j_values = 1:n
  j_values= j_values[j_values!= i]
  f<-function(j){A[i,j]*log(Pi[k, z[j]])+(1-A[i,j])*log(1-Pi[k, z[j]])}
  
  res =sapply(j_values,f)%>%sum()
  return( res - norm(X[i,]-mu[k,],"2")^2/sigma )
}

############################
##### Graph generation #####
############################

# Generate an inhomogeneous Erdös-Renyi graph from a matrix P
sample_IER<-function(P){
  A = apply(P,MARGIN = c(1,2),function(u) rbinom(1,1,prob = min(u,1)))  
  A[lower.tri(A)] = 0
  A =  A + t(A)
  diag(A) = 0
  return(A)
}

# Generate a membership matrix with no overlap
pure_membership<-function(n, K, alpha){
  return(t(rmultinom(n,1,alpha)))
}

###########################
##### Main algorithms #####
###########################

# IR-sLS with knwon sigma
iterCSBM_simp<-function(n,K,A,X, Z_init, iter, sigma){
  cZ = convertZ(Z_init)
  Z = Z_init
  Zn= Z%*% diag(1/colSums(Z))
  Pi = t(Zn)%*% A %*% Zn
  mu = t(Zn)%*% X
  for (t in 1:iter){
    # Estimate n_min and p_max
    n_min = min(colSums(Z))
    p_max = max(Pi)
    lambda = n_min/p_max
    # Refinement of the partition
    for(i in 1:n){
      f<-function(k){ return(lambda*norm_vec2((A[i,]%*%Zn-Pi[k,]))+norm_vec2((X[i,]-mu[k,]))/sigma) } # we assume sigma^2 is known 
      cZ[i] = which.min(sapply(1:K,f))
    }
    # Update model parameters
    Z = convertClust(cZ)
    Zn = Z%*% diag(1/colSums(Z))
    Pi = t(Zn)%*% A %*% Zn
    mu = t(Zn)%*% X
  }
  return(cZ)
}

# IR-LS with known sigma
iterCSBM_snr<-function(n,K,A,X, Z_init, iter,sigma){
  cZ = convertZ(Z_init)
  Z = Z_init
  Zn= Z%*% diag(1/colSums(Z))
  Pi = t(Zn)%*% A %*% Zn
  mu = t(Zn)%*% X
  for (t in 1:iter){
    # Compute intra-class (sqrt) variance for each community 
    intra_varA<-function(k) {
      return(diag(sqrt(Pi[k,]*(1-Pi[k,])/colSums(Z))^(-1)))
    }
    Sigma1 = lapply(1:K,function(k) intra_varA(k)) 
    # Refinement of the partition
    for(i in 1:n){
      f<-function(k){ return(norm_vec2((A[i,]%*%Zn-Pi[k,])%*%Sigma1[[k]])+norm_vec2((X[i,]-mu[k,]))/sigma ) } 
      cZ[i] = which.min(sapply(1:K,f))
    }
    # Update model parameters
    Z = convertClust(cZ)
    Zn = Z%*% diag(1/colSums(Z))
    Pi = t(Zn)%*% A %*% Zn
    mu = t(Zn)%*% X
  }
  return(cZ)
}


#IR-LS with unknown sigma under spherical assumption for the covariates
iterCSBM<-function(n,K,A,X, Z_init, iter){
  cZ = convertZ(Z_init)
  Z_init = Z_init%*% diag(1/colSums(Z_init))
  Pi = t(Z_init)%*% A %*% Z_init
  mu = t(Z_init)%*% X
  Z = Z_init
  for (t in 1:iter){
    # Compute intra-class (sqrt) variance for each community 
    intra_varA<-function(k) {
      return(diag(sqrt(Pi[k,]*(1-Pi[k,]))))
    }
    intra_varX<- function(k){
      subX = X[convertClust(cZ)[,k]%>%as.logical(),]
      return(var(subX))
    }
    Sigma1 = lapply(1:K,function(k) inv(intra_varA(k)+diag(10^-6,K)) ) #numerical stability issues
    Sigma2 = lapply(1:K, function(k) intra_varX(k)%>%inv()%>%sqrtm())
    # Refinement of the partition
    for(i in 1:n){
      f<-function(k){ return(norm((A[i,]%*%Z-Pi[k,])%*%Sigma1[[k]],"2")^2+norm((X[i,]-mu[k,])%*%Sigma2[[k]],"2")^2) }
      cZ[i] = which.min(sapply(1:K,f))
    }
    # Update model parameters
    Z = convertClust(cZ)
    Z = Z%*% diag(1/colSums(Z))
    Pi = t(Z)%*% A %*% Z
    mu = t(Z)%*% X
  }
  return(cZ)
}

#IR-MAP for CSBM
iterCSBM_map<-function(n,K,A,X, Z_init, iter,sigma){
  cZ = convertZ(Z_init)
  Z_init = Z_init%*% diag(1/colSums(Z_init))
  Pi = t(Z_init)%*% A %*% Z_init
  mu = t(Z_init)%*% X
  Z = Z_init
  for (t in 1:iter){
    # Refinement of the partition
    for(i in 1:n){
      score_vect = sapply(1:K,function(k) map_csbm(A,X,k,Pi,mu,cZ,i,sigma))
      cZ[i] = which.max(score_vect)
    }
    # Update model parameters
    Z = convertClust(cZ)
    Z = Z%*% diag(1/colSums(Z))
    Pi = t(Z)%*% A %*% Z
    mu = t(Z)%*% X
  }
  return(cZ)
}


# IR-LS for SBM without covariates
iterSBM<-function(n,K,A, Z_init, iter){
  cZ = convertZ(Z_init)
  Z_init = Z_init%*% diag(1/colSums(Z_init))
  Pi = t(Z_init)%*% A %*% Z_init
  Z = Z_init
  for (t in 1:iter){
    # Compute intra-class (sqrt) variance for each community 
    intra_varA<-function(k) {
      return(diag(sqrt(Pi[k,]*(1-Pi[k,]))))
    }
    Sigma1 = lapply(1:K,function(k) inv(intra_varA(k)+diag(10^-6,K)) )
    # Refinement of the partition
    for(i in 1:n){
      f<-function(k){ return(norm((A[i,]%*%Z-Pi[k,])%*%Sigma1[[k]],"2")^2) }
      cZ[i] = which.min(sapply(1:K,f))
    }
    # Update model parameters
    Z = convertClust(cZ)
    Z = Z%*% diag(1/colSums(Z))
    Pi = t(Z)%*% A %*% Z
  }
  return(cZ)
}

#IR-sLS for SBM without covariates
iterSBM_mod<-function(n,K,A, Z_init, iter){
  cZ = convertZ(Z_init)
  Z_init = Z_init%*% diag(1/colSums(Z_init))
  Pi = t(Z_init)%*% A %*% Z_init
  Z = Z_init
  for (t in 1:iter){
    # Refinement of the partition
    for(i in 1:n){
      f<-function(k){ return(norm((A[i,]%*%Z-Pi[k,]),"2")^2) }
      cZ[i] = which.min(sapply(1:K,f))
    }
    # Update model parameters
    Z = convertClust(cZ)
    Z = Z%*% diag(1/colSums(Z))
    Pi = t(Z)%*% A %*% Z
  }
  return(cZ)
}

#IR-MAP for SBM
iterSBM_map<-function(n,K,A, Z_init, iter){
  # First estimation of the model parameters
  cZ = convertZ(Z_init)
  Z_init = Z_init%*% diag(1/colSums(Z_init))
  Pi = t(Z_init)%*% A %*% Z_init
  Z = Z_init
  for (t in 1:iter){
    # Refinement of the partition
    for(i in 1:n){
      score_vect = sapply(1:K,function(k) map_sbm(A,k,Pi,cZ,i))
      cZ[i] = which.max(score_vect)
    }
    # Update model parameters
    Z = convertClust(cZ)
    Z = Z%*% diag(1/colSums(Z))
    Pi = t(Z)%*% A %*% Z
  }
  return(cZ)
}


# commuLloyd from "Statistical and Computational Guarantees of Lloyd’s Algorithm and Its Variants"
commuLloyd<-function(n,K,A,Z_init,iter) {
  cZ=convertZ(Z_init)
  W = Z_init%*% diag(1/colSums(Z_init))
  B = A%*%W
  for(t in 1:iter){
    for(i in 1:n){
      cZ[i]=which.max(B[i,])
    }
    Z = convertClust(cZ)
    W = Z%*%diag(1/colSums(Z))
    B = A%*%W
  }
  return(cZ)
}

#################################
##### Grid search functions #####
#################################

grid_search_oracle <- function( A, Ker, cZ, K, stop=500, step){
  lambda = seq(0, step*stop, by = step)
  objective_values = c()
  for(i in 1:stop){
    L = A + lambda[i]*Ker
    clust = clust_on_mat_eig(L,K)
    objective_values =c(objective_values,NMI(clust,cZ))
  }
  return(lambda[which.max(objective_values)])
}

grid_search_ker<- function( A, Ker, K, stop=500, step, start = 0){
  lambda = seq(start, start+step*stop, by = step)
  objective_values = c()
  for(i in 1:stop){
    L= A + lambda[i]*Ker
    singular_values =sort(svds(L,K+1)$d,decreasing = T)
    objective_values =c(objective_values, 1-singular_values[K+1]/singular_values[K])
  }
  return(lambda[which.max(objective_values)])
}