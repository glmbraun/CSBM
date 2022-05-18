import numpy as np


def convertClust(clust):
    n = len(clust)
    k= len(np.unique(clust))
    Z = np.zeros((n,k))
    for i in range(n):
        Z[i,clust[i]] = 1
    return np.array(Z,dtype=int)

def convertClust2(clust,k):
    n = len(clust)
    Z = np.zeros((n,k))
    W = np.zeros((n,k))
    for i in range(n):
        Z[i,clust[i]] = 1
    for kk in range(k):
        if Z.sum(0)[kk]>0:
            W[:,kk]= Z[:,kk]/Z.sum(0)[kk] 
    return [np.array(Z,dtype=int),np.array(W)]

def convertZ(Z):
    n = Z.shape[0]
    clust = [ np.where(Z[i,]==1)[0][0] for i in range(n) ]
    return np.array(clust,dtype=int)

# Define clustering methods
def iter_sbm(A,Z_init,k, iter=10):
    n = Z_init.shape[0]
    W = Z_init/Z_init.sum(0)
    clust = convertZ(Z_init)
    for rep in range(iter):
        C = A @ W
        Pi = np.transpose(W)@A@W
        for i in range(n):
            ls = [np.linalg.norm(C[i,]-Pi[kk,]) for kk in range(k)]
            clust[i] = np.argmin(ls)
        W = convertClust(clust)/convertClust(clust).sum(0)
        return clust

#IR-sLS   
def iter_csbm(A,X,Z_init,k, iter=10):
    n = Z_init.shape[0]
    W = Z_init/Z_init.sum(0)
    clust = convertZ(Z_init)
    for rep in range(iter):
        Z = convertClust(clust)
        C = A @ W
        Pi = np.transpose(W)@A@W
        mu = np.transpose(W)@X
        sigma = [np.array([np.linalg.norm(X[i,]-mu[kk,])**2 for i in range(306) if Z[i,kk]==1]).mean() for kk in range(k)]
        lb = Z.sum(0).min()/(Pi.max()) 
        for i in range(n):
            ls = [lb*np.linalg.norm((C[i,]-Pi[kk,]))**2 for kk in range(k)]
            ls2 = [np.linalg.norm(X[i,]-mu[kk,])**2/sigma[kk] for kk in range(k)]
            clust[i] = np.argmin(np.add(ls,ls2))
        W = convertClust(clust)/convertClust(clust).sum(0)
    return clust

#IR-LS
def iter_csbm2(A,X,Z_init,k, iter=10):
    n = Z_init.shape[0]
    W = Z_init/Z_init.sum(0)
    clust = convertZ(Z_init)
    for rep in range(iter):
        Z = convertClust(clust)
        C = A @ W
        Pi = np.transpose(W)@A@W
        mu = np.transpose(W)@X
        sigma = [np.array([np.linalg.norm(X[i,]-mu[kk,])**2 for i in range(306) if Z[i,kk]==1]).mean() for kk in range(k)]
        sbm_var = [np.diag([np.sqrt(sum(Z[:,kk])/np.absolute(Pi[com,kk])) for kk in range(k)]) for com in range(k)]
        lb = Z.sum(0).min()/(Pi.max()) 
        for i in range(n):
            ls = [np.linalg.norm((C[i,]-Pi[kk,])@sbm_var[kk])**2 for kk in range(k)]
            ls2 = [np.linalg.norm(X[i,]-mu[kk,])**2/sigma[kk] for kk in range(k)]
            clust[i] = np.argmin(np.add(ls,ls2))
        W = convertClust(clust)/convertClust(clust).sum(0)
    return clust