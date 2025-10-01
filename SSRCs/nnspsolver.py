# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:02:43 2024

@author: fredy.vides
"""

def nnspsolver(A,Y,L=100,nz=100,tol=1e-2,delta=1e-2):
    from scipy.optimize import nnls
    from numpy.linalg import svd,norm
    from numpy import zeros,dot,argsort,inf
    N = A.shape[1]
    M = Y.shape[1]
    if nz<0:
        nz = N
    X=zeros((N,M))
    X0 = X.copy()
    u,s,_=svd(A,full_matrices=0)
    rk=sum(s>tol)
    u=u[:,:rk]
    A=dot(u.T,A)
    Y=dot(u.T,Y)
    X0[:,0] = nnls(A,Y[:,0])[0]
    for k in range(M):
        w=zeros((N,))
        K=1
        Error=1+tol
        c=X0[:,k]
        x0=c
        ac=abs(c)
        f=argsort(-ac)
        N0=int(min(sum(ac[f]>delta),nz))
        while (K<=L) & (Error>tol):
            ff=f[:N0]
            X[:,k]=w
            c = nnls(A[:,ff],Y[:,k])[0]
            X[ff,k]=c
            Error=norm(x0-X[:,k],inf)
            x0=X[:,k]
            ac=abs(x0)
            f=argsort(-ac)
            N0=int(min(sum(ac[f]>delta),nz))
            K=K+1
    return X