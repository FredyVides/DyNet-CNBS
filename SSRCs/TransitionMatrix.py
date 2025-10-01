# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:06:27 2024

@author: fredy.vides
"""

def TransitionMatrix(data,SM=1,sample_size=39,error = 1e-6,nz = 100,SB = [],SB_Index=[]):
    from numpy import isscalar, zeros, where, diag
    from nnspsolver import nnspsolver
    n0 = data.shape[0]
    E0 = zeros((n0,n0))
    A = E0.copy()
    D = data[:,:sample_size]
    X1 = data[:,1:(sample_size+1)]
    X1 = (X1.T).reshape(sample_size*n0)
    if SB == [] and SB_Index == []:
        Index = E0.copy()
        if isscalar(SM):
            A0 = zeros((n0,n0))
            A0[:,:] = 1
        else:
            A0 = SM
        indices = where(A0.T)
        Index[indices[1],indices[0]] = range(len(indices[0]))
        Index = Index.astype('int')
        n = len(indices[0])
        X0 = zeros((sample_size*n0,n))
        for k in range(n):
                E0[indices[1][k],indices[0][k]] = 1
                X0[:,k] = ((E0@D).T).reshape(sample_size*n0)
                E0[indices[1][k],indices[0][k]] = 0
        C = zeros((n0,n))
        for l in range(n0):
            C[l,Index[where(A0[:,l]),l]] = 1
        Mr = zeros((n+n0,n))
        Mr[:n,:n] = X0.T@X0
        Mr[n:,:n] = C
        rhs = zeros((n+n0,1))
        rhs[:n,0] = X0.T@X1
        rhs[n:,0] = 1
    else:
        Index = SB_Index
        n = len(SB_Index)
        m0 = len(SB)
        X0 = zeros((sample_size*n0,m0))
        k = 0
        for E0 in SB:
           X0[:,k] = ((E0@D).T).reshape(sample_size*n0)
           k+=1
        C = zeros((n,m0))
        for l in range(n):
            C[l,SB_Index[l]] = 1
        Mr = zeros((n+m0,m0))
        Mr[:m0,:m0] = X0.T@X0
        Mr[m0:,:m0] = C
        rhs = zeros((n+m0,1))
        rhs[:m0,0] = X0.T@X1
        rhs[m0:,0] = 1
    a = nnspsolver(Mr,rhs,delta=error,tol=error,nz=len(rhs))
    if SB == []:
        A[indices[1],indices[0]] = a[:n,0]
    else:
        k = 0
        for j in SB:
            A = A + a[k]*j
            k+=1
    sumA=sum(A)
    sumA[where(sumA==0)[0]]=1
    A = A@diag(1/sumA)
    return A