# GETCONSTRAINT --- Find constraint for sum of squares decomposition.

# [A,ZZ] = getconstraint(Z)

# Z is a monomial vector description.
# This function computes the constraint matrix A and the polynomial
# vector ZZ, such that if q satisfies

#    A*q = F'

# Then
#    Z'*Q*Z = F*ZZ      (where Q is the matrix form of q)

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
import symengine
import sympy as sym
from sympy import degree_list,degree,LT,LM,LC,sympify
from sympy import matrix2numpy
from sympy.core.compatibility import as_int
from getequation import getequation
from sympy import ImmutableMatrix as  Matrix
from sympy.matrices import Matrix, eye, zeros, ones, diag, Transpose
from max_kron import max_kron 
from collect import collect
from sympy import poly
from sympy import factor
from monomials import monomials

def getconstraint(Z):
    
    sizeZ=Z.shape[0]
    ZZ=Z + np.matlib.repmat(Z[0],sizeZ,1).T
    ZZ=ZZ[0]
    
    #Creating a structure matrix M which will received the permutation matrices: R1, ..., Rn.
    
    Rn=[]    
    M = [Rn for i in range(sizeZ)]
    M[0]=np.eye(sizeZ, dtype=int)
    
    k=0
    R1=[]
    R2=np.eye(sizeZ, dtype=int)
    
    for i in np.arange(1,sizeZ):
        
        
        R1= np.concatenate((np.eye(sizeZ+k, dtype=int), np.zeros((sizeZ+k,1),dtype=int)), axis=1) 
        R2=np.concatenate((np.zeros((sizeZ,1),dtype=int),R2), axis=1)
        k=k+1
        ZZ=np.append(ZZ,ZZ[len(ZZ)-1]+1)
        for j in range(i):
            M[j] = np.matmul(M[j], R1)
         
        M[i]=R2
   
    Q=np.zeros((sizeZ,sizeZ))
    A=np.zeros((ZZ.shape[0],sizeZ**2))
    c=-1
    for i in range (sizeZ**2):
        r=np.mod(i,sizeZ)
        if r==0:
            c=c+1
        Q[r][c]=1
        j=r
        A[:,i]=(np.array([A[:,i]]).T+np.matmul(M[j].T,np.array([Q[j][:]]).T)).T[0]
        Q[r][c]=0

    
    return A,ZZ