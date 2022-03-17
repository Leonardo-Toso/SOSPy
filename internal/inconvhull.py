#INCONVHULL --- Returns in Z3 the elements of Z1 that are inside the convex hull of Z2

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
import sympy as sym
from sympy import *
from numpy import *
from scipy.linalg import orth
from scipy import linalg
from useconvhulln import useconvhulln
from scipy.linalg import null_space
from sympy import  Matrix

def inconvhull(Z1,Z2):
    
    #First, find the affine subspace where everything lives
    #(for instance, in the homogeneous case)

    nmons=Z2.shape[0]
    
    #Translate so it goes through the origin
    
    mr=np.asarray(Z2.mean(0))
    
    Rzero=Z2 - np.matlib.repmat(mr,nmons,1)
    
    #The columns of N generate the subspace 
    
    N=null_space(Rzero)
    
    
    # Z2*N should be constant 
    cval=np.asarray(np.matmul(Z2,N).mean(0))
    
    #Get only the monomials in the subspace
    
    tol=0.01
    
    sum_ix=np.sum(abs(np.matmul(Z1,N) - np.matlib.repmat(cval,Z1.shape[0],1)),axis=1)
    ix=[]
    for i in range(len(sum_ix)):
        if sum_ix[i]<tol:
            ix=np.concatenate((ix,i),axis=None)
         
    nZ1 = Z1[ix.astype(int),:]
    
    # Now, the inequalities:
    # Find an orthonormal basis for the subspace
    # (I really should do both things at the same time...)
    # Project to the lower dimensional space, so convhull works nicely
    
   
    Q=orth(Rzero.T)
    
    if (matmul(Z2,Q)).shape[1]>1:
        A,B=useconvhulln(matmul(Z2,Q))
        
        
        #Find the ones that satisfy the inequalities, and keep them.
       
        ix_temp= (np.matlib.repmat(B,1,nZ1.shape[0]) -matmul(matmul(A,Q.T),nZ1.T)).min(0)
        ix=[]
        for i in range(len(ix_temp)):
            if ix_temp[i]>-tol:
                ix=np.concatenate((ix,i),axis=None)
                
        Z3=nZ1[ix.astype(int),:]
        
        
    elif (matmul(Z2,Q)).shape[1]==1:
        A=np.array([[1],[-1]])
        B=np.array([max(matmul(Z2,Q)),-min(matmul(Z2,Q))])
        
        
        ix_temp= (np.matlib.repmat(B,1,nZ1.shape[0]) -matmul(matmul(A,Q.T),nZ1.T)).min(0)
        
        ix=[]
        for i in range(len(ix_temp)):
            if ix_temp[i]>-tol:
                ix=np.concatenate((ix,i),axis=None)
             
        Z3=nZ1[ix.astype(int),:]
      
    else:
        Z3=nZ1
        
        
        
    Z3=np.unique(Z3, axis=0)
    
    return Z3
