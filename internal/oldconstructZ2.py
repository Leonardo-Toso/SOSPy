# Old construction II of the vector of monomials Z

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import numpy as np
import sympy as sym
from sympy import degree_list,degree,LT,LM,LC,sympify
from sympy import matrix2numpy
from sympy.core.compatibility import as_int
from getequation import getequation
from sympy import ImmutableMatrix as  Matrix
from sympy.matrices import Matrix, eye, zeros, ones, diag, Transpose
from scipy.sparse import csr_matrix
from max_kron import max_kron 
from collect import collect
from sympy import poly
from sympy import factor
from numpy import matlib

def oldconstructZ2(vartable,d):
    
    

    ZZ=np.zeros((1, vartable)) 
    
    for i in range(vartable):
        ss=ZZ.shape[0]
        
        ZZ = np.matlib.repmat(ZZ,d+1,1)
        
        for j in range(d+1):
            
            for l in np.arange(ss*j,ss*j+ss):
               
                ZZ[l][i]=j
                
        sum_rows=ZZ.sum(axis=1)   
        nelements=len(sum_rows)
        index=-np.ones((1,nelements))
       
        k=0
        for j in range(nelements):    # Throw away invalid monomials
            if (sum_rows[j]<=d):
                index[0][k]=j
                k=k+1
        
        idx=-np.ones((1,k))
        idx=index[0][0:k]
        idx=idx.astype(int)
        
        Ztemp=np.zeros((len(idx),ZZ.shape[1]))
        
        for m in range(len(idx)):
            
            Ztemp[m][:]=ZZ[idx[m]][:]
        
        ZZ=Ztemp
        
        
    sum_rows=ZZ.sum(axis=1)   
    nelements=len(sum_rows)
    index=-np.ones((1,nelements))

    k=0
    for j in range(nelements):    # Throw away invalid monomials
        if (sum_rows[j]==d):
            index[0][k]=j
            k=k+1

    idx=-np.ones((1,k))
    idx=index[0][0:k]
    idx=idx.astype(int)

    Ztemp=np.zeros((len(idx),ZZ.shape[1]))

    for m in range(len(idx)):

        Ztemp[m][:]=ZZ[idx[m]][:]

    Z=Ztemp
       
    return Z