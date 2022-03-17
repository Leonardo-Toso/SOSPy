# FINDMATRIXSOS --- Find a sum of squares decomposition of a given matrix polynomial.

# [Q,Z,decomp] = findsos(P,SOLVER)

# P is a symmetric polynomial matrix.

# SOLVER is an optional argument which specifies the solver
#the default value for options.solver is 'mosek'.

# A positive semidefinite Q and a symbolic monomial vector Z will be
# computed such that

#    (Ir kron Z)' * Q * (Ir kron Z) = P(x)

# If P is not a sum of squares, the function will return empty Q and Z.

# If P is a polynomial with integer coefficients and is represented as a
# symbolic object, then [Q,Z] = findsos(P) will compute a rational
#matrix Q such that Z'*Q*Z = P.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import numpy as np
import sympy as sym 
from numpy import *
from sympy import *
from getequation import getequation
from polystr import polystr
from sossolve import sossolve
from getequationmc import getequationmc
from scipy.linalg import sqrtm
from sympy.physics.quantum import TensorProduct


def findsos(p,options='mosek'):
    p=expand(simplify(p))
    if np.asarray(p).shape==():
        
        symexpr,decvartable,varmat,vartable=polystr(p)
        At,b,Z,sos=getequation(symexpr,decvartable,varmat,vartable)
        sos.polynomial=p
        sos.flag_findsos=1
        sos.find=1
        sol=sossolve(sos,options)

        if len(sol.Q)==0:
            print('No sum of squares decomposition is found')
            decomp=[]
        else:
            print('Sum of squares decomposition found')
            L=sqrtm(sol.Q)
            decomp= L*(TensorProduct(eye(1),sol.Zmon))
        
    else: 
        decvartabletemp=[]
        vartabletemp=[]
        for i in range(p.shape[0]):
            for j in range(p.shape[1]):
                symexpr,decvartable,varmat,vartable=polystr(p[i,j])
                
                decvartabletemp=np.concatenate((decvartabletemp,decvartable),axis=None)
                vartabletemp=np.concatenate((vartabletemp,vartable),axis=None)
        
        vartable=list((Matrix(np.unique(np.array(vartabletemp).astype('str').tolist()))))
        decvartable=list((Matrix(np.unique(np.array(decvartabletemp).astype('str').tolist()))))
        
        At,b,Z,sos=getequationmc(p,decvartable,vartable)
        
        sos.matrixtype=1
        sos.polynomial=p
        sos.flag_findsos=1
        sos.find=1
        sol=sossolve(sos,options)
        if len(sol.Q)==0:
            print('No sum of squares decomposition is found')
            decomp=[]
        else:
            print('Sum of squares decomposition found')
            L=sqrtm(sol.Q)
            decomp= L*(TensorProduct(eye(p.shape[0]),sol.Zmon))
        
        
    
    
    
    return sol.Q,sol.Zmon,decomp


