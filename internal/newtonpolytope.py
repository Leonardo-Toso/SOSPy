#NEWTONPOLYTOPE--- Performing the Newton polytope algorithm to reduce the space of candidate monomials S of a polytope P

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import numpy as np
import mosek
from  mosek.fusion import *
import sys

def newtonpolytope(S,P):
    dim=S.shape[1]
    Sp=np.zeros((1,S.shape[1]))
    for i in range(S.shape[0]):
        ci=(np.append(S[i,:],-1)).T
        c=np.zeros((ci.shape[0],1))
        for k in range(c.shape[0]):
            c[k,0]=ci[k]
        si=c
        c=Matrix.dense(c)
        
        with Model("LP") as M:
            X=M.variable('X',dim+1)
            M.objective(ObjectiveSense.Maximize, Expr.dot(c,X))
            for j in range(P.shape[0]):
                pktemp=np.append(0.5*P[j,:],-1)
                pk=np.zeros((1,pktemp.shape[0]))
                for k in range(pk.shape[1]):
                    pk[0,k]=pktemp[k]
                
                M.constraint(Expr.mul(pk, X),Domain.lessThan(0.0))
            
            M.constraint(Expr.mul(si.T, X),Domain.greaterThan(0.0))
            M.solve()
            info=M.getProblemStatus()
       
            #1) and 3) Cases: Infeasible problem or optimal value equal zero-the monomial must be kept

            if info==ProblemStatus.PrimalAndDualFeasible:
                
                Scol=np.zeros((1,S[i,:].shape[0]))
                for k in range(Scol.shape[1]):
                    Scol[0,k]=(S[i,:])[k]
                if i==0:
                    Sp=Scol
                else:
                    Sp=np.concatenate((Sp,Scol),axis=0)
    
    return Sp