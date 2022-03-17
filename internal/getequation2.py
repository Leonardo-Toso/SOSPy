# GETEQUATION II --- Convert a symbolic expression to At, b, and Z  used in an SOS program.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import numpy as np
from typing import NamedTuple
import symengine
import sympy as sym
from sympy import LT,LM,degree,degree_list
from collect import collect 
from dataclasses import dataclass
from sympy.matrices import Matrix, eye, zeros, ones, diag, Transpose
from scipy.sparse import csr_matrix
from ismember import ismember
from sortNoRepeat import sortNoRepeat
from max_kron import max_kron
from polystr import polystr

def getequation2(symexpr,decvartable,varmat,vartable):
    
     # Handle Polynomial Objects 

    decvarnum=len(decvartable)
    vartable=np.concatenate((vartable,varmat))
    cvartable=np.array(vartable).astype('str').tolist()
    cvartable=np.asarray(cvartable)
    cvartable=np.sort(cvartable)
    dimp=1 
    

   
    g,h,g0=collect(symexpr)
    
    #Then the polynomial expression is decomposed in three terms, such that p=g*h+g0
    
    h0=zeros(1,1)
    h0[0]=1
    if g0[0] != 0:
        g=g.col_join(g0)
        h=h.col_join(h0)
        
    #Define the numbers of mononials and variables in p
    
    nmon=symexpr.nterms
    nvar=symexpr.nvars
    
    #Reorder the monomials matrix with variables in the order listed in cvartable and sorted monomials

    #Create the matrix Z, which will receive the monomials of the expression
    
    Z=zeros(nmon,nvar)
    Z=symexpr.degmat
    
    #Creating the matrix Z in the multivariable case
    
    if len(vartable)!=1:
        polynomial=symexpr.polynomial+symexpr.coeff0
        for i in range(len(decvartable)):
            polynomial=polynomial.subs(decvartable[i],1)
        symexprg,decvartableg,varmatg,vartableg=polystr(polynomial)
        if symexprg.coeff0!=0:
            nterms=symexprg.nterms+1
        else:
            nterms=symexprg.nterms
        Z=np.zeros((nterms,len(vartable)))
        for k in range(symexprg.nterms):
            mon=LT(polynomial)
            for j in range(len(vartable)):
                deg=degree(mon,gen=vartable[j])
                Z[k][j]=deg
            polynomial=polynomial-LT(polynomial)


    #The vector b will receive the coefficients from the symbolic expression

    b=g.T

    
    
    return b,Z