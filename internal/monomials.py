# MONOMIALS --- Construct a vector of monomials with 
#       prespecified degrees.
#
# Z = monomials(VARTABLE,MINDEGREE,MAXDEGREE)
#
# Given a vector of independent variables VARTABLE (can be either
# symbolic or polynomial objects) and values for the minimum and maximum degree,
# this function constructs a column vector 
# containing all possible monomials of degree described in MINDEGREE and MAXDEGREE.
#
# For example, monomials([x1,x2],1,3) with x1 and x2 
# being symbolic variables will return monomials in x1 and 
# x2 of degree 1, 2, and 3.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.





import numpy as np
import symengine
import sympy as sym
from sympy import *
from sympy.polys.monomials import itermonomials
from sympy.polys.orderings import monomial_key

def monomials(vartable,mindeg,maxdeg):
    
   
    if mindeg!=maxdeg:
        M1=sorted(itermonomials(vartable, mindeg, 0), key=monomial_key('grlex', vartable[::-1] ))
        M2=sorted(itermonomials(vartable, maxdeg, mindeg), key=monomial_key('grlex', vartable[::-1] ))
    else: 
        M1=sorted(itermonomials(vartable, mindeg, mindeg), key=monomial_key('grlex', vartable[::-1] ))
        M2=sorted(itermonomials(vartable, maxdeg, mindeg), key=monomial_key('grlex', vartable[::-1] ))
    M1temp=[]
    for i in range(len(M1)):
        if sympify(M1[i]).is_number==False:
            sum_degree=sum(degree_list(M1[i]))
        elif sympify(M1[i]).is_number==True:
            sum_degree=M1[i]
        if sum_degree>=mindeg:
            M1temp=np.append(M1temp,M1[i])
    M1temp=list(M1temp)
    
    M=list(np.append(M1temp,M2))
    M=list(dict.fromkeys(M))
    M=sorted(M, key=monomial_key('grlex', vartable[::-1] ))
    
    
    return M