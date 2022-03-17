#Collect symexpr = g(c)*h(x) + g0 where x are the variables in vartable,
#c is composed of decision variables, and h(x) concerns the vector of monomials.
#This function is used to find unique monomials in the polynomial variables.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
from typing import NamedTuple
import sympy as sym
from sympy import LT,LM,degree,degree_list
from dataclasses import dataclass
from sympy.matrices import Matrix, eye, zeros, ones, diag, Transpose
    
    
def collect(symexpr):

    #Get the matrix h
    
    ptemp=symexpr.polynomial
    h=zeros(symexpr.nterms,1)
    j=0
    for j in range(symexpr.nterms):
        h[j]=LM(ptemp)
        ptemp=ptemp-LT(ptemp)
    
    #Get the matrix g
    
    ptemp=symexpr.polynomial
    g=zeros(symexpr.nterms,1)
    j=0
    for j in range(symexpr.nterms):
        mon_coeff=LT(ptemp)
        mon=LM(ptemp)
        g[j]= mon_coeff/mon
        ptemp=ptemp-LT(ptemp)
        
    #Get the value of g0
    
    g0=zeros(1,1)
    g0[0]=symexpr.coeff0
    
    return g,h,g0