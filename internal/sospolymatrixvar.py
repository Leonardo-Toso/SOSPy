#SOSPOLYMATRIXVAR--Generate symbolic matrices according to the specified dimension, monomials and properties.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import sympy as sym
import numpy as np
from sympy import *
from sospolyvar import sospolyvar

def sospolymatrixvar(prog,monomials,dim,option=[]):
    
    M=-ones(dim[0],dim[1])
    if option!='symmetric' and option!='diagonal':
        for i in range(dim[0]):
            for j in range(dim[1]):
                prog,M[i,j]=sospolyvar(prog,monomials)
                
    elif option=='symmetric' : #Symmetric Matrix Case
        for i in range(dim[0]):
            for j in range(dim[1]):
                if i<=j:
                    prog,M[i,j]=sospolyvar(prog,monomials)
                
        for i in range(dim[0]):
            for j in range(dim[1]):
                if M[i,j]==-1:
                    M[i,j]=M[j,i]
                
    elif option=='diagonal' : #Diagonal Matrix Case
        for i in range(dim[0]):
            for j in range(dim[1]):
                if i==j:
                    prog,M[i,j]=sospolyvar(prog,monomials)
                if M[i,j]==-1:
                    M[i,j]=0
    return prog, M


