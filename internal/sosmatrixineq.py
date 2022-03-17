#SOSMATRIXINEQ --- Creates a SOS constraint from a matrix inequality constraint
# [SOSP] = sosmatrixineq(SOSP,fM)
# SOSP is the sum of squares program.
# fM is a polynomial matrix used to generate the polynomial inequality y'fMy>0

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import numpy as np
import sympy as sym
from sympy import *

def sosmatrixineq(program,p,option=[]):
    program.nconstraint+=1
    (program.matexpr).append(p)
    if option==[]:
        program.typeofmatrixineq=np.concatenate((program.typeofmatrixineq,'quadraticMineq'),axis=None)
    else:
        program.typeofmatrixineq=np.concatenate((program.typeofmatrixineq,'Mineq'),axis=None)
        
    return program