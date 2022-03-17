# PREPROCESS--- Pre-processing the symbolic expression within a given interval.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import sympy as sym
import numpy as np
from sympy import *
import math

def preprocess(symexpr,var,interval):
    
#Get the maximum degree of the independent variable

    maxdeg = 0
    dummy=diff(symexpr,var[0])
    while dummy !=0:
        maxdeg = maxdeg+1
        dummy=diff(dummy,var[0])
        
#Substitute var

    if interval[1]==math.inf:
        newvar = var+interval[0]
        newsymexpr=symexpr.subs(var[0],newvar)
    elif interval[0]==-math.inf:
        newvar = -var+interval[1]
        newsymexpr=symexpr.subs(var[0],newvar)

    else:
        newvar = (interval[1]-interval[0])/2*(1-var[0])/(1+var[0]) + (interval[1]+interval[0])/2
        newsymexpr = symexpr.subs(var[0],newvar)*(1+var[0])**maxdeg
    
    return newsymexpr