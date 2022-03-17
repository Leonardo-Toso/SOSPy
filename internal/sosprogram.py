# SOSPROGRAM --- Initialize a new sum of squares program.
#
# SOSP = sosprogram(VARTABLE,DECVARTABLE)
#
# SOSP is a new sum of squares program. 
# VARTABLE is a vector of independent variables.
# DECVARTABLE is a vector of decision variables (optional).
#
# Both VARTABLE and DECVARTABLE are either symbolic or polynomial objects.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
import sympy as sym
from sympy import *


def sosprogram(var=[],decvar=[]):
 
    class program: 
        var=[]
        decvar=[]
        nvars=0
        ndecvars=0
        Q=[]
        Zmon=[]
        nconstraint=0
        Atcon=[]
        bcon=[]
        ccon=[]
        decvartable=[]
        k=0
        j=1
        ktemp=0
        Ks=[]
        polys=[]
        cons=[]
        extravaridx=[1]
        extravarnum=0
        interval=0
        sos=[]
        find=0
        ineqexpr=[]
        eqexpr=[]
        matexpr=[]
        intervalexpr=[]
        intervalchar=[]
        obj=0
        typeofmatrixineq=[]
        newtonpolytope=1
        
    program.var=var
    program.decvar=decvar
    program.nvars=len(var)
    program.ndecvars=len(decvar)
    program.decvartable=decvar
    program.vartable=var    
    return program
