# SOSINEQ --- Add a new inequality constraint f(x) >= 0
#    to an SOS program 
#
# SOSP = sosineq(SOSP,EXPR)
#
# SOSP is the sum of squares program.
# EXPR is the expression on the left hand side of the constraint, i.e., f(x).
#
# EXPR can be a column vector. In this case, several inequality
# constraints will be added simultaneously to the sum of
# squares program.
#
# SOSP = sosineq(SOSP,EXPR,[a b]) is used to specify the interval 
# a <= x <= b where f(x) has to be non-negative. Currently, this is 
# possible only if the sum of squares program is univariate. The lower
# limit a can be -Inf, and similarly b can be +Inf.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import numpy as np
import sympy as sym
from sympy import *

def sosineq(program,p,interval=[]):
    program.nconstraint+=1
    (program.ineqexpr).append(p)
    if interval!=[]:
        program.intervalexpr=np.concatenate((program.intervalexpr,interval))
        program.intervalchar=np.concatenate((program.intervalchar,'interval'),axis=None)
    else:
        program.intervalchar=np.concatenate((program.intervalchar,'nointerval'),axis=None) 
    return program