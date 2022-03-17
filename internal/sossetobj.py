# SOSSETOBJ --- Set the objective function of an SOS program. 
# SOSP = sossetobj(SOSP,EXPR)
# Given a sum of squares program SOSP and an expression EXPR, SOSSETOBJ
# makes EXPR the objective function of the sum of squares program SOSP, 
# i.e., EXPR will be minimized.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import sympy as sym
import numpy as np
from sosconstraint import sosconstraint

def sossetobj(program,obj):
    program=sosconstraint(program)
    sos=program.sos
    if sos.nconstraint>1:
        for i in range(len(sos.decvartable2)):
            if obj.coeff(sos.decvartable2[i])!=0:
                sos.ccon[i][0]=obj.coeff(sos.decvartable2[i])
    
    elif sos.nconstraint==1:
        sos.obj=np.append(sos.obj,obj)
    program.sos=sos
    program.obj=1
    return program