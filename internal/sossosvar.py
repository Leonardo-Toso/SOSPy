# SOSSOSVAR --- Declare a new sum of squares variable in 
#       an SOS program 
#
# [SOSP,VAR] = sossosvar(SOSP,Z)
#
# SOSP is the sum of squares program.
# VAR is the new sum of squares variable.
# ZSym is the column vector of monomials contained in the Gram 
# decomposition of VAR, i.e.,
#
#       VAR = ZSym' * COEFF * ZSym
#
# where COEFF is a coefficient matrix that is restricted to be 
# positive semidefinite. COEFF and the decision variables contained 
# in it will be constructed automatically by SOSSOSVAR.
#
# Both VAR and Z are either symbolic or polynomial objects.
#
# [SOSP,VAR] = sospolyvar(SOSP,ZSym,'wscoeff') will create
#the decision variables corresponding to VAR (i.e., coeff_xxx)
# also in MATLAB workspace.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.




import sympy as sym
import numpy as np

def sossosvar(prog,Z):
    
    var_holder = {}
    decvartable=[]
    for i in range(len(Z)**2):
        var_holder['coeff_' + str(i)]=sym.Symbol('coeff_%d'%prog.j)
        decvartable=np.append(decvartable,var_holder['coeff_' + str(i)])
        prog.j=prog.j+1
    locals().update(var_holder)
    
    var=0
    k=0
    for i in range(len(Z)):
        for j in range(len(Z)):
            var=var + Z[i]*(Z[j]*var_holder['coeff_' + str(k)])
            k=k+1
            
        
    
    prog.ndecvars=len(decvartable)
    prog.decvartable=np.concatenate((prog.decvartable,decvartable),axis=None)
    prog.decvar=np.concatenate((prog.decvar,decvartable),axis=None)

    
    return prog,var