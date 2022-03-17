# SOSCONSTRAINT---Processing the constraints from the SoS program. Inequality, equality and matrix 
# constraints are processed in this function.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
from ineq import ineq
from eq import eq
from matrixineq import matrixineq

def sosconstraint(program):
    
    #Inequality constraints
    
    k=0
    for i in range(len(program.ineqexpr)):
        if program.intervalchar[i]=='interval':
            program=ineq(program,program.ineqexpr[i],[program.intervalexpr[k],program.intervalexpr[k+1]])
            k=k+2
        else:
            program=ineq(program,program.ineqexpr[i])
            
    
    #Equality constraints
    
    for i in range(len(program.eqexpr)):
        program=eq(program,program.eqexpr[i])
    
    
    #Matrix constraints
    
    for i in range(len(program.matexpr)):
        if program.typeofmatrixineq[i]=='quadraticMineq':
            program=matrixineq(program,program.matexpr[i],[])
        else:
            program=matrixineq(program,program.matexpr[i],'Mineq')
    
   
    return program