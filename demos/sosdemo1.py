#SOSDEMO1 --- Sum of Squares Test

import sympy as sym
import numpy as np
from findsos import findsos
from sosprogram import sosprogram
from sosineq import sosineq
from soseq import soseq
from sossolve import sossolve


x1,x2=sym.symbols('x1 x2')
vartable = [x1, x2]

#====================================
# First, initialize the sum of squares program

program=sosprogram(vartable)   # No decision variables.

#=============================================
# Next, define the inequality

# p(x1,x2) >=  0
p = 2*x1**4 + 2*x1**3*x2 - x1**2*x2**2 + 5*x2**4
program=sosineq(program,p)

# =============================================
# And call solver
sol = sossolve(program)