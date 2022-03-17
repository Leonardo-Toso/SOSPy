#SOSDEMO8 --- Bounds in Probability

import sympy as sym
import numpy as np
from findsos import findsos
from sosprogram import sosprogram
from sosineq import sosineq
from sossolve import sossolve
from sossetobj import sossetobj
from sosgetsol import sosgetsol



#Variables
x,a,b,c=sym.symbols('x a b c')

#The probability adds up to one.
m0 = 1

#Mean
m1  = 1

#Variance
sig = 1/2

# E(x^2)
m2 = sig**2+m1**2

#Support of the random variable
R = [0,5]

#Event whose probability we want to bound
E = [4,5]

# =============================================
# Constructing and solving the SOS program

program=sosprogram([x],[a,b,c])

P = a + b*x + c*x**2 

#Nonnegative on the support
program=sosineq(program,P,R)

#Greater than one on the event
program=sosineq(program,P-1,E)

#The bound 

bnd =  a * m0 + b * m1 + c * m2

# Objective: minimize the bound

program= sossetobj(program, bnd)

# =============================================
# And call solver

sol=sossolve(program)

# =============================================
# Get solution

SOLBND=sosgetsol(sol,bnd)
SOLP=sosgetsol(sol,P)