# SOSDEMO4 --- Matrix Copositivity

import sympy as sym
import numpy as np
from sympy import *
from findsos import findsos
from sosprogram import sosprogram
from sosineq import sosineq
from soseq import soseq
from sossolve import sossolve
from sossetobj import sossetobj
from sospolyvar import sospolyvar
from sosgetsol import sosgetsol
from sosdecvar import sosdecvar 

#Variables
x1,x2,x3,x4,x5=sym.symbols('x1 x2 x3 x4 x5')
vartable=[x1,x2,x3,x4,x5]

#The matrix under consideration

J =Matrix([[1, -1,  1,  1, -1],
          [-1,  1, -1,  1,  1],
          [ 1, -1,  1, -1,  1],
          [ 1,  1, -1,  1, -1],
          [-1,  1,  1, -1, 1]])

#=============================================
# First, initialize the sum of squares program

program=sosprogram(vartable)

# =============================================
# Next, define SOSP constraints

# Constraint : r(x)*J(x) - p(x) = 0

J=Matrix([x1**2, x2**2, x3**2, x4**2, x5**2]).T*J*Matrix([x1**2, x2**2, x3**2, x4**2, x5**2])
r= x1**2 + x2**2 + x3**2 + x4**2 + x5**2

program=sosineq(program,r*J[0])

# =============================================
# And call solver

sol=sossolve(program)