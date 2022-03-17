#SOSDEMO10 --- Set containment

import sympy as sym
import numpy as np
from sympy import *
from findsos import findsos
from sosprogram import sosprogram
from sosineq import sosineq
from sossolve import sossolve
from sossetobj import sossetobj
from sosgetsol import sosgetsol
from monomials import monomials
from sospolymatrixvar import sospolymatrixvar 
from sosmatrixineq import sosmatrixineq

#Variables
x1,x2=sym.symbols('x1 x2')

vartable=[x1,x2]

eps = 1e-6

# =============================================
# This is the problem data

p = x1**2+x2**2
gamma = 1
g0 = Matrix([[2, 0]])*Matrix([[x1],[x2]])
theta = 1

# =============================================
# Initialize the sum of squares program

program=sosprogram(vartable)

# =============================================
# The multiplier

Zmon=monomials(vartable,0,4)

program,s= sospolymatrixvar(program,Zmon,[1,1])


# =============================================
# Term to be added to g0

Zmon1=monomials(vartable,2,3)

program,g1= sospolymatrixvar(program,Zmon1,[1,1])

#=============================================
# The expression to satisfy the set containment

Sc = Matrix([[theta**2-(s*(gamma-p))[0], g0[0] + g1[0]],[g0[0]+g1[0], 1]])
expr=Sc-eps*eye(2)

program=sosmatrixineq(program,expr)
# =============================================
# And call solver
sol=sossolve(program)


# =============================================
# Get solution
s = sosgetsol(sol,s);
g1 = sosgetsol(sol,g1);