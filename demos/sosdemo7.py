#SOSDEMO7 --- Chebyshev polynomials

import sympy as sym
import numpy as np
from findsos import findsos
from sosprogram import sosprogram
from sosineq import sosineq
from sossolve import sossolve
from sossetobj import sossetobj
from sospolyvar import sospolyvar
from sosgetsol import sosgetsol
from sosdecvar import sosdecvar 
from monomials import monomials

ndeg = 8;   # Degree of Chebyshev polynomial

#Variables

x,gam=sym.symbols('x gam')

#=============================================
# First, initialize the sum of squares program

program=sosprogram([x],[gam])
Z=monomials([x],0,ndeg-1)
program,P1 = sospolyvar(program,Z)
P=P1 + gam*x**ndeg

#Imposing the inequalities

program=sosineq(program,1-P,[-1,1])
program=sosineq(program,1+P,[-1,1])

#And setting objective
program= sossetobj(program, -gam)

# =============================================
# And call solver
sol=sossolve(program)

# =============================================
# Finally, get solution

SOLP=sosgetsol(sol,P)
SOLgam=sosgetsol(sol,gam)