# SOSDEMO3 --- Bound on Global Extremum

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
x1,x2,gam=sym.symbols('x1 x2 gam')
vartable=[x1,x2]

# =============================================
# First, initialize the sum of squares program

program=sosprogram(vartable)

#=============================================
# Declare decision variable gam too

program=sosdecvar(program,[gam])

#=============================================
# Next, define SOSP constraints

#Constraint : r(x)*(f(x) - gam) >= 0
# f(x) is the Goldstein-Price function

f1 = x1+x2+1
f2 = 19-14*x1+3*x1**2-14*x2+6*x1*x2+3*x2**2
f3 = 2*x1-3*x2
f4 = 18-32*x1+12*x1**2+48*x2-36*x1*x2+27*x2**2

f = (1+f1**2*f2)*(30+f3**2*f4)

program=sosineq(program,f-gam)

# =============================================
# Set objective : maximize gam
program=sossetobj(program,-gam)

# =============================================
# And call solver

sol=sossolve(program)

# =============================================
# Finally, get solution

SOLgamma=sosgetsol(sol,gam)