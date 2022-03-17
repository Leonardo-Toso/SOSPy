# SOSDEMO2 --- Lyapunov Function Search 

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

#Variables
x1,x2,x3=sym.symbols('x1 x2 x3')
vartable=[x1,x2,x3]

#Constructing the vector field dx/dt = f

f=Matrix([[-x1**3-x1*x3**2],
          [-x2-x1**2*x2],
   [-x3+3*x1**2*x3-3*x3/(x3**2+1)]])

# =============================================
# First, initialize the sum of squares program

program=sosprogram(vartable)

# =============================================
# The Lyapunov function V(x): 

program,V=sospolyvar(program,[x1**2,x2**2,x3**2])

#Constraint 1 : V(x) - (x1^2 + x2^2 + x3^2) >= 0

program=sosineq(program,V-(x1**2+x2**2+x3**2))

#Constraint 2: -dV/dx*(x3^2+1)*f >= 0

expr=-(diff(V,x1)*f[0] + diff(V,x2)*f[1]+ diff(V,x3)*f[2])*(x3**2+1)

program=sosineq(program,expr)

# =============================================
# And call solver
sol=sossolve(program)

# =============================================
# Finally, get solution

SOLV=sosgetsol(sol,V)