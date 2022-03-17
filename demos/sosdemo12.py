# SOSDEM12 --- MAX CUT

import sympy as sym
import numpy as np
import math 
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
from monomials import monomials
from sossosvar import sossosvar

#Variables
x1,x2=sym.symbols('x1 x2 x3 x4 x5')
vartable=[x1,x2,x3,x4,x5]

#Number of cuts
f = 2.5 - 0.5*x1*x2 - 0.5*x2*x3 - 0.5*x3*x4 - 0.5*x4*x5 - 0.5*x5*x1

#Boolean constraints

var_holder = {}
for i in range(5):
    var_holder['bc' + str(i)]=0
locals().update(var_holder)

var_holder['bc' + str(0)]=x1**2 - 1
var_holder['bc' + str(1)]=x2**2 - 1
var_holder['bc' + str(2)]=x3**2 - 1
var_holder['bc' + str(3)]=x4**2 - 1
var_holder['bc' + str(4)]=x5**2 - 1

#=============================================
# First, initialize the sum of squares program

program=sosprogram(vartable)


# =============================================
# Then define SOSP variables

# -- p1(x) -- : sum of squares
# Monomial vector: 5 independent variables, degree <= 1

var_holder2 = {}
for i in range(6):
    var_holder2['p' + str(i)]=0
locals().update(var_holder2)

Z=monomials(vartable,1,0)
program,var_holder2['p' + str(0)]=sossosvar(program,Z)

#-- p2(x) ... p6(x) : polynomials
# Monomial vector: 5 independent variables, degree <= 2

Z=mon(vartable,2,0)
for i in range(5):
    program,var_holder2['p' + str(i+1)] = sospolyvar(program,Z)
    
    
# =============================================
# Next, define SOSP constraints

# Constraint : p1(x)*(gamma - f(x)) +  p2(x)*bc1(x)
#               + ... + p6(x)*bc5(x) - (gamma-f(x))^2 >= 0

gamma = 4;

expr=expand(var_holder2['p' + str(0)]*(gamma-f))

for i in range(1,6):
    expr=expr+expand(var_holder2['p' + str(i)]*var_holder['bc' + str(i-1)])
                     
expr = expr - (gamma-f)**2
expr=expand(expr)

program=sosineq(program,expr)

# ===========================================
# And call solver

sol=sossolve(program)
      