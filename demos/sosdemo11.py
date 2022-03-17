# SOSDEMO11 --- Upper bound for the structured singular value mu

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
x1,x2=sym.symbols('x1 x2 x3 x4 x5 x6 x7 x8')
vartable=[x1,x2,x3,x4,x5,x6,x7,x8]

#The matrix under consideration

alpha = 3 + math.sqrt(3)
beta = math.sqrt(3) - 1
a = math.sqrt(2/alpha)
b = 1/math.sqrt(alpha)
c = b
d = -math.sqrt(beta/alpha)
f = (1 + I)*math.sqrt(1/(alpha*beta))
U = Matrix([[a, 0],[b, b], [c, I*c], [d, f]])
V = expand(Matrix([[0, a],[b, -b],[c, -I*c], [-I*f, -d]]))
M =expand(U*conjugate(V).T)

#Constructing A(x)'s
gam = 0.8724

Z = monomials(vartable,1,1)
Z=Matrix(Z)
var_holder = {}
for i in range(4):
    var_holder['A' + str(i)]=0
locals().update(var_holder)

for i in range(4):
    s=np.zeros((4,4))
    s[i][i]=1
    H=expand(expand(conjugate(M[i,:]).T*M[i,:])- (gam**2)*s)
    H = expand(sym.BlockMatrix([[re(H), -im(H)],[im(H), re(H)]]))
    H=Matrix(H)
    var_holder['A' + str(i)]=expand(Z.T*H*Z)[0] 
    
#=============================================
# First, initialize the sum of squares program

program=sosprogram(vartable)

# -- Q(x)'s -- : sums of squares
# Monomial vector: [x1; ... x8]
var_holder2 = {}
for i in range(4):
    var_holder2['Q' + str(i)]=0
locals().update(var_holder2)

for i in range(4):
    program,var_holder2['Q' + str(i)]=sossosvar(program,Z)   
    
#-- r's -- : constant sum of squares

Z=monomials(vartable,0,0)
r=zeros(4,4)
for i in range(4):
    for j in range(i+1,4):
        program,r[i,j]=sossosvar(program,Z)   
        
# =============================================
# Next, define SOSP constraints
#Constraint : -sum(Qi(x)*Ai(x)) - sum(rij*Ai(x)*Aj(x)) + I(x) >= 0

expr=0
#Adding term

for i in range(4):
    expr = expr - expand(var_holder['A' + str(i)]*var_holder2['Q' + str(i)])
for i in range(4):
    for j in range(i+1,4):
        expr = expr - expand(var_holder['A' + str(i)]*var_holder['A' + str(j)]*r[i,j])

#Constant term: I(x) = -(x1^4 + ... + x8^4)
Ix=Matrix(vartable)
Ixsum=0
for i in range(len(vartable)):
    Ix[i,0]=-Ix[i,0]**4
    Ixsum=Ixsum+Ix[i,0]
expr=expr+Ixsum
expr=expand(expr)
    
program=sosineq(program,expr)

# =============================================
# And call solver

sol=sossolve(program)