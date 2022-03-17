# SOSDEMO6 --- Control design for Linear Hybrid System with periodic jump.


import sympy as sym
import numpy as np
from sympy import *
from findsos import findsos
from sosprogram import sosprogram
from sosineq import sosineq
from sosmatrixineq import sosmatrixineq
from soseq import soseq
from sossolve import sossolve
from monomials import monomials 
from sospolyvar import sospolyvar
from sosgetsol import sosgetsol
from sospolymatrixvar import sospolymatrixvar


#Variables 

xi,x1,x2,x3,x4=sym.symbols('xi x1 x2 x3 x4')
vartable = [xi,x1,x2,x3,x4]

#State represetation for a given Linear Hybrid System

A1=Matrix([[-1,0],[0,1]])
B=Matrix([[1],[0]])
E1=Matrix([[2,0],[1,0.5]])
F=Matrix([[0],[0]])
v=Matrix([[x1],[x2]])
v1=Matrix([[x1],[x2],[x3],[x4]])

#Convex weighting factor

alpha=0.1

#Jump period
tau=1


#Tolarance
eps=1e-4

#=============================================
# First, initialize the sum of squares program

program=sosprogram(vartable)
program.newtonpolytope=1

# P matrix

Z1=monomials([xi],0,0)
program,P= sospolymatrixvar(program,Z1,[2,2],'symmetric')

# L1 matrix

Z2=monomials([xi],0,0)
program,L1= sospolymatrixvar(program,Z2,[1,2])

#L2 matrix 

Z3=monomials([xi],0,2)
program,L2= sospolymatrixvar(program,Z3,[2,2])

# N matrix

Z4=monomials([xi],0,0)
program,N= sospolymatrixvar(program,Z4,[4,4],'symmetric')

#Lambda matrix

Lambda1=np.concatenate(((2*P*A1.T)+(2*alpha*L1.T*B.T)+(1/(tau))*E1*P*E1.T+(((1-alpha)))*(E1*L2.T+L2*E1.T)-(1/tau)*P,(1-alpha)*L2.T),axis=1)
Lambda2=np.concatenate(((1-alpha)*L2,-(P)),axis=1)
Lambda=Matrix(np.concatenate((Lambda1,Lambda2),axis=0))

#=============================================
# Next, define the inequalities

program=sosineq(program,(v.T*(P-eps*eye(2))*v)[0])

program=sosineq(program,(v1.T*(-Lambda-N*(tau-xi)*xi)*v1)[0])

program=sosineq(program,(v1.T*(N-eps*eye(4))*v1)[0])

                       
# =============================================
# And call solver
sol = sossolve(program)