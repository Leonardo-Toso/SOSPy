# SOSDEMO5 --- Conditional Convex Polynomial Fitting

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

#Variables
x,s,theta=sym.symbols('x s theta')
vartable = [x,s]
decvartable=[theta]

#=============================================
# First, initialize the sum of squares program
program=sosprogram(vartable,decvartable)
program.newtonpolytope=1

Z1=monomials([x],0,2)
program,f= sospolyvar(program,Z1)

g=x**2 -24

Z2=monomials(vartable,0,4)
program,w= sospolyvar(program,Z2)

#Generating the data
u=np.arange(0,7,1)

#Cost function 
a=[]
for i in range(len(u)):
     a=np.concatenate((a,f.subs(x,u[i]) - np.exp(u[i])),axis=None)
    
at=zeros(len(u),1)
for i in range(at.shape[0]):
    at[i]=a[i]


#Constraints
GDf=diff(diff(f,x),x)

C1=s*GDf*s - w*(1-g)
C2=w

C3_1=np.concatenate((theta,at.T),axis=None)
C3_1t=zeros(1,len(u)+1)
for i in range(C3_1t.shape[1]):
    C3_1t[i]=C3_1[i]


C3_2=np.concatenate((at,np.eye(len(u))),axis=1)

C3=Matrix(np.concatenate((C3_1t,C3_2),axis=0))
#=============================================
# Next, define the inequality

program=sosineq(program,C1)
program=sosineq(program,C2)
program=sosmatrixineq(program,C3)

#Objective 
program=sossetobj(program,theta)

# =============================================
# And call solver
sol = sossolve(program)

# =============================================
# Finally, get solution

f=sosgetsol(sol,f)