#SOSDEMO9 --- Matrix SOS decomposition

import sympy as sym
import numpy as np
from sympy import *
from findsos import findsos
from sympy.physics.quantum import TensorProduct


#Variables
x1,x2,x3=sym.symbols('x1 x2 x3')
#=============================================
#Consider the following candidate sum of squares matrix P(x)

P=Matrix([[x1**4+x1**2*x2**2+x1**2*x3**2, x1*x2*x3**2-x1**3*x2-x1*x2*(x2**2+2*x3**2)],
          [x1*x2*x3**2-x1**3*x2-x1*x2*(x2**2+2*x3**2), x1**2*x2**2+x2**2*x3**2+(x2**2+2*x3**2)**2]])


#Test if P(x1,x2,x3) is an SOS matrix and return H so that P = H.'*H

Q,Z,H=findsos(P)

#Verify that P - H'*H = 0 and P- kron(I,Z)'*Q*kron(I,Z)= 0 to within  numerical tolerance.

Pd_1=expand(TensorProduct(eye(2),Z).T*Q*TensorProduct(eye(2),Z))
Pd_2=expand(H.T*H)