#SOSMATRIXINEQ --- Creates a SOS constraint from a matrix inequality constraint
# [SOSP] = sosmatrixineq(SOSP,fM)
# SOSP is the sum of squares program.
# fM is a polynomial matrix used to generate the polynomial inequality y'fMy>0

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
import sympy as sym
from sympy import *
from ineq import ineq
from polystr import polystr
from getequationmc import getequationmc
from sossolve2 import sossolve2

def matrixineq(program,p,option=[]):

    m,n=p.shape

    if m!=n:
        print('Error: Matrix fM in inequality fM>0 must be square.')
    else: 
        if option==[]:
            option='quadraticMineq' # sets the default: asks for sos expression v'M(x)v

        if option=='quadraticMineq':
            
            #Create the vector of variables Mvar to generate the quadratic expression M_var'*fM*Mvar

            var_holder = {}
            vartable=[]
            varMconst=zeros(n,1)
            for i in range(n):
                var_holder['Mvar_' + str(i)]=sym.Symbol('Mvar_%d'%(i+1))
                varMconst[i,0]=var_holder['Mvar_' + str(i)]
                vartable=np.append(vartable,var_holder['Mvar_' + str(i)])
            locals().update(var_holder)

            expr=expand((varMconst.T*p*varMconst)[0,0])


            vartabletemp=program.vartable
            program.nvars=len(vartable)
            program.vartable=np.concatenate((program.vartable,vartable),axis=None)
            program.var=np.concatenate((program.var,vartable),axis=None)

            program=ineq(program,expr)
            sos=program.sos
            sos.exprtype='sparsemultipartite'
            sos.exprmultipart=Matrix([[vartabletemp,varMconst]])
            sos.symvartable=vartabletemp
            sos.varmatsymvartable=varMconst

        if option=='Mineq':
            program.ktemp+=1
            At,b,Z,sos=getequationmc(p,program.decvartable,program.vartable)
            sos.polynomial=p
            sos.type='ineq'
            index=[]

            if program.nconstraint>1:#Multiconstraints
                sos.k=program.k
                sos.extravarnum2=program.extravarnum
                sos.extravaridx2=program.extravaridx
                sos.newtonpolytope=program.newtonpolytope
                At,b,c=sossolve2(sos)

                if program.k==0:
                    program.Ks=[int(np.sqrt(At.shape[0]))]
                else:
                    program.Ks=np.concatenate((program.Ks, int(np.sqrt(At.shape[0]))), axis=None)

                if len(index)!=0:
                    for l in range(len(index)):
                        At=np.insert(At, int(index[l]), 0, axis=0)

                if program.k==0:

                    program.cons='ineq'
                    program.polys=p
                    program.Atcon=At
                    program.bcon=b
                    program.k=program.k+1
                else:
                    program.cons=np.concatenate((program.cons, 'ineq'), axis=None)
                    program.polys=np.concatenate((program.polys, p), axis=None)

                    #Constructing the matrix At for the multiconstraints case
                    
                    nrows1=abs(program.Atcon.shape[0]-len(program.decvartable))
                    nrows2=abs(At.shape[0]-len(program.decvartable))
                    for m in range(int(nrows1)):
                         At=np.insert(At, len(program.decvartable)+m, 0, axis=0)

                    Arow=np.zeros((nrows2,program.Atcon.shape[1]))
                    program.Atcon=np.append(program.Atcon, Arow, axis=0)

                    program.Atcon=np.concatenate((program.Atcon, At), axis=1)
                    program.bcon=np.concatenate((program.bcon, b), axis=None)
                    program.k=program.k+1

                sos.Atcon=program.Atcon
                sos.bcon=program.bcon
                sos.ccon=np.zeros((sos.Atcon.shape[0],1))
                sos.Ks=program.Ks
                sos.polys=program.polys
                sos.cons=program.cons
                program.extravarnum=sos.extravarnum2
                program.extravaridx=sos.extravaridx2
            else:
                program.k=program.k+1
            sos.nconstraint=program.nconstraint
            sos.decvartable2=program.decvartable
  
    sos.newtonpolytope=program.newtonpolytope    
    program.sos=sos
    
    return program