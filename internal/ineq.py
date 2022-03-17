# INEQ --- Process a new inequality constraint f(x) >= 0
#    to an SOS program 
#
# SOSP = ineq(SOSP,EXPR,INTERVAL)
#
# SOSP is the sum of squares program.
# EXPR is the expression on the left hand side of the constraint, i.e., f(x).
#
# EXPR can be a column vector. In this case, several inequality
# constraints will be added simultaneously to the sum of
# squares program.
#
# SOSP = sosineq(SOSP,EXPR,[a b]) is used to specify the interval 
# a <= x <= b where f(x) has to be non-negative. Currently, this is 
# possible only if the sum of squares program is univariate. The lower
# limit a can be -Inf, and similarly b can be +Inf.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



import numpy as np
import symengine
import sympy as sym
from sympy import *
from sympy import degree_list,degree,matrix2numpy,sympify,poly
from sympy import ImmutableMatrix as  Matrix
from polystr import polystr
from getequation import getequation
from sossolve2 import sossolve2
from preprocess import preprocess
from expand_simplify_sos import expand_simplify_sos

def ineq(program,p,interval=[]):
    
    
    
    p=expand_simplify_sos(p,0)
    
    
    program.ktemp+=1
    
    #Analysing whether all the decision vars are using in this problem
    
    decvartemp=[]
    symb=list(p.free_symbols)
    index2=np.isin(program.decvartable,symb)
    index=[]
    for i in range(len(index2)):
        if index2[i]==True:
            decvartemp=np.concatenate((decvartemp,program.decvartable[i]), axis=None) 
        else:
            index=np.concatenate((index,i), axis=None) 
    program.decvar=decvartemp



    #Analysing whether if all the variables are using in this problem
    
    vartemp=[]
    index2=np.isin(program.vartable,symb)
    for i in range(len(index2)):
        if index2[i]==True:
            vartemp=np.concatenate((vartemp,program.vartable[i]), axis=None) 

    program.var=vartemp


    program.nvars=len(program.var)

    flag_interval=[]
    if program.nvars==1 and len(interval)!=0:
        flag_interval='interval'
        program.interval+=1
        
        p=expand_simplify_sos(str(simplify(preprocess(p,program.var,interval))),1)
        
       

    if program.nvars!=1 and len(interval)!=0:
        print("Error: Not a univariate polynomial.")



    if program.nvars==1: #Single variable

        symexpr, decvartable, varmat, vartable=polystr(p)

        if len(program.decvar)!=0: #For this case, the program is composed of decision vars
            
            decvartable=program.decvar
            vartable=program.var
            At,b,Z,sos=getequation(symexpr,decvartable,varmat,vartable)

            sos.Z=np.transpose([sos.Z])
            if flag_interval!='interval':
                sos.b=Matrix(sos.b)
                sos.b=(sos.b).T
            sos.polynomial=p
            sos.type='ineq'
            sos.interval=flag_interval
            sos.flag_interval=program.interval
        else:
            At,b,Z,sos=getequation(symexpr,decvartable,varmat,vartable)
            sos.polynomial=p
            sos.type='ineq'
            sos.interval=flag_interval
            sos.flag_interval=program.interval

    else: #Multivariable

        symexpr, decvartable, varmat, vartable=polystr(p)
        if len(program.decvar)!=0: #In this case, the program has decision vars
            decvartable=program.decvar
            vartable=program.var
            At,b,Z,sos=getequation(symexpr,decvartable,varmat,vartable)
            sos.polynomial=p
            sos.type='ineq'
            sos.interval=flag_interval
            sos.flag_interval=program.interval
        else:
            decvartable=program.decvar
            vartable=program.var
            At,b,Z,sos=getequation(symexpr,decvartable,varmat,vartable)
            sos.polynomial=p
            sos.type='ineq'
            sos.interval=flag_interval
            sos.flag_interval=program.interval


    if program.nconstraint>1:#Multiconstraints
        sos.k=program.k
        sos.extravarnum2=program.extravarnum
        sos.extravaridx2=program.extravaridx
        sos.newtonpolytope=program.newtonpolytope
        At,b,c=sossolve2(sos)


        if program.k==0:
            if program.nvars==1:
                program.Ks=[int(np.sqrt(At.shape[0]))]
            else:
                program.Ks=[int(np.sqrt(sos.extravaridx[1]-sos.extravaridx[0]))]
        else:
            if program.nvars==1:
                program.Ks=np.concatenate((program.Ks, int(np.sqrt(At.shape[0]))), axis=None)
            else:
                program.Ks=np.concatenate((program.Ks,int(np.sqrt(sos.extravaridx[1]-sos.extravaridx[0]))), axis=None)

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
        sos.flag_interval=program.interval


    elif len(program.decvar)==0:
        sos.flag_findsos=1
    sos.nconstraint=program.nconstraint
    sos.decvartable2=program.decvartable

    if program.nconstraint==1:
        program.k=program.k+1

       
    sos.newtonpolytope=program.newtonpolytope
    program.sos=sos 
   
    return program