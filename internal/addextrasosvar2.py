# Adding slack SOS variables type II 

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
import symengine
import sympy as sym
from sympy import degree_list,degree,LT,LM,LC,sympify
from sympy import matrix2numpy
from sympy.core.compatibility import as_int
from sympy import ImmutableMatrix as  Matrix
from sympy.matrices import Matrix, eye, zeros, ones, diag, Transpose
from scipy.sparse import csr_matrix
from scipy import sparse
from max_kron import max_kron 
from collect import collect
from sympy import poly
from sympy import factor
from monomials2 import monomials2
from monpolytope import monpolytope
from  getconstraint2 import getconstraint2
from findcommonZ import findcommonZ
from inconvhull import inconvhull
from newtonpolytope import newtonpolytope
from diaginconsistency import diaginconsistency

def addextrasosvar2(sos,I):
    sos.At=np.flip(sos.At, axis=1)
    for i in range(I):
        
       
        
        numstates=sos.Z.shape[1]
        
        # Creating extra variables
        
        maxdeg=np.amax(sos.Z)
        mindeg = np.amin(sos.Z)
        
        
        floor_deg=np.floor(maxdeg/2) 
        ceil_deg=np.ceil(mindeg/2)
        d=[ceil_deg,floor_deg + 1]
        
        
        
        #Creating the candidate monomials
        
        Z=monpolytope(numstates,d)
         
            
        #Applying the Newton polytope algorithm to reduce the set of candidate monomials
        
        if sos.find==1:
            
                #Newton Polytope 
                Zhull=newtonpolytope(Z,sos.Z)
                Z=Zhull

               
        else: 
            
            if sos.newtonpolytope==1:
                
                #Newton Polytope 
                Zhull=newtonpolytope(Z,sos.Z)
                Z=Zhull

              

        #Discarding unnecessary monomials
        
        
        maxdegree = np.max(sos.Z, axis=0)/2 
        mindegree = np.min(sos.Z, axis=0)/2 
        
        j=0
        
        while j+1<=Z.shape[0]:
            Zdummy1 = maxdegree-Z[j][:]
            Zdummy2 = Z[j][:]-mindegree
            Zdummy=[Zdummy1,Zdummy2]
            idx=[]
            for i in range(len(Zdummy)):
                if Zdummy[i][0]<0:
                    idx=np.concatenate((idx,i),axis=None)
            if len(idx)!=0:
                Z = np.array([[Z[0:j,:]], [Z[j:][:]]])
            else:
                j = j+1
            
        #Adding slack variables
        for k in range(2):
            sos.extravarnum=sos.extravarnum+1
            var=sos.extravarnum
            sos.extravarZ=Z
            
           
            sos.extravarnum2=sos.extravarnum2+1
            var2=sos.extravarnum2
            
            #Getting the constraint expression


            T,ZZ=getconstraint2(Z) 
            sos.extravarZZ=ZZ
            sos.extravarT=T.T
            sos.extravaridx=np.append(sos.extravaridx,sos.extravaridx[var-1]+(Z.shape[0])**2)  
            sos.extravaridx2=np.append(sos.extravaridx2,sos.extravaridx2[var2-1]+(Z.shape[0])**2)  

            #Completing sos.At

            sos.At=np.concatenate((sos.At,np.zeros((sos.extravarT.shape[0],sos.At.shape[1]))),axis=0)

         

            #Modifying expression
            degoffset = k*np.ones((sos.extravarZZ.shape[0],1))
            
            #Ensure correct size
            class  pc: 
                F=[]
                Z=[]
            
            
            
            
            pc.Z=sos.extravarZZ+degoffset
            pc.F=-np.eye(pc.Z.shape[0], dtype=int)
            
            R1,R2,newZ=findcommonZ(sos.Z,pc.Z)
            
            if len(sos.At)==0:
                sos.At=np.zeros((sos.At.shape[0],R1.shape[0]))
            
            sos.At=np.matmul(sos.At,R1)
            lidx=sos.extravaridx[var-1]
            uidx=sos.extravaridx[var]-1
            sos.At[lidx-1:uidx][:]=sos.At[lidx-1:uidx][:]-np.matmul(np.matmul(sos.extravarT,pc.F),R2)
            sos.b=np.matmul(R1.T,sos.b)
            sos.Z=newZ
            
            findidx=[]
            for j in range(Z.shape[0]):
                if Z[j]<maxdegree:
                    findidx=np.concatenate((findidx,j),axis=None)
            Z=Z[findidx.astype(int)]       
            
    return sos