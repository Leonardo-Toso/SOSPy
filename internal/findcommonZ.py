# FINDCOMMONZ --- Find common Z and permutation matrices R1, R2

# [R1,R2,Z] = findcommonZ(Z1,Z2)

# Given two vectors of monomials Z1 and Z2, this 
# function will compute another vector of monomials Z
# containing all the monomials of Z1 and Z2, and
# permutation matrices R1, R2 such that

#  Z1 = R1*Z
#  Z2 = R2*Z

# Assumption: all the monomials in Z1, as well as
# the monomials in Z2, are DISTINCT --- but Z1 and Z2 may 
# have common monomials.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
import sympy as sym
from sympy import degree_list,degree,LT,LM,LC,sympify
from sympy import matrix2numpy
from sympy.core.compatibility import as_int
from getequation import getequation
from sympy import ImmutableMatrix as  Matrix
from sympy.matrices import Matrix, eye, zeros, ones, diag, Transpose
from max_kron import max_kron 
from collect import collect
from sympy import poly
from sympy import factor


def findcommonZ(Z1,Z2):
    
    
    if (Z1.shape[0] + Z2.shape[0]) <=1:
        Z=np.concatenate((Z1,Z2),axis=0)
        
        #Creating the matrix R1
        
        R1=np.eye(Z1.shape[0],Z.shape[0])
        
        #Creating the matrix R2
        
        R1=np.eye(Z2.shape[0],Z.shape[0])
        
    
    else:
        sizeZ1=Z1.shape[0]
        Ind1=np.zeros((sizeZ1,1))
        for i in range(sizeZ1):
            Ind1[i][0]=i+1
        
        sizeZ2=Z2.shape[0]
        Ind2=np.zeros((sizeZ2,1))
        for i in range(sizeZ2):
            Ind2[i][0]=i+1
                
     
        Ind12=(sizeZ2+1)*np.ones((Ind1.shape[0],1))
        Ind1block=np.concatenate((Ind1,Ind12),axis=1)
        Ind21=(sizeZ1+1)*np.ones((Ind2.shape[0],1)) 
        Ind2block=np.concatenate((Ind21,Ind2),axis=1)
        Ind=np.concatenate((Ind1block,Ind2block),axis=0)
        
        
        #Constructing Z
        
        ZZ=np.concatenate((Z1,Z2),axis=0)
        IndSort=np.lexsort(np.fliplr(ZZ).T)
        ZZ=ZZ[np.lexsort(np.fliplr(ZZ).T)]
        
        sizeZZ=ZZ.shape[0]
       
        Ztemp1=np.zeros((1,ZZ.shape[1]))
        for i in range(Ztemp1.shape[1]):
            Ztemp1[0][i]=ZZ[sizeZZ-1][i]
        
        Ztemp2=np.zeros((sizeZZ-1,ZZ.shape[1]))
        for i in range(Ztemp2.shape[0]):
            for j in range(Ztemp2.shape[1]):
                Ztemp2[i][j]=ZZ[i][j]
        
        Ztemp=ZZ- np.concatenate((Ztemp1,Ztemp2),axis=0)
        
        sum_Ztemp=np.sum(abs(Ztemp),axis=1)
        
        I=[]
        
        for i in range(len(sum_Ztemp)):
            if sum_Ztemp[i]>0:
                
                I=np.concatenate((I,i),axis=None)
    
        INull=[]
        
        for i in range(len(sum_Ztemp)):
            if sum_Ztemp[i]==0:
                INull=np.concatenate((INull,i),axis=None)
        
        if len(I)==0:
            I=1
            INull=2
    
        Z=np.zeros((len(I),ZZ.shape[1]))
        for i in range(len(I)):
            Z[i][:]=ZZ[int(I[i])][:]
        
        #Constructing permutation matrix
        
        Ind=Ind[IndSort,:]
       
       
        
        for i in INull:
            Ind[int(i-1),1] = Ind[int(i),1]
            Ind[int(i),1] = sizeZ2+1;
        
        I=I.astype(int)
        Ind=Ind[I,:]
        Ind=Ind.astype(int)
        
        #Constructing the matrix R1
        
        R1temp1=np.eye(sizeZ1)
        R1temp2=np.zeros((sizeZ1,len(I)-sizeZ1))
        R1=np.concatenate((R1temp1,R1temp2),axis=1)
        R1=R1[:,Ind[:,0]-1]
        
        #Constructing the matrix R2
        
        R2temp1=np.eye(sizeZ2)
        R2temp2=np.zeros((sizeZ2,len(I)-sizeZ2))
        R2=np.concatenate((R2temp1,R2temp2),axis=1)
        R2=R2[:,Ind[:,1]-1]
        
       
        
    return R1,R2,Z