# GETEQUATIONMC --- Convert a symbolic matrix to At, b, and Z  used in an SOS program.

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import numpy as np
from typing import NamedTuple
import symengine
import sympy as sym
from sympy import *
from sympy import LT,LM,degree,degree_list,expand,degree
from collect import collect 
from dataclasses import dataclass
from sympy.matrices import Matrix, eye, zeros, ones, diag, Transpose
from scipy.sparse import csr_matrix
from ismember import ismember
from sortNoRepeat import sortNoRepeat
from max_kron import max_kron
from polystr import polystr
from getequation  import getequation
from getequation2  import getequation2
from sympy import sympify
import sys

def getequationmc(symexpr,decvartable,vartable):
    
    
    if len(decvartable)==0:
        decvarnum=0
        
    else:
        decvarnum=len(decvartable)
        
    
    FPexpr=symexpr
    
    dimp=FPexpr.shape[1]
    
    if decvarnum==0:
        
        Zfull=[]
        k=0
        for i in range(dimp):
            for j in range(dimp):
                if sympify(FPexpr[i,j]).is_number==False:
                    symexpr,decvartablei,varmat,vartablei=polystr(FPexpr[i,j])
                    b,Z=getequation2(symexpr,decvartable,varmat,vartable)
                elif sympify(FPexpr[i,j]).is_number==True and FPexpr[i,j]!=0 :
                    Z=np.zeros((len(vartable),1))
                    b=np.zeros((1,len(vartable)))
                    b[0][b.shape[1]-1]=FPexpr[i,j]
                    b=Matrix(b)
                
                if k==0:
                    Zfull=Z
                    k=1
                else:
                    Zfull=np.concatenate((Zfull,Z),axis=0)
    
        Zfull=np.flip(np.flip(np.unique(Zfull,axis=0)),axis=1)
        
        Z=Zfull
        
        nmon,nvar=Zfull.shape
        
        coeffnts=np.zeros((nmon*dimp,dimp))
        
        for i in range(dimp):
            for j in range(dimp):
                if sympify(FPexpr[i,j]).is_number==False:
                    symexpr,decvartablei,varmat,vartablei=polystr(FPexpr[i,j])
                    b,Z=getequation2(symexpr,decvartable,varmat,vartable)
                    
                    for k in range(b.shape[1]):
                        dummyvar=[b[0,k],Z[k,:]]
                        s_ijk =dummyvar[0]
                        mon_k= dummyvar[1]

                        val_sum=np.sum((Zfull==np.kron(np.ones((nmon,1)),mon_k)).astype(int),axis=1)
                        val=np.amax(np.sum((Zfull==np.kron(np.ones((nmon,1)),mon_k)).astype(int),axis=1))
                        ind_k=[]
                        for l in range(len(val_sum)):
                            if val_sum[l]==val:
                                ind_k=np.concatenate((ind_k,l),axis=None)
                        coeffnts[(ind_k.astype(int))*dimp+i,j] = s_ijk
                        coeffnts[(ind_k.astype(int))*dimp+j,i] = s_ijk
                
                elif sympify(FPexpr[i,j]).is_number==True and FPexpr[i,j]!=0 :
                    Z=np.zeros((len(vartable),1))
                    b=np.zeros((1,Z.shape[0]))
                    b[0][b.shape[1]-1]=FPexpr[i,j]
                    b=Matrix(b)
                    
                    for k in range(b.shape[1]):
                        dummyvar=[b[0,k],Z[k,:]]
                        s_ijk =dummyvar[0]
                        mon_k= dummyvar[1]

                        val_sum=np.sum((Zfull==np.kron(np.ones((nmon,1)),mon_k)).astype(int),axis=1)
                        val=np.amax(np.sum((Zfull==np.kron(np.ones((nmon,1)),mon_k)).astype(int),axis=1))
                        ind_k=[]
                        for l in range(len(val_sum)):
                            if val_sum[l]==val:
                                ind_k=np.concatenate((ind_k,l),axis=None)
                        coeffnts[(ind_k.astype(int))*dimp+i,j] = s_ijk
                        coeffnts[(ind_k.astype(int))*dimp+j,i] = s_ijk
                
                
    else:
        Zfull=[]
        
        k=0
        for i in range(dimp):
            for j in range(dimp):
                
                
                if sympify(FPexpr[i,j]).is_number==False:
                    symexpr,decvartablei,varmat,vartablei=polystr(FPexpr[i,j])
                    b,Z=getequation2(symexpr,[],varmat,vartable)
                elif sympify(FPexpr[i,j]).is_number==True and FPexpr[i,j]!=0 :
                    Z=np.zeros((len(vartable),1))
                    b=np.zeros((1,Z.shape[0]))
                    b[0][b.shape[1]-1]=FPexpr[i,j]
                    b=Matrix(b)
                
                
                if k==0:
                    Zfull=Z
                    k=1
                else:
                    Zfull=np.concatenate((Zfull,Z),axis=0)
                
        Zfull=np.flip(np.flip(np.unique(Zfull,axis=0)),axis=1)
        
        Z=Zfull
        
        nmon,nvar=Zfull.shape
        
        coeffnts=np.zeros((nmon*dimp,dimp))
        coeffnts_decvar=zeros(nmon*dimp,dimp)
        
        for i in range(dimp):
            for j in range(dimp):
                if sympify(FPexpr[i,j]).is_number==False:
                     
                    symexpr,decvartablei,varmat,vartablei=polystr(FPexpr[i,j]) 
                    b2,Z=getequation2(symexpr,[],varmat,vartable)
                    ZLM=Z
                    
                    Z=np.flip(np.flip(np.unique(Z,axis=0)),axis=1)
                    
                    FPexprtemp=FPexpr[i,j]
                    indices=[]
                    indp=[]
                    for n in range(ZLM.shape[0]):
                        ind=[]
                        for m in range(ZLM.shape[0]):   
                            value_temp=np.asarray((ZLM[n,:]==ZLM[m,:]).astype(int))
                            value=np.sum(value_temp)
                            if value==ZLM.shape[1]:
                                ind=np.concatenate((ind,m),axis=None)
                        
                        if list(ind)!=list(indp):
                            indices=np.concatenate((indices,Matrix([[ind]])),axis=None)
                        indp=ind
                              
                            
                    
                    coeffs_temp=zeros(ZLM.shape[0],1)
                    for m in range(symexpr.nterms):
                        FPexprLM=LT(FPexprtemp)/prod(vartable**(ZLM[m,:]))
                        FPexprtemp=FPexprtemp-LT(FPexprtemp)
                        coeffs_temp[m,0]=FPexprLM
                    
                    coeffs=zeros(Z.shape[0],1)
                    
                    
                    for m in range(len(indices)):
                        val_temp=0
                        if len(indices[m])>1:
                            for n in range(len(indices[m])):
                                val_temp=+val_temp + coeffs_temp[int(indices[m][n])]
                            
                            coeffs[m]=val_temp
                        else:
                            coeffs[m]=coeffs_temp[int(indices[m][0])]
                    
                   
                    if list(Z[Z.shape[0]-1,:])==list((np.zeros((1,Z.shape[1])))[0,:]):
                        coeffs[coeffs.shape[0]-1,0]= coeffs[coeffs.shape[0]-1,0] +symexpr.coeff0 
                    
                    coeff_temp=coeffs
                    for a in range(len(decvartable)):
                        coeff_temp=coeff_temp.subs(decvartable[a],0)
                    
                   
                    
                    
                    b=coeff_temp.T
                    
                    for k in range(b.shape[1]):
                        dummyvar=[b[0,k],Z[k,:]]
                        dummy_decvar=[coeffs[k,0],Z[k,:]]
                        s_ijk_decvar = dummy_decvar
                        s_ijk =dummyvar[0]
                        mon_k= dummyvar[1]

                        val_sum=np.sum((Zfull==np.kron(np.ones((nmon,1)),mon_k)).astype(int),axis=1)
                        val=np.amax(np.sum((Zfull==np.kron(np.ones((nmon,1)),mon_k)).astype(int),axis=1))
                        ind_k=[]
                        for l in range(len(val_sum)):
                            if val_sum[l]==val:
                                ind_k=np.concatenate((ind_k,l),axis=None)
                        
                        coeffnts[(ind_k.astype(int))*dimp+i,j] = s_ijk
                        coeffnts[(ind_k.astype(int))*dimp+j,i] = s_ijk
                        
                        coeffnts_decvar[int((ind_k.astype(int))*dimp+i),j] = s_ijk_decvar[0]
                        coeffnts_decvar[int((ind_k.astype(int))*dimp+j),i] = s_ijk_decvar[0]
                
                elif sympify(FPexpr[i,j]).is_number==True:
                    Z=np.zeros((len(vartable),1))
                    nterms=0
                    ZLM=Z
                    FPexprtemp=FPexpr[i,j]
                    indices=[]
                    indp=[]
                    
                    for n in range(ZLM.shape[0]):
                        ind=[]
                        for m in range(ZLM.shape[0]):   
                            value_temp=np.asarray((ZLM[n,:]==ZLM[m,:]).astype(int))
                            value=np.sum(value_temp)
                            if value==ZLM.shape[1]:
                                ind=np.concatenate((ind,m),axis=None)
                        
                        if list(ind)!=list(indp):
                            indices=np.concatenate((indices,Matrix([[ind]])),axis=None)
                        indp=ind
                              
                            
                    
                    coeffs_temp=zeros(ZLM.shape[0],1)
                    for m in range(nterms):
                        FPexprLM=LT(FPexprtemp)/prod(vartable**(ZLM[m,:]))
                        FPexprtemp=FPexprtemp-LT(FPexprtemp)
                        coeffs_temp[m,0]=FPexprLM
                    
                    coeffs=zeros(Z.shape[0],1)
                    
                    
                    for m in range(len(indices)):
                        val_temp=0
                        if len(indices[m])>1:
                            for n in range(len(indices[m])):
                                val_temp=+val_temp + coeffs_temp[int(indices[m][n])]
                            
                            coeffs[m]=val_temp
                        else:
                            coeffs[m]=coeffs_temp[int(indices[m][0])]
                    
                   
                    if list(Z[Z.shape[0]-1,:])==list((np.zeros((1,Z.shape[1])))[0,:]):
                        coeffs[coeffs.shape[0]-1,0]= coeffs[coeffs.shape[0]-1,0] + FPexpr[i,j]
                    
                    coeff_temp=coeffs
                    for a in range(len(decvartable)):
                        coeff_temp=coeff_temp.subs(decvartable[a],0)
                    
                   
                    
                    
                    b=coeff_temp.T
                    
                    for k in range(b.shape[1]):
                        dummyvar=[b[0,k],Z[k,:]]
                        dummy_decvar=[coeffs[k,0],Z[k,:]]
                        s_ijk_decvar = dummy_decvar
                        s_ijk =dummyvar[0]
                        mon_k= dummyvar[1]

                        val_sum=np.sum((Zfull==np.kron(np.ones((nmon,1)),mon_k)).astype(int),axis=1)
                        val=np.amax(np.sum((Zfull==np.kron(np.ones((nmon,1)),mon_k)).astype(int),axis=1))
                        ind_k=[]
                        for l in range(len(val_sum)):
                            if val_sum[l]==val:
                                ind_k=np.concatenate((ind_k,l),axis=None)
                        
                        coeffnts[(ind_k.astype(int))*dimp+i,j] = s_ijk
                        coeffnts[(ind_k.astype(int))*dimp+j,i] = s_ijk
                        
                        coeffnts_decvar[int((ind_k.astype(int))*dimp+i),j] = s_ijk_decvar[0]
                        coeffnts_decvar[int((ind_k.astype(int))*dimp+j),i] = s_ijk_decvar[0]
        
    
    #Constructing At and b matrices
    
    Z=Zfull

    At=np.zeros((decvarnum,Z.shape[0]*dimp**2))

    if decvarnum!=0:
        for i in range(Z.shape[0]):
            Mivec = Matrix(np.reshape(coeffnts_decvar[(i)*dimp:(i+1)*dimp,:],(dimp**2,1)))
            At[:,(i)*dimp**2:(i+1)*dimp**2]=-(Mivec.jacobian(decvartable)).T
            
    
    
    
    
    b=coeffnts                
            
        
            
    class sos:
        At=0
        b=0
        Z=0
        extravarZ=[]
        extravarZZ=[]
        extravarT=[]
        extravaridx=[]
        extravarnum=[]
        Zget=[]
        vartable=[]
        coeff0=[]
        decvartable=[]
        polynomial=[]
        soldecvar=[]
        type=[]
        Atcon=[]
        bcon=[]
        ccon=[]
        nconstraint=1
        decvartable2=[]
        Ks=[]
        polys=[]
        cons=[]
        varidx=[]
        varnum=[]
        exprnum=[]
        flag_findsos=[]
        obj=[]
        interval=[]
        k=0
        extravaridx2=[]
        extravarnum2=0
        flag_interval=0
        exprtype=[]
        exprmultipart=[]
        symvartable=[]
        varmatsymvartable=[]
        
        
    sos.At=At
    sos.b=b
    sos.Z=Z
    sos.extravaridx=[1]
    sos.Zget=Z
    sos.vartable=vartable
    sos.coeff0=symexpr.coeff0
    sos.decvartable=decvartable
    sos.varnum=-1
    sos.extravarnum=0
    sos.exprnum=1
    
   
    return At,b,Z,sos