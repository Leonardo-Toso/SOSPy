#POLYSTR--- Construct a polynomial structure composed of a variety of information from the polynomial structure

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
from sympy import degree_list
from sympy import degree
from sympy import LT
from dataclasses import dataclass
from sympy import sympify
import time
from sympy.parsing.sympy_parser import parse_expr
from sympy.polys.orderings import monomial_key


def polystr (polynomial):
    

    #Given a polynomial we would like to extract some informations of them.
    
    #Extract the coefficients   
    
    coeff=polynomial.as_coefficients_dict()
    
 
    expr_mon=list(coeff.items())
    expr_mon=np.array(expr_mon)
    
    coeff= [ v for v in coeff.values() ]
    coeff=np.asarray(coeff)
    coeff=np.transpose(coeff)
    
    #Extract the number of terms 
    
    nterms=len(coeff)
    
    #Exctract the independent coefficient- coeff0
    
    ptemp=polynomial
    coeff0=0
    
    
    pol=str(ptemp)

    #Spliting the polynomial
    
    terms=pol.split()
    terms=list(filter(('+').__ne__, terms))
    terms_complete=[]
    j=0
    flag_terms=0
    for i in range(len(terms)):
        if terms[i]=='-':
            terms_complete=np.concatenate((terms_complete,terms[i] + terms[i+1]),axis=None)
            flag_terms=1
        else:
            if flag_terms==0:
                terms_complete=np.concatenate((terms_complete,terms[i]),axis=None)
            flag_terms=0
    terms=list(terms_complete)
    
 
    flag_coeff0=sympify(parse_expr(terms[len(terms)-1])).is_number 
    ptemp=parse_expr(terms[len(terms)-1])
    
    
    
    
    if ptemp!=0 and flag_coeff0==True:
        flag_coeff0=1
    else:
        flag_coeff0=0
    
    if flag_coeff0==1:
        coeff0=ptemp
        polynomial=polynomial-ptemp
        
    #Extract the coefficients   
    
    coeff=polynomial.as_coefficients_dict()
    coeff= [ v for v in coeff.values() ]
    coeff=np.asarray(coeff)
    coeff=np.transpose(coeff)
   
    
    #Extract the number of terms 
    
    nterms=len(coeff)
    
    
    #Extract the variables
    
    varname=polynomial.free_symbols
    varname=list(varname)
    varname=np.asarray(varname)
    
    
    #Extract the number of variables
    
    nvars=len(varname)
    
    cvartable=np.array(varname).astype('str').tolist()
    cvartable=np.asarray(cvartable)
    cvartable=np.sort(cvartable)
    
    # Create the matrix of polynomial degree
    
    degmat=np.zeros((nterms,nvars))
    ptemp=polynomial
    
    i=0
    k=0
    while(i<=nterms-1 and k<=nvars-1):
        monomial=LT(ptemp,gens=sym.Symbol(cvartable[k]))
        if i<nterms-1:
            ptemp=ptemp-monomial
        if ptemp!=0:
            
            for j in range(nvars):
                degree_var=degree(monomial,gen=sym.Symbol(cvartable[j]))
                if degree==' ' :
                    degmat=[i][j]=0
                else: 
                    degmat[i][j]=degree_var
                    
            i=i+1
        if ptemp==0:
            if k<=nvars-1:
                k=k+1
                ptemp=monomial 
            else:
                ptemp=monomial 
    
    
    if coeff0!=0:
        
        extra_row=np.zeros((1,nvars))
        degmat=np.concatenate((degmat,extra_row), axis=0)
        
      
    #Create maxdeg and mindeg
    
    max_min=np.zeros((1,nterms))
    for i in range(nterms):
        for j in range(nvars):
            max_min[0][i]+=degmat[i][j]
            
    max_min=np.sort(max_min)
    maxdeg=max_min[0][nterms-1]
    mindeg=max_min[0][0]
    
    #Create decvartable
    
    decvartable=[]  
    
    #Create varmat
    
    varmat=[]
    
    
    vartable=[]
    for i in range(len(cvartable)):
        vartable=np.append(vartable,sym.Symbol(cvartable[i]))
    
    vartable=list(vartable)
    
    #Create the structure symexpr
    
    class symexpr:
        polynomial=0
        coefficient=0
        degmat=0 
        varname=0
        nterms=0
        nvars=0
        maxdeg=0 
        mindeg=0 
        coeff0=0
        
    #Define the valeus in symexpr  
    
    symexpr=symexpr() 
    symexpr.polynomial=polynomial
    symexpr.coefficient=coeff
    symexpr.nterms=nterms
    symexpr.varname=varname
    symexpr.nvars=nvars
    symexpr.degmat=degmat
    symexpr.maxdeg=maxdeg
    symexpr.mindeg=mindeg
    symexpr.coeff0=coeff0
    
    
       
    return symexpr, decvartable, varmat, vartable