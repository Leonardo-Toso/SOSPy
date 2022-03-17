#Block diagonalization function for performing the block diagonalization algorithm between the 
#polytope P and the set of candidate monomials S

# This file is part of SOSTOOLSPy - Sum of Squares Toolbox for Python language ver 0.1.

# Copyright (C) 2021  L.F.Toso, G.Valmorbida,

# Send bug reports and feedback to: sostoolspy@centralesupelec.fr

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import numpy as np
import sympy as sym
import itertools

def blockdiagonalization(P,S):
    
    
    Binpos = np.asarray(list(itertools.product([0, 1], repeat=P.shape[0])))
    
    R=[]
    for i in range(Binpos.shape[0]):
        
        if int(np.sum(Binpos[i,:].dot(P),axis=0)%2)==0:
            R=np.concatenate((R,Binpos[i,:]),axis=None)
    
    R=(np.reshape(R, (int(len(R)/Binpos.shape[1]),int(Binpos.shape[1])))).T
    
    W=(R.T).dot(S)
    
    
    return W