# SOSPy
A free Python toolbox that allows Python users to formulate and solve sums of squares (SOS) optimization programs.

## Introduction
**SOSPy** is a free Python toolbox that allows Python users to formulate and solve sums of squares (SOS) optimization programs. SOSPy can be used to specify and solve sum of squares polynomial problems using a very simple, flexible, and intuitive high-level notation. The SOS programs can be solved using [Mosek](https://github.com/MOSEK). This well-known semidefinite programming solver, with SOSPy handling internally all the necessary reformulations and data conversion.

## What is a "sum of squares optimization program"? Why would I want such a thing?

A *sum of squares (SOS) program*, in the simplest case, has the form:
```
minimize: c_1 * u_1 + ... + c_n * u_n

subject to constraints: 

P_i(x) := A_i0(x) + A_i1(x) * u_1 + ... + A_in(x) * u_n

are sums of squares of polynomials (for i=1..n).
```

Here, the `A_ij(x)` are multivariate polynomials, and the decision variables `u_i` are scalars. This is a convex optimization problem, since the objective function is linear and the set of feasible `u_i` is convex.

While this looks quite nice, perhaps you are actually interested in more concrete problems such as:

* Constrained or unconstrained optimization of polynomial functions.
* Mixed continuous-discrete optimization.
* Finding Lyapunov or Bendixson-Dulac functions for nonlinear dynamical systems (with polynomial vector fields).
* Deciding copositivity of a matrix.
* Inequalities in probability theory.
* Distinguishing separable from entangled states in quantum systems.

Or, more generally, problems that deal with basic semialgebraic sets (sets defined by polynomial equalities and inequalities).


## Distribution and release information

SOSPy is freely available under the GNU public license v3.0.


## System requirements

To install and run SOSTOOLS, you need:

* [Anaconda Python](https://www.anaconda.com/products/individual).
* Python Symbolic toolbox (https://www.sympy.org/pt/index.html).
* An SDP solver, [Mosek](https://github.com/MOSEK).
* SOSPy can easily be run on Windows or macOS machines.


## Authors

The software has been written and is maintained by:

* [Leonardo Felipe Toso]
* [Giorgio Valmorbida](https://www.l2s.centralesupelec.fr/perso/giorgio.valmorbida)


## References
For a detailed explanation of the theory and applications of sums of squares programming, as well as references to related work, please see:

* *Structured Semidefinite Programs and Semialgebraic Geometry Methods in Robustness and Optimization*.
California Institute of Technology, Pasadena, CA, May 2000. [Abstract](https://www.mit.edu/~parrilo/pubs/files/Thesis_abstract.html), [Pdf](http://www.mit.edu/~parrilo/pubs/files/thesis.pdf), [CaltechTHESIS](https://resolver.caltech.edu/CaltechETD:etd-05062004-055516). 

* *Semidefinite programming relaxations for semialgebraic problems*.
P. A. Parrilo, Mathematical Programming Ser. B, Vol. 96, No.2, pp. 293-320, 2003.  [Abstract](https://www.mit.edu/~parrilo/pubs/files/SDPrelax_abstract.html), [pdf](http://www.mit.edu/~parrilo/pubs/files/SDPrelaxations.pdf).

* *Minimizing polynomial functions*
P. A. Parrilo, B. Sturmfels. [arXiv](https://arxiv.org/abs/math.OC/0103170).

* *SOSTOOLS MATLAB
For more information and examples concerning sum-of-squares programming, please see [Getting started with Sum of Squares](https://sums-of-squares.github.io/sos/index.html#matlab) and the [SOSTOOLS user's guide](docs/sostools.pdf).


 
## Feedback
For comments, bug reports, encouragement, suggestions, complaints, etc., please send email to: leonardo-felipe.toso@student-cs.fr.

If you use SOSPy for research purposes, we would be happy to hear about it and mention it in the reference guide. Please drop us a line, to leonardo-felipe.toso@student-cs.fr.

## Citing 

Please use the following when citing SOSPy:

```
@manual{sostools,
author = {L. F. Toso and G. Valmorbida},
title = {{SOSPy}: Sum of squares optimization toolbox for {Python}},
note = {},
year = {2022},
address = {},
}
```
