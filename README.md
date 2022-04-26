# LinCoef
Compute the linearization coefficients of classical orthogonal polynomials

## PRESENTATION

The main function `LinCoef.m` computes the linearization coefficients for products of classical orthogonal polynomials
of the __Jacobi__ (including Gegenbauer, Legendre and Chebyshev), __generalized Hermite__, and __generalized Laguerre__ families.
These linearization coefficients correspond to the third-order moments of polynomial sets.
Codes are also provided to compute Gauss quadrature sets and evaluate those polynomials and their derivatives at specific
nodes. The main function to construct quadrature sets is `GNodeWt.m`.

These codes have been developed at [ONERA](http://www.onera.fr) ([Computational Fluid Dynamics Dept.](https://www.onera.fr/en/daaa)) for [Matlab](http://www.mathworks.com/products/matlab/). They have not been extensively tested, but should run on Matlab 
versions R2013a and newer.

* contact : [Éric Savin](mailto:eric.savin@onera.fr)
* contributors (by order of first commit): É. Savin

## REFERENCE

[É. Savin, B. Faverjon. Higher-order moments of generalized polynomial chaos expansions for intrusive and non-intrusive uncertainty quantification. AIAA Paper \#2017-0597 (2017);](https://doi.org/10.2514/6.2017-0597) [arXiv:1607.01914.](https://arxiv.org/abs/1607.01914)

[É. Savin, B. Faverjon. Computation of higher-order moments of generalized polynomial chaos expansions. *Int. J. Num. Methods Engng.* __111__(12), 1192-1200 (2017).](https://doi.org/10.1002/nme.5505)

## USE

The main function for computing linearization coefficients is `LinCoef.m`.

The main function for computing Gauss quadrature sets is `GNodeWt.m`.

`TestHermite.m` allows to compare the linearization formulas with Gauss-Hermite summation for generalized Hermite polynomials.

`TestJacobi.m` allows to compare the linearization formulas with Gauss-Jacobi summation for Jacobi, Gegenbauer (including
Legendre) and Chebyshev polynomials.

`TestLaguerre.m` allows to compare the linearization formulas with Gauss-Laguerre summation for generalized Laguerre
polynomials.
