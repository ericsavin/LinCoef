# LinCoef
Compute the linearization coefficients of classical orthogonal polynomials

## PRESENTATION

The main function `LinCoef.m` computes the linearization coefficients for products of classical orthogonal polynomials
of the Jacobi (including Gegenbauer, Legendre and Chebyshev), generalized Hermite, and generalized Laguerre families.
These linearization coefficients correspond to the third-order moments of these polynomial sets.
Codes are also provided to compute Gauss quadrature sets and evaluate those polynomials and their derivatives at specific
nodes. The main function to construct quadrature sets is `GNodeWt.m`.

These codes have been developed at [Onera]{http://www.onera.fr} (Computational Fluid Dynamics Dept.). They are developed
for [MATLAB](http://www.mathworks.com/products/matlab/). The codes have not been extensively tested, but should run on Matlab 
versions R2013a and newer.

* contact : [Eric Savin](mailto:eric.savin@onera.fr)
* contributors (by order of first commit): E. Savin

## REFERENCE

B. Faverjon, E. Savin. Computation of higher-order moments of generalized polynomial chaos expansions. Preprint, 2016.

## USE

The main function for computing linearization coefficients is `LinCoef.m`.
The main function for computing Gauss quadrature sets is `GNodeWt.m`.
`TestHermite.m` allows to compare the linearization formulas with Gauss-Hermite summation for generalized Hermite polynomials.
`TestJacobi.m` allows to compare the linearization formulas with Gauss-Jacobi summation for Jacobi, Gegenbauer (including
Legendre) and Chebyshev polynomials.
`TestLaguerre.m` allows to compare the linearization formulas with Gauss-Laguerre summation for generalized Laguerre
polynomials.
