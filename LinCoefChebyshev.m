function [C] = LinCoefChebyshev(j,k,alpha)
% Compute the linearization coefficients for Chebyshev polynomials
%
% INPUT ---
% j,k         : orders of the polynomials
% alpha       : = -0.5
%
% OUTPUT --
% C           : non-normalized linearization coefficients
%
% ESavin / June 2016 / eric.savin@onera.fr
%
% Licensing: This code is distributed under the CeCILL-C license.
%

alpha1   = alpha+1;
jmk      = abs(j-k);
jpk      = j+k;
norm     = pochhammer(alpha1,j).*pochhammer(alpha1,k)./factorial(j)./factorial(k)./2;

[C]      = zeros(jpk+1,1);
C(jmk+1) = norm.*factorial(jmk)./pochhammer(alpha1,jmk);
C(jpk+1) = norm.*factorial(jpk)./pochhammer(alpha1,jpk);

return;

