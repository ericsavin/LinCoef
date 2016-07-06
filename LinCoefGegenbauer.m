function [C] = LinCoefGegenbauer(j,k,alpha)
% Compute the linearization coefficients for Gegenbauer (Jacobi) polynomials
% After Chaggara-Koepf (2010) Eq. (28)
%
% INPUT ---
% j,k         : orders of the polynomials
% alpha       : parameter of the probability law
%
% OUTPUT --
% C           : non-normalized linearization coefficients for 0 =< n =< j+k
%               NOTE: C(n) = 0 for n < abs(j-k)
%
% ESavin / April 2016 / eric.savin@onera.fr
%
% Licensing: This code is distributed under the CeCILL-C license.
%

alphap = alpha + 1.0;
alpha1 = alpha + 0.5;
alpha2 = 2.*alpha1;
nmax   = min(j,k);
jpk    = j+k;

C      = zeros(jpk+1,1);

for n = 0:nmax
   jpkmn       = jpk - n;
   jpkm2n      = jpk - 2.*n;
   num1        = pochhammer(alpha2,jpkmn).*pochhammer(alpha1,j-n).*pochhammer(alpha1,k-n).*pochhammer(alpha1,n);
   den1        = pochhammer(alpha1,jpkmn).*pochhammer(alphap,jpkm2n);
   coef1       = (jpkm2n + alpha1).*factorial(jpkm2n)./(alpha1+jpkmn)./factorial(n)./factorial(j-n)./factorial(k-n);
   coef2       = pochhammer(alphap,j).*pochhammer(alphap,k)./pochhammer(alpha2,j)./pochhammer(alpha2,k);
   C(jpkm2n+1) = coef1.*coef2.*num1./den1;
end

return;

