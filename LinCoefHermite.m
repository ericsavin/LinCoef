function [C] = LinCoefHermite(j,k,alpha)
% Compute the linearization coefficients for generalized Hermite polynomials
% After Chaggara 2010, Chaggara-Koepf 2011 Eq. (3.5) with the
% standardization of Chihara 1978.
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

jf   = floor(j/2);
kf   = floor(k/2);
nmax = min(j,k);
jpk  = j+k;

% Rosenblum (1994) standardization
%   Gjk = factorial(j).*factorial(k);
% Chihara (1978) standardization
Gjk  = GaMu(j,alpha).*GaMu(k,alpha);

C    = zeros(jpk+1,1);

for n = 0:nmax
   jpkm2n = jpk - 2.*n;
% Rosenblum (1994) standardization
%    coef   = Gjk./factorial(jpkm2n)./factorial(n);
% Chihara (1978) standardization
   coef   = Gjk./GaMu(jpkm2n,alpha)./factorial(n);
   for ip = 0:jf
      coefj = GaMu(j-2.*ip,alpha).*factorial(ip);
      for iq = 0:kf
         coefk = GaMu(k-2.*iq,alpha).*factorial(iq);
         C(jpkm2n+1) = C(jpkm2n+1) + coef.*GaMu(jpk - 2.*(ip+iq),alpha).*pochhammer(-n,ip+iq)./coefj./coefk;  
      end
   end
end

return;

