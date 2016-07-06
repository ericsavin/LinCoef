function [C] = LinCoefLaguerre(j,k,alpha,varargin)
% Compute the linearization coefficients for Laguerre polynomials
% method = 'KDF' : after Chaggara (2007) Eq. (3.24)
% method = 'HYG' : after  Watson (1938). DEFAULT choice
%
% INPUT ---
% j,k         : orders of the polynomials
% alpha       : parameter of the probability law
% varargin    : method of computation
%
% OUTPUT --
% C           : non-normalized linearization coefficients for 0 =< n =< j+k.
%               NOTE: C(n) = 0 for n < abs(j-k)
%
% ESavin / April 2016 / eric.savin@onera.fr
%
% Licensing: This code is distributed under the CeCILL-C license.
%

method = 'HYG'; % DEFAULT choice
if (nargin >=4)
   method = varargin{1};
end

jpk    = j+k;
jmk    = abs(j-k);
nmax   = jpk-jmk;
nmin   = min(j,k);

C      = zeros(jpk+1,1);

switch method
   
   case 'KDF' % Terminating Kampe de Feriet 2:0(F)1:2 series of Chaggara (2007) Eq. (3.24)
       
      coefjk = nchoosek(jpk,k);
      bbb    = [-alpha-j; -j];
      ccc    = [-alpha-k; -k];
      ddd    = [-alpha-jpk; -jpk];
   
      for n = 0:nmax
         jpkmn      = jpk - n;
         coefn      = pochhammer(-alpha-jpk,n)./factorial(n);
         C(jpkmn+1) = coefjk.*coefn.*KdFLaguerre(n,bbb,ccc,ddd);
      end
      
   case 'HYG' % Terminating hypergeometric 3(F)2 formula of Watson (1938)
       
      for n = jmk:jpk
         jpkmn      = jpk - n;
         coef       = power(-2.0,jpkmn).*factorial(n)./factorial(jpkmn);
         C(n+1)     = HyGLaguerre(n,j,k,alpha+1).*coef;
      end
      
   otherwise
        
      error('LinCoefLaguerre: this option is not supported');
        
end

return;

