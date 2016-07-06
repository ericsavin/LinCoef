function [C] = LinCoefJacobi(j,k,alpha,beta,varargin)
% Compute the linearization coefficients for Jacobi polynomials
% method = 'KDF' : after Chaggara-Koepf (2010) Eq. (12)
% method = 'HYG' : after Abd-Elhaamed (2015) Eq. (10) with standardization of Eq. (3)
% method = 'REC' : after Gasper (1970) & Hylleraas (1962). DEFAULT choice
%
% INPUT ---
% j,k         : orders of the polynomials
% alpha, beta : parameters of the probability law
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

method = 'REC'; % DEFAULT choice
if (nargin >=5)
   method = varargin{1};
end

apb    = alpha + beta;
apbp1  = apb + 1.0;
alpha1 = alpha + 1;
beta1  = beta + 1;
jmk    = abs(j-k);
jpk    = j+k;
nmax   = jpk-jmk;
jmin   = min(j,k);
jmax   = max(j,k);

switch method
   
   case 'KDF' % Kampe de Feriet 2:0(F)1:2 formula of Chaggara-Koepf (2010) Eq. (12)
       
      C      = zeros(jpk+1,1);

      coefjk = pochhammer(alpha1,jpk).*factorial(jpk)./factorial(j)./factorial(k);
      doefjk = pochhammer(apbp1,2*j).*pochhammer(apbp1,2*k)./pochhammer(apbp1,j)./pochhammer(apbp1,k);
      bbb    = [-j; -alpha-j];
      ccc    = [-k; -alpha-k];
      ddd    = [-jpk; -alpha-jpk];
      eee    = -2*j-apb;
      fff    = -2*k-apb;

      for n = 0:nmax
         jpkmn      = jpk - n;
         num1       = pochhammer(apbp1,jpkmn).*(2.*jpkmn + apbp1).*power(-1,n);
         den1       = factorial(n).*pochhammer(alpha1,jpkmn).*pochhammer(apbp1,jpk+jpkmn+1);
         aaa        = -apbp1-jpkmn-jpk;
         C(jpkmn+1) = coefjk.*doefjk.*num1./den1.*KdFJacobi(n,aaa,bbb,ccc,ddd,eee,fff);
      end

   case 'HYG' % Hypergeometric 4(F)3 formula of Abd-Elhaamed (2015) Eq. (10) with standardization of Eq. (3)
       
      coefjk = gamma(alpha1).*gamma(beta+jmax+1).*gamma(apbp1+2*jmin);
      doefjk = factorial(jmk).*gamma(apbp1+1+2*jmax).*gamma(apbp1+jmin).*gamma(jmk+alpha1).*gamma(jmin+alpha1);
      
      coefjk = coefjk.*pochhammer(alpha1,jmin).*pochhammer(alpha1,jmax)./factorial(jmin); % standardization
       
      for n = 0:2*jmin
         jmkpn       = jmk+n;
         coefn       = (2*jmkpn+apbp1).*gamma(jmkpn+alpha1).*gamma(2*jmk+n+apbp1);
         doefn       = factorial(n).*gamma(jmkpn+beta+1);
         
         coefn       = coefn.*factorial(jmkpn)./pochhammer(alpha1,jmkpn); % standardization
         
         C(jmkpn+1)  = 0;
         
         for p = 0:n
            numerator   = [p+2*jmk+1; -jmin; -jmin-beta];
            denominator = [jmk+1; jmk+alpha1; -2*jmin-apb];
            coef        = pochhammer(-n,p).*pochhammer(n+2*jmk+apbp1,p)./factorial(p)./pochhammer(2*k+apbp1+1,p);
            C(jmkpn+1)  = C(jmkpn+1) + HyGJacobi(p,numerator,denominator).*coef;
         end
         
         C(jmkpn+1) = C(jmkpn+1).*coefjk.*coefn./doefjk./doefn;
      
      end
       
    case 'REC' % Induction formula of Gasper (1970) & Hylleraas (1962)
       
       C        = zeros(jpk+1,1);
        
       d(1)     = 0;
       d(2)     = gamma(2*jmin+apbp1).*gamma(jmax+alpha1).*gamma(jmax+beta1).*gamma(2*jmk+apbp1+1) ...
                  ./gamma(jmin+apbp1)./factorial(jmk)./gamma(jmk+alpha1)./gamma(2*jmax+apbp1+1)./factorial(jmin);
       
       C(jmk+1) = d(2).*factorial(jmk)./gamma(jmk+beta1);
       
       nn = 2;
       for n = (jmk+1):jpk
          nn     = nn + 1;
          cnm1   = 2.*(power(jpk+apbp1,2)-power(n-2+apbp1,2)).*(power(n-2+apbp1,2)-power(jmk,2)).*(n-2+beta1) ...
                   ./(2*(n-2)+apbp1)./(2*(n-2)+apbp1+1);
          cnp1   = 2.*(power(jpk+apbp1,2)-power(n,2)).*(power(n,2)-power(jmk,2)).*(n+alpha) ...
                   ./(2*n+apbp1)./(2*n+apb);
          cn     = (power(jpk+apbp1,2)-power(n+apb,2)).*(power(n,2)-power(jmk,2))./(2*n+apb) ...
                   - (power(jpk+apbp1,2)-power(n-1+apb,2)).*(power(n-1,2)-power(jmk,2))./(2*(n-1)+apb);
          d(nn)  = (cnm1.*d(nn-2) + (beta-alpha).*cn.*d(nn-1))./cnp1;
          C(n+1) = d(nn).*power(-1,n+jpk).*factorial(n)./gamma(n+beta1);
       end
       
    otherwise
        
       error('LinCoefJacobi: this option is not supported');
        
end

return;

