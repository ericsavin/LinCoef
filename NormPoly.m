function [norm] = NormPoly(n,alpha,beta,TYPE)
% Compute the normalization constant of classical orthogonal polynomials
% by the usual formulas
%
% INPUT ---
% n           : order of the polynomial
% alpha, beta : parameters of the probability law
% TYPE        : type of polynomials
%               [JAC: Jacobi, GEG: Gegenbauer, CHE: Chebyshev, LAG: Laguerre, HER: Hermite]
%
% OUTPUT --
% norm        : normalization constant
%
% ESavin / April 2016 / eric.savin@onera.fr
%
% Licensing: This code is distributed under the CeCILL-C license.
%

n1 = n+1;

switch TYPE
    
    case 'JAC' % Jacobi polynomials
        apbp1 = alpha+beta+1.0;
        apbpn = alpha+beta+n;
        norm  = power(2,apbp1).*gamma(alpha+n1).*gamma(beta+n1)./factorial(n)./gamma(apbpn+1)./(apbpn+n1);
        
    case 'GEG' % Gegenbauer polynomials
        apbp1 = alpha+alpha+1.0;
        apbpn = alpha+alpha+n;
        norm  = power(2,apbp1).*gamma(alpha+n1).*gamma(alpha+n1)./factorial(n)./gamma(apbpn+1)./(apbpn+n1);
        
    case 'CHE' % Chebyshev polynomials
        if (n==0)
           norm = pi;
        else
           norm = pi./2.*power(pochhammer(alpha+1,n)./factorial(n),2);
        end
    
    case 'LAG' % Generalized Laguerre polynomials
        norm  = gamma(n1+alpha)./factorial(n);
        
    case 'HER' % Generalized Hermite polynomials
        norm  = power(2,2*n).*gamma(floor(n/2)+1).*gamma(floor(n1/2)+alpha+0.5);
        
    otherwise
        error('NormPoly: this option is not supported');
        
end

return;

