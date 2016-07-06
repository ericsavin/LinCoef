function [Bjkn] = LinCoef(j,k,alpha,beta,TYPE)
% Compute the linearization coefficients of classical orthonormal polynomials
%
% INPUT ---
% j,k         : orders of the polynomials
% alpha, beta : parameters of the probability law
% TYPE        : type of polynomials
%               [JAC: Jacobi, GEG: Gegenbauer, CHE: Chebyshev, LAG: Laguerre, HER: Hermite]
%
% OUTPUT --
% Bjkn        : normalized linearization coefficients for 0 =< n =< j+k
%               NOTE: Bjkn(n) = 0 for n < abs(j-k)
%
% ESavin / April 2016 / eric.savin@onera.fr
%
% Licensing: This code is distributed under the CeCILL-C license.
%

switch TYPE
    
    case 'JAC' % Jacobi polynomials
        [C]  = LinCoefJacobi(j,k,alpha,beta);
        
    case 'GEG' % Gegenbauer polynomials
        if (alpha==-.5)
           [C]  = LinCoefChebyshev(j,k,alpha); 
        else
           [C]  = LinCoefGegenbauer(j,k,alpha);
        end
        beta = alpha;
        
    case 'CHE' % Chebyshev polynomials
        [C]  = LinCoefChebyshev(j,k,alpha);
        beta = alpha;
    
    case 'LAG' % Generalized Laguerre polynomials
        [C]  = LinCoefLaguerre(j,k,alpha);
        beta = alpha;
        
    case 'HER' % Generalized Hermite polynomials
        [C]  = LinCoefHermite(j,k,alpha);
        beta = alpha;
        
    otherwise
        error('LinCoef: this option is not supported');
        
end

% Normalization
jpk    = j+k;
gammaj = NormPoly(j,alpha,beta,TYPE);
gammak = NormPoly(k,alpha,beta,TYPE);
for l = 0:jpk
   gammal  = NormPoly(l,alpha,beta,TYPE);
   Bjkn(l+1) = C(l+1).*sqrt(gammal./gammaj./gammak);
end

return;

