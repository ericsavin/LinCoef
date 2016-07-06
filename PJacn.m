function [P,Q] = PJacn(x,NN,alpha,beta,filter)
% Evaluate Jacobi polynomials of type (alpha,beta) > -1 and their
% derivatives at points x up to the order NN
% Apply a filter to compute filtered Jacobi polynomials
%
% Note: the computation of Q=dP/dx uses Hesthaven, Warburton: 
%       Nodal Discontinuous Galerkin Methods,
%       Springer, New York (2008), Eq.(A.2) p.445.
% See also Chihara (1978) Eq. (2.23) p. 149
%
% ESavin / April 2016 / eric.savin@onera.fr
%

nx = length(x);
N1 = NN+1;
P(1,1:nx) = filter(1).*JacobiP(x,alpha,beta,0);
Q(1,1:nx) = zeros(1,nx);
for iN = 2:N1
%%%   norm       = sqrt((iN+alpha+beta)*(iN-1)); % Normalized polynomials
   norm       = (iN+alpha+beta)./2; % Non-normalized polynomials
   P(iN,1:nx) = filter(iN).*JacobiP(x,alpha,beta,iN-1);
   Q(iN,1:nx) = filter(iN).*norm.*JacobiP(x,alpha+1,beta+1,iN-2);
end
return;
