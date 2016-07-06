function [P,Q] = PGLagn(x,NN,alpha)
% Evaluate generalized Laguerre polynomials of type alpha > -1 and their
% derivatives at points x up to the order NN
%
% Note: the computation of Q=dP/dx uses L_{n}^{alpha}(x)'=-L_{n-1}^{alpha+1}(x)
% after Chihara (1978) Eq. (2.22) p. 149
%
% ESavin / April 2016 / eric.savin@onera.fr
%

nx = length(x);
N1 = NN+1;
P(1,1:nx) = LaguerreP(x,alpha,0);
Q(1,1:nx) = zeros(1,nx);
for iN = 2:N1
   P(iN,1:nx) =  LaguerreP(x,alpha,iN-1);
   Q(iN,1:nx) = -LaguerreP(x,alpha+1,iN-2);
end

return;
