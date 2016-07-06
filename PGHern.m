function [P,Q] = PGHern(x,NN,alpha)
% Evaluate generalized Hermite polynomials of type alpha > -0.5 and their
% derivatives at points x up to the order NN
%
% Note: the computation of Q=dP/dx uses Chihara (1978) Eq.(2.47) p.157
%
% ESavin / April 2016 / eric.savin@onera.fr
%

nx = length(x);
N1 = NN+1;
P(1,1:nx) = HermiteP(x,alpha,0);
Q(1,1:nx) = zeros(1,nx);
if (NN > 0)
   gamma1    = NormPoly(1,alpha,alpha,'HER');
   P(2,1:nx) = HermiteP(x,alpha,1);
   Q(2,1:nx) = 2.*ones(1,nx)./sqrt(gamma1);
   for iN = 3:N1
      P(iN,1:nx) =  HermiteP(x,alpha,iN-1);
      Q(iN,1:nx) =  2.*(iN-1).*HermiteP(x,alpha,iN-1) + 2.*(iN-2).*(1+power(-1,iN)).*alpha.*HermiteP(x,alpha,iN-2)./x;
   end
end

return;
