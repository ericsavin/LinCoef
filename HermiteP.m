function [P] = HermiteP(x,alpha,N)
% Evaluate generalized Hermite polynomial of type (alpha) > -0.5
% at points x for order N using recurrence formula
% After Chihara 1978, pp.156-157.
%
% Note: they are normalized to be orthonormal
%
% ESavin / April 2016 / eric.savin@onera.fr
%

% Turn points into row if needed
xp = x; dims = size(xp);
if (dims(2) == 1)
   xp = xp';
end

PL = zeros(N+1,length(xp));

gamma0  = NormPoly(0,alpha,alpha,'HER');
PL(1,:) = 1.0;
if (N==0)
   P = PL'./sqrt(gamma0);
   return;
end

gamma1  = NormPoly(1,alpha,alpha,'HER');
PL(2,:) = 2.*xp;
if (N==1)
   P = PL(N+1,:)'./sqrt(gamma1);
   return;
end

% Forward recurrence
gammaN  = NormPoly(N,alpha,alpha,'HER');
for i = 1:N-1
   coef = (1-power(-1,i)).*alpha;
   PL(i+2,:) = 2.*xp.*PL(i+1,:) - 2.*(i+coef).*PL(i,:);
end

P = PL(N+1,:)'./sqrt(gammaN);

return;
