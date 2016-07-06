function [P] = LaguerreP(x,alpha,N)
% Evaluate generalized Laguerre polynomial of type (alpha) > -1
% at points x for order N using recurrence formula
%
% Note: they are normalized to be orthonormal

% ESavin / April 2016 / eric.savin@onera.fr
%

% Turn points into row if needed
xp = x; dims = size(xp);
if (dims(2) == 1)
   xp = xp';
end

PL = zeros(N+1,length(xp));

gamma0  = NormPoly(0,alpha,alpha,'LAG');
PL(1,:) = 1.0;
if (N==0)
   P = PL'./sqrt(gamma0);
   return;
end

gamma1  = NormPoly(1,alpha,alpha,'LAG');
PL(2,:) = (alpha+1-xp);
if (N==1)
   P = PL(N+1,:)'./sqrt(gamma1);
   return;
end

% Forward recurrence
gammaN  = NormPoly(N,alpha,alpha,'LAG');
for i = 1:N-1
   PL(i+2,:) = (-(i+alpha).*PL(i,:) + (2.*i+1+alpha-xp).*PL(i+1,:))./(i+1);
end
P = PL(N+1,:)'./sqrt(gammaN);

return;
