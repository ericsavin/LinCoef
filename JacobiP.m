function [P] = JacobiP(x,alpha,beta,N)
% Evaluate Jacobi polynomial of type (alpha,beta) > -1
% at points x for order N
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

% Initial values P_0(x) ans P_1(x)
if (alpha == -0.5 & beta == -0.5)
   gamma0 = pi;
else
   gamma0 = 2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)* ...
         gamma(beta+1)/gamma(alpha+beta+1);
end
PL(1,:) = 1.0/sqrt(gamma0);
if (N==0)
   P = PL';
   return;
end

gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
PL(2,:) = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/sqrt(gamma1);
if (N==1)
   P = PL(N+1,:)';
   return;
end

% Repeat value in recurrence
aold = 2/(alpha+beta+2)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

% Forward recurrence using the symmetry of the recurrence
for i = 1:N-1
   h1 = 2*i + alpha + beta;
   anew = 2/(h1+2)*sqrt((i+1)*(i+1+alpha+beta)*(i+1+alpha)* ...
          (i+1+beta)/(h1+1)/(h1+3));
   bnew = -(alpha^2-beta^2)/h1/(h1+2);
   PL(i+2,:) = 1/anew*(-aold*PL(i,:) + (xp-bnew).*PL(i+1,:));
   aold = anew;
end
P = PL(N+1,:)';
return;
