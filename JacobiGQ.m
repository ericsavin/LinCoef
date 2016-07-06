function [x,w] = JacobiGQ(N,alpha,beta)
% Evaluate N-th order Gauss quadrature points x and weights w
% associated with Jacobi polynomials of type (alpha,beta) > -1
% Note: if N is the order, length(x) = N+1
%
% ESavin / April 2016 / eric.savin@onera.fr
%

if (N==0)
   x(1) = (beta-alpha)/(alpha+beta+2);
   w(1) = 2.^(alpha+beta+1).*gamma(alpha+1).*gamma(beta+1)./gamma(alpha+beta+2);
   return;
end

% Form symmetric matrix for recurrence
J  = zeros(N+1);
h1 = 2*(0:N) + alpha + beta;
J1 = [-1/2*(alpha-beta)./(h1(1)+2)  -1/2*(alpha^2-beta^2)./(h1(2:(N+1))+2)./h1(2:(N+1))];
J  = diag(J1) + ...
     diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha+beta).* ...
     ((1:N)+alpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)),1);
J = J+J';

% Compute quadrature by eigenvalue solve
[V,D] = eig(J); x = diag(D);
w = (V(1,:)').^2*2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)* ...
    gamma(beta+1)/gamma(alpha+beta+1);

return;
