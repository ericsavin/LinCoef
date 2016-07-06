function [x,w] = HermiteGQ(N,alpha)
% Evaluate N-th order Gauss-Hermite quadrature points x and weights w
% associated with generalized Hermite polynomials of type alpha > -0.5
% Note: if N is the order, length(x) = N+1
%
% ESavin / April 2016 / eric.savin@onera.fr
%

if (N==0)
   x(1) = 0;
   w(1) = gamma(alpha+0.5);
   return;
end

% Form symmetric Jacobi matrix
J     = zeros(N+1);
h1    = (1:N) + alpha.*(1 - power(-1,1:N));
J     = diag(sqrt(h1./2),1);
J     = J+J';

% Compute quadrature by eigenvalue solve
[V,D] = eig(J);
x     = diag(D);
w     = (V(1,:)').^2.*gamma(alpha+0.5);

return;
