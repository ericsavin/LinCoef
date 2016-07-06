function [x,w] = LaguerreGQ(N,alpha)
% Evaluate N-th order Gauss-Laguerre quadrature points x and weights w
% associated with generalized Laguerre polynomials of type alpha > -1.0
% Note: if N is the order, length(x) = N+1
%
% ESavin / April 2016 / eric.savin@onera.fr
%

if (N==0)
   x(1) = 0;
   w(1) = gamma(alpha+1.0);
   return;
end

% Form symmetric Jacobi matrix
J = zeros(N+1);
J = diag(2.*[0:N]+alpha+1)+diag(-sqrt([1:N].*([1:N]+alpha)),1)+diag(-sqrt([1:N].*([1:N]+alpha)),-1);

% Compute quadrature by eigenvalue solve
[V,D] = eig(J);
x     = diag(D);
w     = (V(1,:)').^2.*gamma(alpha+1.0);

return;
