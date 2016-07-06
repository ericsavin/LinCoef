function [gamu] = GaMu(n,alpha)
% Evaluate a special function for the linearization coefficients
% of Hermite polynomials
%
% ESavin / April 2016 / eric.savin@onera.fr
%

ip   = mod(n,2);
m    = floor(n/2);
gamu = power(2,n).*factorial(m).*pochhammer(alpha+0.5,m+ip);

return;

