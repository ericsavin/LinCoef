function [HG] = HyGJacobi(n,a,b)
% Evaluate the generalized hypergeometric function pFq(a,b,x) at x
% for a(1:p) and b(1:q) by direct summation
%
% ESavin / April 2016 / eric.savin@onera.fr
%

Aj = 1.0;
HG = Aj;

for j = 0:n
   Aj = Aj.*(-n+j).*prod(a+j)./prod(b+j)./(j+1);
   HG = HG + Aj;
end
  
return;

