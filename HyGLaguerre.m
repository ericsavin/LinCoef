function [HG] = HyGLaguerre(n,j,k,alp)
% Evaluate the terminating hypergeometric series pFq(a,b,1)
% for a(1:p) and b(1:q) by direct summation
%
% ESavin / April 2016 / eric.savin@onera.fr
%

nn = j + k - n;
lm = max(max(j-n,k-n),0);
lM = floor(nn/2);
Aj = pochhammer(-(nn-1)/2,lm).*pochhammer(-nn/2,lm).*pochhammer(alp+n,lm) ...
     ./factorial(n-j+lm)./factorial(n-k+lm)./factorial(lm);
HG = Aj;

for ll = (lm+1):lM
   Aj = Aj.*(-(nn-1)/2+ll-1).*(-nn/2+ll-1).*(alp+n+ll-1)./(n-j+ll)./(n-k+ll)./ll;
   HG = HG + Aj;
end
  
return;

