function [KdF] = KdFJacobi(n,a,b,c,d,e,f)
% Evaluate the terminating Kampe de Feriet series p:qFr:s
% at x=1, y=1 by direct summation
% Linearization coefficient of Jacobi polynomials 
%
% ESavin / June 2016 / eric.savin@onera.fr
%
aa = [-n,a];
KdF  = 0.0;
for ir = 0:n
   ad   = ProdPoch(aa,ir)./ProdPoch(d,ir);
   AJK  = 0.;
   for is = 0:ir
      bcef = ProdPoch(b,is).*ProdPoch(c,ir-is)./ProdPoch(e,is)./ProdPoch(f,ir-is)./factorial(is)./factorial(ir-is);
      AJK  = AJK + ad.*bcef;
   end
   KdF  = KdF + AJK;
end
  
return;

