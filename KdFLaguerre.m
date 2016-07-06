function [KdF] = KdFLaguerre(n,b,c,d)
% Evaluate the terminating Kampe de Feriet series p:0Fr:s
% at x=1, y=1 by direct summation
% Linearization coefficient of generalized Laguerre polynomials
%
% ESavin / June 2016 / eric.savin@onera.fr
%

KdF = 0.;
for ir = 0:n
   ad   = pochhammer(-n,ir)./ProdPoch(d,ir);
   AJK  = 0.;
   for is = 0:ir
      bc  = ProdPoch(b,is).*ProdPoch(c,ir-is)./factorial(is)./factorial(ir-is);    
      AJK = AJK + ad.*bc;
   end
   KdF  = KdF + AJK;
end
  
return;

