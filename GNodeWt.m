function [Xinter,Winter] = GNodeWt(NN,alpha,beta,TYPE)
% Compute Gauss quadrature sets from Jacobi, Hermite and Laguerre polynomials
%
% INPUT ---
% NN          : order of interpolation. N1 = NN+1 is the number of points
% alpha, beta : parameters of the probability law
% TYPE        : type of polynomials
%               [JAC: Jacobi, LAG: Laguerre, HER: Hermite]
%
% OUTPUT --
% Xinter      : N1 quadratures nodes 
% Winter      : N1 quadrature weights
%
% ESavin / March 2016 / eric.savin@onera.fr
%
% Licensing: This code is distributed under the CeCILL-C license.
%

N1 = NN+1;
N2 = 2.*NN;

if (alpha==-0.5 & beta==alpha) % Chebyshev quadratures using exact formula
 
   IJ = (0:NN)';
   switch TYPE
      
      case 'JAC' % Gauss quadrature (exact for orders up to 2*NN+1)
         Xinter = -cos((2.*IJ+1).*pi./(N2+2));
	     Winter = pi./N1.*ones(N1,1);
         mu0 = 2.^(alpha+beta+1).*gamma(alpha+1).*gamma(beta+1)./gamma(alpha+beta+2);
         
      otherwise
         error('GNodeWt: this option is not supported');
   end

else

   switch TYPE
      
      case 'JAC' % Gauss-Jacobi quadrature (exact for orders up to 2*NN+1)
         [Xinter,Winter] = JacobiGQ(NN,alpha,beta);
         mu0 = 2.^(alpha+beta+1).*gamma(alpha+1).*gamma(beta+1)./gamma(alpha+beta+2);
      
      case 'HER' % Gauss-Hermite quadrature (exact for orders up to 2*NN+1)
         [Xinter,Winter] = HermiteGQ(NN,alpha);
         mu0 = gamma(alpha+0.5);
      
      case 'LAG' % Gauss-Laguerre quadrature (exact for orders up to 2*NN+1)
         [Xinter,Winter] = LaguerreGQ(NN,alpha);
         mu0 = gamma(alpha+1);
      
      otherwise
         error('GNodeWt: this option is not supported');
   end
   
end

%%% Zero-th moment
if (abs(sum(Winter) - mu0) > 1.0E-10)
   error('GNodeWt: sum of the weights %f is not zero-th moment mu0 = %f.\n',sum(Winter),mu0);
end

return;
      
