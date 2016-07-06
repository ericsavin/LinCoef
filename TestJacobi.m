clear all;
close ('all');

% Choose parameters
alpha = -.4;
beta  = 10.0;
TYPE  = 'JAC';
JJ    = 17;
KK    = 8;

JPK   = JJ+KK;
NJ    = JPK; % Max order of the polynomials
nh    = JPK; % We integrate exactly polynomials of max order 2*nh+1
nh1   = nh+1;

% Test Chebyshev
if (alpha==-.5 & beta==-.5)
   TYPE = 'CHE';
end

% Test orthonormality of Jacobi polynomials
% Check that diag(YJYJ) is one
[XJac,WJac] = GNodeWt(nh,alpha,beta,'JAC');
sigma = ones(NJ+1,1);
[YJ(1:NJ+1,:),ZJ(1:NJ+1,:)] = PJacn(XJac,NJ,alpha,beta,sigma);
YJYJ = YJ*diag(WJac)*YJ';
clear ZJ;

% 3rd moment by Gauss-Jacobi quadrature
for LL = 0:JPK
   momJA(LL+1) = sum(YJ(JJ+1,1:nh1).*YJ(KK+1,1:nh1).*YJ(LL+1,1:nh1).*WJac(1:nh1)');
end
% 3rd moment by linearization formula
momJB = LinCoef(JJ,KK,alpha,beta,TYPE);

% Test orthonormality of Gegenbauer polynomials
% Check that diag(YGYG) is one
[XGeg,WGeg] = GNodeWt(nh,alpha,alpha,'JAC');
sigma = ones(NJ+1,1);
[YG(1:NJ+1,:),ZG(1:NJ+1,:)] = PJacn(XGeg,NJ,alpha,alpha,sigma);
YGYG = YG*diag(WGeg)*YG';
clear ZG;

% 3rd moment by Gauss-Jacobi quadrature
for LL = 0:JPK
   momGA(LL+1) = sum(YG(JJ+1,1:nh1).*YG(KK+1,1:nh1).*YG(LL+1,1:nh1).*WGeg(1:nh1)');
end
% 3rd moment by linearization formula
% WARNING *** Does not work with Chebyshev polynomials
momGB = LinCoef(JJ,KK,alpha,alpha,'GEG');