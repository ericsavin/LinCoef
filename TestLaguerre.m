clear all;
close ('all');

% Choose parameters
alpha = -.9;
TYPE  = 'LAG';
JJ    = 12;
KK    = 13;

JPK   = JJ+KK;
NJ    = JPK; % Max order of the polynomials
nh    = JPK; % We integrate exactly polynomials of max order 2*nh+1
nh1   = nh+1;

% Test orthonormality of generalized Laguerre polynomials
% Check that diag(YLYL) is one
[XLag,WLag] = GNodeWt(nh,alpha,alpha,TYPE);
[YL(1:NJ+1,:),ZL(1:NJ+1,:)] = PGLagn(XLag,NJ,alpha);
YLYL = YL*diag(WLag)*YL';
clear ZL;

% 3rd moment by Gauss-Laguerre quadrature
for LL = 0:JPK
   momLA(LL+1) = sum(YL(JJ+1,1:nh1).*YL(KK+1,1:nh1).*YL(LL+1,1:nh1).*WLag(1:nh1)');
end
% 3rd moment by linearization formula
momLB = LinCoef(JJ,KK,alpha,alpha,TYPE);

