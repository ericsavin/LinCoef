clear all;
close ('all');

% Choose parameters
alpha = 4.7;
TYPE  = 'HER';
JJ    = 20;
KK    = 10; % max about 20

JPK   = JJ+KK;
NJ    = JPK;   % Max order of the polynomials
nh    = JPK;   % We integrate exactly polynomials of max order 2*nh+1
nh1   = nh+1;

% Test orthonormality of generalized Hermite polynomials
% Check that diag(YHYH) is one
[XHer,WHer] = GNodeWt(nh,alpha,alpha,TYPE);
[YH(1:NJ+1,:),ZH(1:NJ+1,:)] = PGHern(XHer,NJ,alpha);
YHYH = YH*diag(WHer)*YH';
weight = sum(WHer);
clear ZH;

% 3rd moment by Gauss-Hermite quadrature
for LL = 0:JPK
   momHA(LL+1) = sum(YH(JJ+1,1:nh1).*YH(KK+1,1:nh1).*YH(LL+1,1:nh1).*WHer(1:nh1)');
end
% 3rd moment by linearization formula
momHB = LinCoef(JJ,KK,alpha,alpha,TYPE);
