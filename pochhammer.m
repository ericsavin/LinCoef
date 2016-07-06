function [S] = pochhammer(x,n)
% Compute the Pochhammer symbol (x)_n
% Implementation of the formula (x)_n = x(x+1)(x+2)...(x+n-1)
% for n integer solely with due consideration of some special cases
% Vectorized version of pochhammer(x,n)
%
% ESavin / April 2016 / eric.savin@onera.fr
%

if (n < -1)
   error('pochhammer: the case of a negative integer input is not supported');
end

nx = length(x);

if (n==0)
   S = ones(size(x));
   return;
end

if (n==1)
   S = x;
   return;
end

if (n==-1)
   S = 1./(1-x);
   return;
end

XX = repmat(reshape(x,nx,1),1,n);
NN = repmat((0:n-1),nx,1);
S  = reshape(prod(XX'+NN'),size(x));

return;

