function [S] = ProdPoch(x,n)
% Compute the Pochhammer symbol (x)_n for a sequence (x_j)
% Implementation of the formula (x)_n = x(x+1)(x+2)...(x+n-1)
% for n integer solely with due consideration of some special cases
%
% ESavin / April 2016 / eric.savin@onera.fr
%

S = prod(pochhammer(x,n));

return;

