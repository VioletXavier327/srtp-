function z = optfield(x,y)
% Z = optfield(X,Y) copies fields in the structure X to Z and those fields
% in Y but not in X to Z.
%
% This function is useful for combining mandatory and optional fields from
% provided from input. For this purpose, X is the input structure, Y is the
% default structure, Z is the resulting structure.

% Written by SK Au, CityU
% 100727 - original

z = y;

if ~isstruct(x)
  return;
end

fx = fieldnames(x);
nx = length(fx(:));
for i=1:nx
  z.(fx{i}) = x.(fx{i});
end

