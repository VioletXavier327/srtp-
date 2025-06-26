function nstr = num2str2(n,nd)
% num2str2(n,nd) returns a string of length nd the number n padded with zeros.
% e.g., num2str(13,3) gives '013'
% when n requires more than nd digits num2str2 will return simply the
% string of n
%
% This function is useful for generating a sequence of file names

% Written by Siu-Kui Au, CityU
% Ver 1.0 original code, 1 Sep 08

nstr = num2str(n);

m = length(nstr);

if m<nd
  nstr = [repmat('0',1,nd-m),nstr];
end