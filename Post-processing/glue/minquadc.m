function [x,v0,d0] = minquadc(A,B,c)
% x = minquadc(A,B) computes the unique minimizing solution of 
% Q = 0.5*x.'*A*x + x.'*B, subjected to the constraint x.'*x = c
% where x = (n,1), A = (n,n) symmetric positive definite, B = (n,1), c>0.
% B can be given as either a column or row vector.
% 
% The solution x satisfies the Normal equation:
% A*x + B - d*x = 0 with x.'*x = 1, where d is an eigenvalue
% 
% x = minquadc(A,B) assumes c = 1
%
% The solution of x is // to the first half of the eigenvector
% corresponding to the smallest eigenvalue of the matrix
% D = [A B*B.'/c; eye(n), A]
%
% [x,v,d]=... also returns the real eigenvalues (d) and eigenvectors (v) of D
% 
% For the theory of this problem, See Gander 1981, 1989

% Written by SK Au, CityU
% 15 Feb 09, original

if nargin<3
  c = 1;
elseif length(c(:))~=1
  error('c must be a scalar');
elseif ~(c>0)
  error('c must be positive');
end

[nA,mA]=size(A);
if nA~=mA
  error('A must be a square matrix');
end

nB = length(B(:));
if nB~=nA
  error('The length of B must be equal to the first dimension of A');
end

D = [A, B(:)*B(:).'/c; eye(nA), A];

[v0,d0] = eig(D);

v = v0;
d = d0;

d = diag(d);

% retain only modes with real eigenvalues
imagd = imag(d);
I = find(~imagd); % index of modes with real eigenvalues
d = d(I);
v = v(:,I);

% sort modes
[d,I] = sort(d);
v = v(:,I);

x = v(1:nA,1);
y = v(nA+1:end,1);
x = x/(B.'*y);

