function r = mac(x,y)
% r = mac(x,y) calculates the modal assurance criteria ratio
%
% x = (n,m1), y = (n,m2), => r = (m1,m2)
normx = sqrt(sum(x.^2,1)); % (1,m1)
normy = sqrt(sum(y.^2,1)); % (1,m2)
r = (x.'*y)./(normx.'*normy); % (m1,m2)
