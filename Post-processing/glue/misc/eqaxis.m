function eqaxis(nrow,ncol,I,opt)
% EQUAXIS(NROW,NCOL,I) sets the limits of the subplots I(:) to be the min. 
% required to hold all the entries in the subplots.
%

xyz = [1 1 1];
if nargin>3
  switch lower(opt)
  case 'x',
    xyz = [1 0 0];
  case 'y',
    xyz = [0 1 0];
  case 'z',
    xyz = [0 0 1];
  case 'xy',
    xyz = [1 1 0];
  case 'yz',
    xyz = [0 1 1];
  case 'xz',
    xyz = [1 0 1];
  otherwise
    error('Unknow option');
  end
end

n = length(I);

subplot(nrow,ncol,I(1));
v = axis;

for k=2:n
  subplot(nrow,ncol,I(k));
  axisv = axis;
  if xyz(1)
    v(1) = min(v(1),axisv(1));
    v(2) = max(v(2),axisv(2));
  end
  if xyz(2)
    v(3) = min(v(3),axisv(3));
    v(4) = max(v(4),axisv(4));
  end

  if xyz(3)
    if length(v)>4
      v(5) = min(v(5),axisv(5));
      v(6) = max(v(6),axisv(6));
    end
  end
end 

for k=1:n
  subplot(nrow,ncol,I(k));
  axis(v);
end
