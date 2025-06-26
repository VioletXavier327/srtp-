function [phi,c,measdof,lam,J] = glueref(dof,phic,Iref)
% phi = glueref(dof,phic,Iref) glues mode shape using traditional method by
% calculating the multipliers by least square. Iref gives the setup number
% such that all other setups have at least one common dof with it.
%
% Note:
% 1. If a dof is measured by more than one channels in the same setup, only
% the information from only the last channel is used

% Written by SK Au, CityU
% 100315 - original

ns = length(dof);

if length(phic(:))~=ns
  error('phic must have the same no. of elements as dof.');
end

% determine measured dofs, following user numbering
%---------------------------------------------------
measdof = [];
for i=1:ns
  measdof = union(measdof,dof{i}(:).');
end

n = length(measdof(:));

% convert to dof to consecutive numbers
dof_user = dof;
for i=1:ns
  for j=1:n
    I = dof_user{i}(:)==measdof(j);
    dof{i}(I) = j;
  end
end

% normalize all experimental mode shapes
for i=1:ns
  phic{i} = normc(phic{i});
end

% check that all setups have the same no. of modes
m = size(phic{1},2); % no. of modes, same for all setups
nc = zeros(1,ns); % no. of meas channels in each setup
for i=1:ns
  [nc(i),mi] = size(phic{i});
  if m~=mi
    error('all data sets must have the same no. of modes.');
  end
end

% determine scaling factor c(:)
%=======================================================================
c = ones(ns,m);
Ic = cell(1,ns);
for i = 1:m
  for ii=1:ns
    [Ic{ii},I1,I2] = intersect(dof{Iref},dof{ii});
    c(ii,i) = (phic{Iref}(I1,i).'*phic{ii}(I2,i))/norm(phic{ii}(I2,i))^2;
  end
end

p = cell(1,m);
for i=1:m
  p{i} = zeros(n,ns);
  for ii=1:ns
    p{i}(dof{ii},ii) = phic{ii}(:,i)*c(ii,i);
  end
end

% determine the no. of times a dof appears in different setups
% must be based on dof_unique, otherwise will double count meas. dofs that
% are measured by more than one channels in a setup.

dof_unique = cell(size(dof));
for ii=1:ns
  dof_unique{ii} = unique(dof{ii});
end
nn = zeros(n,1);
for i=1:n
  nn(i) = nnz(find([dof_unique{:}]==i));
end

% nn = zeros(n,1);
% for i=1:n
%   nn(i) = nnz(find([dof{:}]==i));
% end


% glue mode shape
phi = zeros(n,m);
for i=1:m
  phi(:,i) = sum(p{i},2)./nn;
end

phi = normc(phi); % normalize

  % back calc. c(:) and lam(:)
  c = zeros(ns,m); lam = zeros(ns,m);
  for i=1:ns
    ptmp = sum(phi(dof{i},:).*phic{i},1); % inner product
    c(i,:) = sign(ptmp).*sqrt(sum(phi(dof{i},:).^2,1));
    lam(i,:) = 1 - abs(ptmp./c(i,:));
  end

    % calc. measure of fit from first principle
  J = zeros(1,m);
  for i=1:ns
    ntmp = length(dof{i}(:));
    J = J + sum((phi(dof{i},:)-repmat(c(i,:),[ntmp,1]).*phic{i}).^2,1);
  end
