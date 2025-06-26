function [dof1_user,dof2_user,I1,I2] = refdof(dof)
% [dof1,dof2,I1,I2] = refdof(DOF)
% separates the dofs in DOF into reference dofs (that appeard more than
% once in setups and normal dofs. 
%
% dof = cell(1,ns), dof{i} = (1,ni) gives the measured dofs in setup i
% dof1 = cell(1,ns), dof1{i} = (1,n1i) gives the reference dofs in setup i
% dof2 = cell(1,ns), dof2{i} = (1,n2i) gives the normal dofs in setup i
% I1 and I2 are such that
% dof{i}(I1{i}) = dof1{i}, dof{i}(I2{i}) = dof2{i}

% Written by SK Au, CityU
% 100312 - original code

ns = length(dof(:));

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

% determine ref dofs that appear more than once in setups
J = zeros(n,1); % counts the no. of times dof appears
for i=1:n
  for ii=1:ns
    J(i) = J(i) + sum(dof{ii}==i);    
  end
end
J1 = find(J>1); % dofs that appear more than once in setups
J2 = setdiff(1:n,J1); % dofs that appear only once
% J2 = find(J<=1); % dofs that appear only once

% condense dof and phic to repeated dofs only
dof1 = cell(1,ns); dof2 = cell(1,ns);
I1 = cell(1,ns); I2 = cell(1,ns);
for ii=1:ns
 [dum,I1{ii}] = intersect(dof{ii},J1);
 dof1{ii} = dof{ii}(I1{ii});

 [dum,I2{ii}] = intersect(dof{ii},J2);
 dof2{ii} = dof{ii}(I2{ii});
end

% convert indices to user defined ones
dof1_user = dof1;
dof2_user = dof2;
for i=1:ns
  for j=1:n
    II = dof1{i}(:)==j;
    dof1_user{i}(II) = measdof(j);
    
    II = dof2{i}(:)==j;
    dof2_user{i}(II) = measdof(j);    
  end
end
