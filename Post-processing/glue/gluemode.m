function [phi,c,measdof,lam,J,niter] = gluemode(dof,phic,options)
% phi = gluemode(dof,phic) combines the modeshapes in phic to obtain the
% best least-square fit modeshape phi. 
%
% dof = (1,ns) cell, dof{i} = (1,nci) array giving measured dofs in the
% i-th data set; nci = no. of channels of i-th data set
%
% phic = (1,ns) cell, phic = (nci,m), each column gives a modeshape; m is
% the no. of modes, which must be the same for all data sets.
%
% phi = (n,m) modeshape where n is the largest dof no. given among all
% dof{:}. If the i-th dof is not measured, then phi(i,:) = nan
%
% [phi,c]=gluemode(...) also returns the scaling factors c = (ns,m)
% 
% [phi,c,measdof]=gluemode(...) also returns the user-defined dofs
%
% [phi,c,measdof,lam,J]=gluemode(...) also returns the user-defined dofs
%
% By default, the algorithm solves the condensed system of ref dofs
% To solve the original full system (which takes much longer time), use
% [...] = gluemode(...,1)
%
% notes:
% 1. the dofs need not be consecutive
% 2. the meas dofs need not be unique
% 3. the mode shapes need not be normalized
% 4. a given dof can be measured by multiple channels
% 5. The min. least square measure of fit is 
%      J=2*sum(lam.*c.^2,1)
%    The fractional measure of fit is 
%      J./sum(c.^2,1)=2*sum(lam.*normc(c).^2,1)

% Written by SK Au, CityU
%
% 090316 - original, gluing modeshapes requires iteration
% 090511 - correct theory of glueing modeshapes, no iteration needed
% 091226 - corrected theory of glueing modeshapes, iteration needed to make
%          sure empirical modeshapes in differetn setups and theoretically
%          implied ones have the same norm 
% 100119 - handled trivial case of ns = 1
% 100227 - fixed bug of L when meas dof not unique. 
%          now meas dof need not be unique 
% 100228 - fixed bug of converting user dof numbering to consecutive ones
% 100309 - used optimality to update lam(:); convergence based on lam(:)
%          speed up the calculation of A; cleaned up code
% 100310 - found problem re calc. of A. make phic{ii}(:,j).'/(1-lam(:,j))
%          to be calc first; otherwise may have convergence problem in some
%          cases
% 100318 - debugged the case when a dof is measured by more than one
%          channel in the same setup
% 100320 - included condensed algorithm, used when options(1)==1
% 100321 - used new algorithm. should be robust now
% 100311 - debugged the expression of AA; needs to update during iteartion
% 120424 - debugged Ic when single setup and 'lsref'

maxiter = 500;
tol = 1e-6; % convergence tolerance in mac over successive mode shapes

if nargin<3
  options = 0;
end

ns = length(dof);

if length(phic(:))~=ns
  error('phic must have the same no. of elements as dof.');
end

% determine measured dofs, following user numbering
%---------------------------------------------------
measdof = [];
for i=1:ns
%   measdof = union(measdof,dof{i}(:).');
  measdof = union(measdof,dof{i}(:));
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

if ~options(1) % condensed algorithm
  % separate into ref dofs and normal dofs
  [dof1,dof2,Ic1,Ic2] = refdof(dof);
  
  % Ic collects the ref dofs and ind. dofs
  dof1 = unique([dof1{:}]); dof1 = dof1(:).';
  Ic = {dof1,dof2{:}}; % (1,1+ns) cell
  nr = length(Ic{1}(:)); % no. of ref. dofs    
%   Ic = {unique([dof1{:}]),dof2{:}}; % (1,1+ns) cell
%   nr = length(Ic{1}(:)); % no. of ref. dofs  
end

% compute phi
%-------------------
phi = nan(n,m);
c = zeros(ns,m);
lam = zeros(ns,m);

L = cell(1,ns);
LtL = cell(1,ns);
% LtL = zeros(n,n,ns);
LtP = cell(1,ns);
% LtP = zeros(n,m,ns);
for ii=1:ns
 L{ii} = zeros(nc(ii),n);
 for ic = 1:nc(ii)
   L{ii}(ic,dof{ii}(ic)) = 1;
 end
 LtL{ii} = L{ii}.'*L{ii}; % (n,n)
%   LtL(:,:,ii) = L{ii}.'*L{ii}; % (n,n)
 LtP{ii} = L{ii}.'*phic{ii};
%   LtP(:,:,ii) = L{ii}.'*phic{ii}; % (n,m)
end

% determine optimal mode shapes for every mode
%========================================================================
J = zeros(1,m); c = zeros(ns,m); lam = zeros(ns,m);
niter = zeros(1,m);
for j=1:m % different modes

  LPPL = cell(1,ns);
  for ii=1:ns
    LPPL{ii} = phic{ii}(:,j).'*L{ii};
    LPPL{ii} = LPPL{ii}.'*LPPL{ii}; % (n,n)
  end
  
  % obtain initial guess for phi(:,j)
  A = zeros(n,n);
  for ii=1:ns
    A = A + LtL{ii} - LPPL{ii}; 
%     A = A + LtL(:,:,ii) - LPPL{ii}; 
  end
  [v,d]=eig(A);
  d = diag(d);
  [d,I]=sort(d);
  phi(:,j) = v(:,I(1));

  % back calc. c(:) from first principle
  for ii=1:ns
%     keyboard
    ptmp = sum(phi(dof{ii},j).*phic{ii}(:,j),1); % inner product
    c(ii,j) = sign(ptmp).*norm(phi(dof{ii},j));
  end

  if ~options(1) % condensed algorithm
    % compute Lc  = (n,nc);  Bc = Lc.'*Lc = (nc,nc)
    Lc = eye(nr);
    for ii=1:ns
      if ~isempty(dof2{ii})
        Lc = blkdiag(Lc,normc(phic{ii}(Ic2{ii},j)));
      end
    end
  end
  
  
  
  % iterate to get optimal phi(:,j)
  %==============================
  for k=1:maxiter
    
    % calc. optimal c(:)
    for ii=1:ns
      ptmp = sum(phi(dof{ii},j).*phic{ii}(:,j),1); % inner product
      c(ii,j) = sign(ptmp).*norm(phi(dof{ii},j));
      lam(ii,j) = 1 - abs(ptmp./c(ii,j));          
    end

    % update AA
%     AA = sum(repmat(reshape(1+lam(:,j),[1 1 ns]),[n,n,1]).*LtL,3); % (n,n)
    AA = zeros(n,n);
    for ii=1:ns
      AA = AA + (1+lam(ii,j))*LtL{ii}; 
    end

    % update BB
%     BB = -sum(repmat(c(:,j).',[n 1]).*squeeze(LtP(:,j,:)),2); % (n,1)  
    BB = zeros(n,1);
    for ii=1:ns
      BB = BB - LtP{ii}(:,j)*c(ii,j);
    end
    
    if ~options(1) % condensed algorithm       
%       keyboard
      I = cell2mat(Ic); % (1,n)
      AAc = Lc.'*AA(I,I)*Lc;        
      BBc = Lc.'*BB(I);
%       keyboard
      v = minquadc(AAc,BBc);
      phi([Ic{:}],j) = Lc*v;
    else % full dim.
       phi(:,j) = minquadc(AA,BB);
    end
    
    % check convergence of phi(:,j)
    %================================
    if k==1
      phi_last = phi(:,j);
%       c_last = c(:,j);
    else
      eq = 1-abs(mac(phi(:,j),phi_last));
      
      if abs(eq)<tol
        break;
      else
        if k<maxiter
          phi_last = phi(:,j);
%           c_last = c(:,j);      
        end
      end
    end

  end % k
  
  niter(j) = k;
  
  if (k==maxiter)
    warning(['gluemode does not converge within maxiter: ','mode no.',num2str(j)]);
  end    

  % back calc. c(:) from first principle
  for ii=1:ns
    ptmp = sum(phi(dof{ii},j).*phic{ii}(:,j),1); % inner product
    c(ii,j) = sign(ptmp).*norm(phi(dof{ii},j));
    lam(ii,j) = 1 - abs(ptmp./c(ii,j));    
  end
  
  
  % calc. measure of fit from first principle
  J(j) = 0;
  for ii=1:ns
    ntmp = length(dof{ii}(:));
    J(j) = J(j) + ...
      sum((phi(dof{ii},j)-c(ii,j)*phic{ii}(:,j)).^2,1);
  end
  
end % for j


%====================================================================
% private functions
%====================================================================
    