function o = rptsetup(varargin)
% o = RPTSETUP prompts the user to select setup file (mat)
% and combine their information to give overall modal id results. 
%
% Help glue/contents for the format of setup files
%
% The followings are performed:
% 1. glue mode shape by global least square method
% 2. scale PSD of modal of modal force of each setup to be consistent w/
%    scaling of global mode shape
% 3. combine individual stastistics of modal properties to give an
%    overall statistics, assuming time-invariance over setups
%
% o = RPTSETUP(IN) supplies the structure IN with the following optional
% fields:
% - xls = char of setup file (mat) name. [] prompts the user
% - imode = array of modes to process. [] means all modes
% - isetup = array of setup no. to process. [] plots undeformed shape only
%            if not provided, all valid setups will be used
%            The setup no. is the one stated in xls
% - flag_disp = 1 to display steps; 0 not to
% - gluemethod = 'lsall', 'lsref','ref'
% - flag_plt = 1 to plot setup results; 0 not to
% - ndigits
%
% By default, 
% xls = [], imode = [], flag_disp = 0, 
% gluemethod='lsref', isetup not exist, flag_plt = 0, ndigits = 2
% 
% The output structure o contains the following fields:
% 1. setup = structure containing fields .mpv and .cov of modal
%  properties f,z,S,rms,phi of individual setups:
%  [n = total #dofs, m= #modes]
%    .f = (1,m), natural frequency (Hz)
%    .z = (1,m), damping ratio 
%    .S = (1,m), PSD of modal force scaled
%    .rms = (1,m), rms of modal response
%    .phi = (n,m)
%  
%  2. stats = structure containing .mpv and .cov of overall statistics
%    .mpv is the average weighted by inverse of variance
%    .cov is the reciprocal of the sum of inverse of variance
%
%  3. phi0 = (n,m), global mode shape, 
%     dof0 = (n,1) array of dof no. consistent w/ phi
%     c = (ns,m), lam = (ns,m) see GLUEMODE
%
%  Notes:
%  1. all quantitites (S, rms, phi) are scaled so that global mode
%     has unit norm
%  2. for setup.mpv.phi, when a dof is meas. by more than one ch. in a 
%     given setup, the result is averaged
%  3. If the result file for a given setup cannot be found, the setup will
%     be automatically removed
%  4. In setup (mat) file, if mat = {} or missing, null results will be
%     returned with warning.

% Written by SK Au, TCU, Liverpool
% 130524 - updated to be consistent w/ bayoma format
% 121016 - debugged scaling of rms
% 121012 - debugged nan stats
% 120425 - debugged setups w/o result files
% 120420 - original code, adapted from some part of gui_pltmode3

in = sinput(varargin{:});

% contants and defaults
tag.prefix = '[FILE]';
tag.setup = {'[SETUP]','[ENDSETUP]'};
tag.loc = {'[LOC]','ENDLOC'};
tag.line = {'[LINE]','ENDLINE'};

in_default.xls = [];
in_default.imode = [];
in_default.flag_disp = 1;
in_default.gluemethod = 'lsref';
in_default.flag_plt = 1;
in_default.ndigits = 2;
in_default.xls_result = 'result.xls';

o = [];
var_name = {'f','z','S','rms'};
var_unit = {'Hz','','g^2/Hz','g'};

if nargin<1
  in = in_default;
end

in = optfield(in,in_default);

imode = in.imode;

if isempty(in.xls)
  in.xls = uigetfile2('*.xls','Pick setup xls file.'); 
end
if isempty(in.xls)
  return;
end
if iscell(in.xls)
  in.xls = in.xls{1}; % char
end

switch in.gluemethod
  case {'lsall','lsref'}
  case 'ref'
    if ~isfield(in,'Iref')
      warning('Iref not provided; gluemethod forced to ''lsref''.');
      in.gluemethod = 'lsref';
    end
  otherwise
    error('Unknown GLUEMETHOD');
end

% result file prefix
if ~isfield(in,'prefix')
  [dum,in.prefix] = xlstag(in.xls,tag.prefix,[1 1]);
  if iscell(in.prefix)
    in.prefix = in.prefix{1}; % char
  end
end

if isempty(in.prefix)
  error('Tag ''[FILE]'' cannot be found in XLS');
end

% load setup info.
%=================
setupinfo = xlstag(in.xls,tag.setup{1},tag.setup{2},'cleanrow','cleancol');
setup_num = setupinfo(:,1).';
dof = setupinfo(:,2:end);

% generate mat file for all setups
mat = cell(1,length(setup_num(:)));
% ndigits = ceil(log10(length(mat(:))))+1;
% ndigits = max(ndigits,2);

% ==============================================================================
% 211009 Revise: allow the users not nessarily name the result file as
% "result01" when single setup
% ===============================================================================
if length(mat(:))==1
    mat = uigetfile2('*.mat','Pick result file to view.');
else
    for i=1:length(mat(:))
        mat{i} = [in.prefix,num2str2(setup_num(i),in.ndigits)];
    end
end
% =================================================================================

% old format:
% =================================================================================
%     for i=1:length(mat(:))
%         mat{i} = [in.prefix,num2str2(setup_num(i),in.ndigits)];
%     end
% =================================================================================

% keyboard

% in.isetup and ns
%------------------
if ~isfield(in,'isetup')
  in.isetup = setup_num(:).';
end
isetup_user = in.isetup; % stored for display only
in.isetup = myfun_renum(in.isetup,setup_num); % consecutively numbered
%
I = find(isnan(in.isetup));
if ~isempty(I)
  disp('Warning: The following setups not defined in xls are removed:');
  disp([num2str(isetup_user(I))]);
  in.isetup(I) = [];
end
ns = length(in.isetup(:));

% check if result file exist; delete the setup if not
Imiss = [];
for i=1:ns
  [pp,name] = fileparts(mat{in.isetup(i)});
  if ~exist([fullfile(pp,name),'.mat'],'file')
%     Imiss = [Imiss,setup_num(in.isetup(i))];    
    Imiss  = [Imiss,i];
%     in.isetup(i) = nan;
  end
end
% if any(isnan(in.isetup))
if ~isempty(Imiss)
  disp('Warning: Result file for the following setups are missing:');
  disp(num2str(setup_num(Imiss)));
%   disp(mat(in.isetup(Imiss)));
%   for i=1:length(Imiss)
%     fprintf('%s ',mat{in.isetup(Imiss(i))});
%   end
%   fprintf('\n');
end
% in.isetup(isnan(in.isetup))=[]; % retain only setups wih results
in.isetup(Imiss) = [];
ns = length(in.isetup(:));
mat = mat(in.isetup); % retain only selected setups

% keyboard
% load dof info; cellref.dof
%--------------------------------------------------------------------
dof = num2cell(dof,2); % cell w/ 1 col.
dof = dof(in.isetup); % retain only selected setups
for i=1:ns
  dof{i}(isnan(dof{i}))=[]; % remove iddle channels
end

n = length(unique([dof{:}])); % no. of measured dofs

% load results of different setups
%==================================
s = cell(1,ns);
for i=1:ns
  s{i} = load(mat{i});
end

if ns>0
  m = length(s{1}.f(:)); % no. of modes, same for all setups
  
  % check that all setups have the same no. of modes
  for i=1:ns
    if length(s{i}.f(:))~=m
      error('f must have the same no. of modes in all setups');
    end
    if length(s{i}.z(:))~=m
      error('z must have the same no. of modes in all setups');
    end
    if size(s{i}.phi,2)~=m
      error('phi must have the same no. of modes in all setups');
    end  
    if length(s{i}.rms(:))~=m
      error('rms must have the same no. of modes in all setups');
    end
  end
  
  if isempty(imode)
    imode = [1:m];
  end
  
  if any(imode<0|imode>m)
    error('index in IMODE must be from 1 to max. no. identified.');
  end
  m = length(imode(:));
  
else % nominal value for ns=0
  m = 1;
end

% m = 0;

% initialize
phi0 = nan(n,m); dof0 = nan(n,m); c = nan(ns,m); lam = nan(ns,m);
[setup.mpv.f,setup.mpv.z,setup.mpv.S,setup.mpv.rms,...
 setup.cov.f,setup.cov.z,setup.cov.S,setup.cov.rms] = deal(nan(ns,m));
setup.mpv.phi = cell(ns,1);
[stats.mpv.f,stats.mpv.z,stats.mpv.S,stats.mpv.rms,...
 stats.cov.f,stats.cov.z,stats.cov.S,stats.cov.rms] = deal(nan(1,m));

if ns>0

  % set mpv to nan if offband
  %-----------------------
  for i=1:ns
    if isfield(s{i},'offband')
      offband = s{i}.offband;
      s{i}.f(offband) = nan;
      s{i}.z(offband) = nan;
      s{i}.S(offband,:) = nan; s{i}.S(:,offband) = nan;
      s{i}.Se(offband) = nan;
      s{i}.rms(offband) = nan;
      s{i}.sn(offband) = nan;      
    end
  end 
  
  % fill up setup.mpv & setup.cov
  %-------------------------------
  for i=1:ns
    setup.mpv.f(i,:) = s{i}.f(imode); % (1,m)
    setup.mpv.z(i,:) = s{i}.z(imode); % (1,m)

    if ~iscell(s{i}.S) % format before 130524
      setup.mpv.S(i,:) = diag(s{i}.S(imode,imode)).'; % (1,m)
    else % format after 130524
      setup.mpv.S(i,:) = s{i}.Sii(imode); % (1,m)
    end 
    
    setup.mpv.rms(i,:) = s{i}.rms(imode); % (1,m)
    setup.mpv.Se(i,:) = s{i}.Se(imode); % (1,m)
    setup.mpv.phi{i} = s{i}.phi(:,imode); % (ni,m)

    if isfield(s{i},'coefv')
      if iscell(s{i}.coefv) % format before 130524
        setup.cov.f(i,:) = s{i}.coefv{1}(imode); % (1,m)
        setup.cov.z(i,:) = s{i}.coefv{2}(imode); % (1,m)
        setup.cov.S(i,:) = s{i}.coefv{3}(imode); % (1,m)
        setup.cov.V(i,:) = s{i}.coefv{5}(imode); % (1,m)
        
        if isfield(s{i},'rmscov') % format before 130524
          setup.cov.rms(i,:) = s{i}.rmscov(imode); % (1,m)
        end
        
      else % format after 130524
        setup.cov.f(i,:) = s{i}.coefv.f(imode); % (1,m)
        setup.cov.z(i,:) = s{i}.coefv.z(imode); % (1,m)
        setup.cov.S(i,:) = s{i}.coefv.Sii(imode); % (1,m)
        setup.cov.Se(i,:) = s{i}.coefv.Se(imode); % (1,m)
        
        setup.cov.rms(i,:) = s{i}.coefv.rms(imode); % (1,m)
      end
    end
        
  end
  
  % glue mode shapes
  %==================
  if in.flag_disp
    fprintf('Glueing mode shapes ... ');
  end

  for i=1:m
    if in.flag_disp
      fprintf('%2i/%2i ',[i,m]);
    end
  
    % mode shapes for the i-th mode
    phic_tmp = cell(1,ns);
    for j=1:ns
      phic_tmp{j} = setup.mpv.phi{j}(:,i);
    end

    switch in.gluemethod
      case 'lsall' % glue based on all dofs, slower
        [phi0(:,i),c(:,i),dof0,lam(:,i)] = ...
          gluemode(dof,phic_tmp,1); % least sq.
      case 'lsref' % condense dofs, equivalent but faster
%         keyboard
        [phi0(:,i),c(:,i),dof0,lam(:,i)] = ...
          gluemode(dof,phic_tmp); % least sq.
      case 'ref'
        [phi0(:,i),c(:,i),dof0,lam(:,i)] = ...
          glueref(dof,phic_tmp,in.Iref); % least sq.
      otherwise
        error('Unknown GLUEMETHOD');     
    end
  end % for i=1:m
  if in.flag_disp
    fprintf('done\n');
  end

  % scale results to be consistent w/ global mode shape scaling
  setup.mpv.S = setup.mpv.S./c.^2; % (ns,m)
%   setup.mpv.rms = setup.mpv.rms./c.^2; % (ns,m)
  setup.mpv.rms = setup.mpv.rms./abs(c); % (ns,m)
  
  for i=1:ns
    ni = size(setup.mpv.phi{i},1); % no. of measured dofs
    setup.mpv.phi{i} = setup.mpv.phi{i}.*repmat(c(i,:),[ni,1]);
  end

  % convert setup.mpv.phi to (n,ns,m)
  % when a dof is meas. by more than 1 ch. in a setup, result is averaged
  %----------------------------------------------------------------------
  p = zeros(n,ns,m); % tmp array for phi
  cnt = zeros(n,ns,m); % counts # ch. a dof is measured in a setup
  for i=1:ns
    ni = length(dof{i}(:));
    for j=1:ni
      I = find(dof0==dof{i}(j));
      p(I,i,:) = p(I,i,:) + reshape(setup.mpv.phi{i}(j,:),[1,1,m]);
      cnt(I,i,:) = cnt(I,i,:) + 1;
    end
  end
  I=(cnt~=0);
  p(I) = p(I)./cnt(I);
  p(~I)=nan;
  setup.mpv.phi = p;

  % overall statistics (weighted by reciprocal of variance)
  %--------------------------------------------------------
  tmpname = fieldnames(setup.cov);
  for i=1:length(tmpname(:))
%     xi = setup.mpv.(tmpname{i}); % (ns,m), mpv
%     si = xi.*setup.cov.(tmpname{i}); % (ns,m), std.
%     if any(isnan(si(:)))
%       wi = ones(ns,m)/ns;
%     else
%       wi = 1./si.^2; 
%       wi = wi./repmat(sum(wi,1),[ns,1]); % (ns,m), weights
%     end
%     mpvi = nansum(wi.*xi,1); % weighted mean
%     covi = sqrt(1./nansum(1./si.^2,1))./mpvi; % (1,m), c.o.v.  
%     stats.mpv.(tmpname{i})=mpvi;
%     stats.cov.(tmpname{i})=covi;

    xi = setup.mpv.(tmpname{i}); % (ns,m), mpv
    si = xi.*setup.cov.(tmpname{i}); % (ns,m), std.
    
    wi = 1./si.^2; 
    wi = wi./repmat(nansum(wi,1),[ns,1]); % (ns,m), weights
    mpvi = nansum(wi.*xi,1); % weighted mean
    covi = sqrt(1./nansum(1./si.^2,1))./mpvi; % (1,m), c.o.v.  
    stats.mpv.(tmpname{i})=mpvi;
    stats.cov.(tmpname{i})=covi;


  end

% plot results
%================
if in.flag_plt
  figure;
  nvar = length(var_name(:));
  nrow = nvar;
  ncol = 1;
  for i=1:nvar
    subplot(nrow,ncol,i);
    mpv = setup.mpv.(var_name{i});
    stdv = mpv.*setup.cov.(var_name{i});   
    errorbar(repmat(in.isetup(:),[1,m]),mpv,2*stdv,'.-');
    [dum,I] = sort(in.isetup);
    set(gca,'xtick',sort(in.isetup),'xticklabel',setup_num(in.isetup(I)));
    if isequal(var_name{i},'S')|isequal(var_name{i},'rms')
      set(gca,'yscale','log');
    end
    ylabel([var_name{i},' [',var_unit{i},']']);
    % force lower y limit to 0 if negative
    v = axis;
    if v(3)<0
      v(3)=0;
      axis(v);
    end
  end
  xlabel('Setup No.');
  aspect(1,1);
end
  
  % write results to xls
  %======================
%   tab = cell(ns+1,5);
%   titlecell = {'Setup','Freq.[Hz]','Damp.[%]',...
%     'S [(micro-g)^2/Hz]','RMS [micro-g]'}; 
%   for i=1:m
%     tab = [[1:ns].',setup.mpv.f(:,i),setup.mpv.z(:,i)*1e2,...
%       setup.mpv.S(:,i)*1e12,setup.mpv.rms(:,i)*1e6];
%     tab = mat2cell(tab,ones(ns,1),ones(1,5));
%     tab = [titlecell;tab];
%     xlswrite(in.xls_result,tab,['mode',num2str2(i,2)]);
%   end
%   disp(['->',in.xls_result]);
end % if ns>0



% set output
%===============
o.isetup = setup_num(in.isetup);
o.setup = setup;
o.stats = stats;
o.phi = phi0; % global mode shape
o.dof = dof0; % dof of global mode shape
o.c = c;
o.lam = lam;

%======================================================================
% private functions below
%======================================================================
% function x = myfun_xlsread(xls,ref)
% % similar to xlsread, but ref can contain sheet information
% I = strfind(ref,'!');
% if length(I(:))>1
%   error('cell reference name has more than one ''!''.');
% end
% if isempty(I)
%   x = xlsread(xls,ref);
% else
%   x = xlsread(xls,ref(1:I-1),ref(I+1:end));
% end
%------------------------------------------------------------------------
function y = myfun_renum(x,I)
% y = myfun_renum(x,I) gives y = same size as x with 
% y(i) = j if I(j) = x(i), 
% i.e., y is a consecutively numbered version of x
% If x(i) is not in I(:), then y(i) = nan

y = nan(size(x));
nI = length(I(:));
for i=1:nI
  y(x(:)==I(i)) = i;
end

