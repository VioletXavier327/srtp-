function o = viewsetup(varargin)
% VIEWSETUP combines info from multiple setups and views mode shape
% animation. 
%
% VIEWSETUP
% prompts the user for a setup mat file. Help glue/contents for the format
% of setup files
%
% O = VIEWSETUP returns a structure containing the following fields:
%   M = movie matrix which can be used for producing movie .avi file. E.g.,
%   movie2avi('mymovie',o.M) generates mymove.avi
%
% O = VIEWSETUP(IN) supplies the structure IN with the following optional fields
%  xls = char of setup (mat) file name
%  imode = array of mode indices to show; [] means all modes
%  isetup = array of setup no. to process. [] plots undeformed shape only
%           if not provided, all valid setups will be used
%           The setup no. refers to those stated in xls file
%  gluemethod = 'lsall', 'lsref' or 'ref'; for 'ref', need to provide Iref
%  flag_movie = 0 to just plot; 1 to show movie
%
% The following fields are related to appearance of plots used by MVMODE
%  hfig = figure handle; [] to generate new figure
%  flag_deform = 1 to plot deformed shape; active when given line info
%  flag_undeform = 1 to plot undeformed shape; active when given line info
%  flag_arrow = 0 to plot arrows
%  flag_visible = 1 to show axis; 0 not to
%  view = 'xyz','xy','xz','yz' or (1,3) array or 'solo'
%         If 'solo', only the first mode will be shown
%  linewidth = line width for deformed/undeformed shape
%  aspect = 'auto' aspect rato of fig. is changed automatifcally
%                to nrow by ncol, with height remaining unchanged;                
%  marksize = size of dots indicating locations
%  stretch = scalar multiplier to adjust magnitude of deformed shape
%
% By default,
%  xls = [],
%  imode = [], isetup = [], gluemethod = 'lsref', flag_movie = 1,
%  flag_deform = 1, flag_undeform = 1, flag_visible = 1,
%  view = 'xyz', linewidth = 1, aspect = '', markersize = 6, stretch=1
%
% Notes:
% 1. Only when line info is given and flag_deform = 1 will deformed shapes be
%    plotted. Otherwise arrow is plotted.
% 2. If the result file for a given setup cannot be found, the setup will
%    be automatically removed
% 3. In setup (mat) file, if mat = {} or missing, null results will be
%    returned with warning.

% Overloaded methods:
%  gui_pltmode3, pltmode3

% 121012 - swapped w/ viewmode
% 120827 - set phi(.) to 0 if loc_dof(.)=0
% 120423 - simplified code
% 120421 - substantial revision to combine solo and mult mode plots
% 120420 - substantial revision on inputs and move all procedures for
%          dof to RPTMODE; no longer allow selection of cell range
% 110309 - modeified rmssetup & rms so that it is wrt global mode shape
%          normalized to unity
% 090921 - allow local orientation to be specifed optionally
%          strike - right hand grid rule, about z
%          dip - right hand grid rule, about x'
% 091102 - allow pressing ok in channel data to assume default
% 091213 - use own scaling, allows deformed shapes to be optionally plotted
% 100112 - included solo mode for movie
% 100116 - included flag_undeform
% 100119 - included flag_aspect
% 100213 - major revision made on dof numbering
%           now numbering can be arbitary and non-consecutive
%           set of meas dof can be a subset of total no. of dof in loc
%           dof{i} can have different no. of channels
% 100214 - introduced ISETUP; required MATFILES & DOF to have the same dim. 
% 100227 - debuggged conversion of user defined dof to consecutive ones
%          debuggged conversion of user defined nodes to consecutive ones
% 100828 - allowed cell address to contain sheet name information; created
%          myfun_xlsread
% 100830 - fixed bug on measured dofs not defined in loc

in = sinput(varargin{:});

% defaults and contants
largenumber = 9999;
tag.setup = {'[SETUP]','[ENDSETUP]'};
tag.loc = {'[LOC]','[ENDLOC]'};
tag.line = {'[LINE]','[ENDLINE]'};

in_default.xls = [];
in_default.flag_movie = 1;
o = [];

if nargin<1
  in = in_default;
end
in = optfield(in,in_default);

if isempty(in.xls)
  in.xls = uigetfile2('*.xls','Pick setup xls file.');
end
if isempty(in.xls)
  return;
end
if iscell(in.xls)
  in.xls = in.xls{1}; % char
end

% line data
%-------------
in.lineI = xlstag(in.xls,tag.line{1},tag.line{2},'cleanrow','cleancol');
in.lineI = in.lineI.'; % each col gives a line

% location data
%---------------
% loc = xlstag(in.xls,tag.loc{1},tag.loc{2},'cleanrow','cleancol');
loc = xlstag(in.xls,tag.loc{1},tag.loc{2},'cleanrow');
if isempty(loc)
  error('Cannot find ''[LOC]'' tag in XLS');
end
% make sure loc has ncol_loc columns (fill with nan if not given)
if size(loc,2)<4
  error('Location data requires at least 4 col. (node #, x,y,z-coord)');
end
if size(loc,2)>9
  warning('Only the first 9 col.s of location data will be used.');
  loc = loc(:,1:9);
else
  loc = [loc,nan(size(loc,1),9-size(loc,2))];
end

% remove rows that do not supply node no., x,y,or z coord; blank rows
loc(any(isnan(loc(:,1:4)),2),:)=[];
loc_num = loc(:,1); % loc number, need not be consecutive
loc_xyz = loc(:,2:4); % (# loc, 3)
loc_dof = loc(:,5:7); % (# loc, 3)
loc_rot = loc(:,8:9); % local orientation info.
loc_rot(isnan(loc_rot)) = 0; % set all nan to zero

% convert loc_dof to consecutive numbering
loc_dof_user = loc_dof; % save for future use
userdof = unique(loc_dof(~isnan(loc_dof(:)))); % all dof excl. nan
loc_dof = myfun_renum(loc_dof,userdof);

% add dofs to fill up those not defined to get 3-D effect
I = find(isnan(loc_dof(:)));
n = length(userdof(:)) + length(I(:));
loc_dof(I) = length(userdof(:)) + [1:length(I(:))];

% Process line data
%=====================================================
if ~isempty(in.lineI)
%   % remove col. w/ nan; they may be gaps
%   in.lineI(:,any(isnan(in.lineI),1)) = []; 
  
  % convert to consecutive location numbering, entries w/ undefined nodes are nan
  lineI_user = in.lineI;
  in.lineI = myfun_renum(in.lineI,loc_num);

  I = find(isnan(in.lineI(:))); % col. index w/ undefined nodes
  if ~isempty(I)
    tmpstr = ['Warning: Removed lines with undefined nodes:'];
    disp(tmpstr); 
    disp(num2str(unique(lineI_user(I(:)).')));
  end    
  in.lineI(:,any(isnan(in.lineI),1)) = []; % remove col. containing nan
end
if isempty(in.lineI)
  in.flag_line = 0;
end


% glue mode shape
% 1. convert dof0, phi0 to consecutive numbering
% 2. remove entries w/ dof not defined in loc_dof
%=====================================================================
o_glue = rptsetup(setfield(in,'flag_plt',0));
dof0 = o_glue.dof; 
phi0 = o_glue.phi;

m = max(size(phi0,2),1); % no. of modes, min. 1

dof0_user = dof0; % store for display only
dof0 = myfun_renum(dof0,userdof); % consecutiving numbering
%
I = find(isnan(dof0));

if ~isempty(I)
  tmpstr = ['Warning: Removed measured dofs not defined ',...
    'in location data from mode shape:'];
  disp(tmpstr);
  disp(num2str(dof0_user(I(:)).'));
end
%
% remove meas. dofs not defined
phi0(isnan(dof0),:) = [];
dof0(isnan(dof0)) = [];

if isempty(phi0) & in.flag_movie
  in.flag_movie = 0;
  warning('No mode shape to plot. FLAG_MOVIE forced to 0.');
end

% Expand mode shape to cover meas. and unmeas. dofs
%======================================================
phi = nan(n,m);
phi(dof0,:) = phi0;

% at any loc, if a dof is meas, then unmeas dof are set to 0
phi = myfun_setzero(phi,loc_dof);


% set to zero if userdof(.) is zero
phi(loc_dof(loc_dof_user(:)==0),:)=0;


% transform from local to global coordinate
phi = myfun_local2global(phi,loc_dof,loc_rot);


% call MVMODE to generate frames
%======================================================================
tmpin = in; 
tmpin.f = o_glue.stats.mpv.f; tmpin.z = o_glue.stats.mpv.z;
%======================================================================
% 211128 Revised by Wei:add the c.o.v. of f and z
%======================================================================
tmpin.f_cov = o_glue.stats.cov.f;
tmpin.z_cov = o_glue.stats.cov.z;
%======================================================================
if in.flag_movie==0
  tmpin.Nf = 0.5; % 1/2 frame per half cycle,ie 1 frame per cycle
end 
M = mvmode(phi,loc_xyz,loc_dof,tmpin);


% play movie
%---------------
if in.flag_movie
  fprintf('Movie being played.\n');
  movie(gcf,M,largenumber,6);
end

% set output
%==============
o.M = M;

%========================================================================
% private function
%========================================================================
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
%-------------------------------------------------------------------------
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

%-------------------------------------------------------------------------
function y = myfun_local2global(x,loc_dof,loc_rot)
% transform from local to global coordinate
% note that z-rotation is applied first, followed by x-rotation
% (sequence matters here)

nloc = size(loc_dof,1);
y = x;
for i=1:nloc
  a1 = loc_rot(i,1)*pi/180; % strike angle in rad
  a2 = loc_rot(i,2)*pi/180; % dip angle in rad
  R = [cos(a1) -sin(a1) 0;sin(a1) cos(a1) 0;0 0 1]; % (3,3), rot. about z
  R = R*[1 0 0; 0 cos(a2) -sin(a2);0 sin(a2) cos(a2)]; % rot. about x'
  y(loc_dof(i,:),:) = R*x(loc_dof(i,:),:); % global orient.

end
%-------------------------------------------------------------------------
function y = myfun_setzero(x,loc_dof)
% set remaining dofs to zero if at least one dof is measured; 

y = x;
if isempty(x)
  return;
end

nloc = size(loc_dof);

% condition 1
for i=1:nloc
  if any(~isnan(x(loc_dof(i,:),1)))
    J = find(isnan(x(loc_dof(i,:),1)));    
    y(loc_dof(i,J),:) = 0;
  end  
end

% % condition 2
% I = loc_dof(:)==0;
% y(loc_dof(I),:) = 0;