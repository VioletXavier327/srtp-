function o = chksetup(varargin)
% CHKSETUP 
% reads setup information and produces plots of measured dofs in different
% setups for checking.
%
% O = CHKSETUP(IN) supplies the structure IN with the following optional fields
%  xls = char of setup (mat) file name
%  isetup = array of setup no. to process. [] do nothing
%           Note that it is the number appearing in xls file, i.e., need
%           not be consecutively numbered
%
%  The following fields are related to appearance of plots
%  flag_visible = 1 to show axis; 0 not to
%  view = 'xyz','xy','xz','yz' or (1,3) array
%  linewidth = line width for mesh
%  markersize_node = size of dots indicating locations
%  markersize_arrow = size of dot representing dof
%  stretch = scalar multiplier to adjust magnitude of deformed shape
%  flag_movie = 1 to produce movie; 0 to provide static figures
%               avi file can be produced from o.M, e.g.,
%               movie2avi(o.M,'mymovie')
%  flag_avi = 1 to produce avi file
%  avi = char of avi file name
%  fps = frames per second, used only when flag_movie = 1
%  color_ch = char or cell array of color for arrows showing channels
%             char => same for all channels; if cell dim. is less than no.
%             of channels, it will be replicated. 
%             Remember to put {color_ch} so that the input structure has a 
%             cell as its field rather than being an array of structures
%  
%  color_rep = no. of row-wise replications of color_ch
%            E.g., suppose there are 3 biaxial sensors, with channels
%            x1,y1,x2,y2,x3,3. Then setting
%            color_rep = 2 will show the two dofs measured by sensor 1 in 
%            blue, sensor 2 in green and senosr 3 in red.
%  arrow_width = line width of arrow showing measured dofs
%  aspect = cell, input to aspect(:)
%  flag_ch_txt = 1 to show channel no.; 0 not to
%  flag_node_txt = 1 to show node number next to node; 0 not to
%  flag_setup_txt = 1 to show setup number next to node measured; 0 not to
%  flag_holdon = 1 to accumulate setups (hold on); 0 not to (hold off)
%                if flag_movie = 0, this also means that all setups will be
%                shown in a single setup
%  flag_tick = 1 to show ticks on axis; 0 not to
%  fontsize_ch = font size of ch. no.
%  fontsize_node = fontsize of node number
%  fontsize_setup = fontsize of setup text
%
% By default,
%  xls = [], isetup = [], flag_visible = 1, flag_ch_txt = 1;
%  view = 'xyz', linewidth = 1, markersize_node = 6, markersize_arrow=6
%  stretch=1, flag_movie = 1
%  fps = 1, color_ch = {'b','g','r','c','m','k'}, color_rep = 1
%  arrow_width = 2, flag_avi = 0, avi = xls, aspect=[]
%  flag_node_txt = 0, flag_setup_txt = 0, flag_holdon = 0,
%  fontsize_node = 8, fontsize_setup = 8, flag_tick = 1
%
% Examples:
% 1. To check geometry and dof/channel correspondence, simply
% >> chksetup
% 2. To plot layout plan with node no. and setup no. to use on site
% >> in.flag_movie=0, in.flag_holdon=1,in.flag_node_txt=1,
%    in.flag_setup_txt=1, chksetup(in); 

% Written by SK Au, TCU
% 150914 - added color_undeform
% 120510 - added misc. to show setup no. and node no.
% 120426 - original code adapted from MVMODE

in = sinput(varargin{:});

% defaults and contants
largenumber = 999;
tag.setup = {'[SETUP]','[ENDSETUP]'};
tag.loc = {'[LOC]','[ENDLOC]'};
tag.line = {'[LINE]','[ENDLINE]'};

in_default.xls = [];
in_default.flag_visible = 1;
in_default.view = 'xyz'; % default isometric view
in_default.linewidth = 1;
in_default.markersize_node = 6;
in_default.markersize_arrow = 6;
in_default.stretch =1;
in_default.fontsize_ch = 10;
in_default.flag_line = 1;
in_default.flag_movie = 1;
in_default.fps = 1;
in_default.color_ch = {'b','g','r','c','m','k'};
in_default.color_rep = 1;
in_default.arrow_width = 2;
in_default.flag_avi = 0;
in_default.avi = [];
in_default.aspect = [];
in_default.flag_ch_txt = 1;
in_default.flag_node_txt = 0;
in_default.flag_setup_txt = 0;
in_default.flag_holdon = 0;
in_default.fontsize_node = 8;
in_default.fontsize_setup = 8;
in_default.flag_tick = 1;
in_default.color_undeform = 'r';

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

if isempty(in.avi)
  [tmpp,tmpf] = fileparts(in.xls);
  in.avi = [fullfile(tmpp,tmpf),'.avi'];
end

% line data
%-------------
in.lineI = xlstag(in.xls,tag.line{1},tag.line{2},'cleanrow','cleancol');
in.lineI = in.lineI.'; % each col gives a line

%=========================================================================
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
userdof = unique(loc_dof(~isnan(loc_dof(:)))); % all dof excl. nan
loc_dof = myfun_renum(loc_dof,userdof);

% add dofs to fill up those not defined to get 3-D effect
I = find(isnan(loc_dof(:)));
n = length(userdof(:)) + length(I(:));
loc_dof(I) = length(userdof(:)) + [1:length(I(:))];

% Process line data
%=====================================================
if ~isempty(in.lineI) 
  % convert to consecutive numbering, entries w/ undefined nodes are nan
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

% dof_setup, isetup
%--------------------------------------------------------------------
dof_setup = xlstag(in.xls,tag.setup{1},tag.setup{2},'cleanrow','cleancol');
in.setup_num = dof_setup(:,1);
dof_setup(:,1) = []; % can contain nan if ch. not used in a setup

if ~isfield(in,'isetup')
  in.isetup = in.setup_num(:).';
end
isetup_user = in.isetup; % stored for display only
in.isetup = myfun_renum(in.isetup,in.setup_num); % consecutively numbered
if any(isnan(in.isetup))
  disp('Warning: The following setups not defined in xls are removed:');
  disp([num2str(isetup_user(isnan(in.isetup)))]);
  in.isetup = in.isetup(~isnan(in.isetup));
end
dof_setup = dof_setup(in.isetup,:); % retain only selected setups

% consecutive numbering
dof_setup = myfun_renum(dof_setup,userdof);

[ns,nch] = size(dof_setup); % ns=#setups, nch=#channels

% if ns==0
%   return;
% end
%==========================================================================

% ns>0
%======================================================================

% fill up color for channels
if ~iscell(in.color_ch)
  in.color_ch = {in.color_ch};
end
in.color_ch = repmat(in.color_ch(:).',[in.color_rep 1]);
ncolor = length(in.color_ch(:));
if ncolor<nch
  in.color_ch = repmat(in.color_ch(:).',[1 ceil(nch/ncolor)]);
end

% convert view to numeric
if ~isnumeric(in.view)
  switch in.view
    case {'xy','yx'}
      in.view = [0 0 1];
    case {'yz','zy'},
      in.view = [1 0 0];
    case {'xz','zx'},
      in.view = [0 -1 0];
    case {'xyz','xzy','yxz','yzx','zxy','zyx'},
      in.view = [1 -1 1];
    otherwise,
      error('Unknow in.view');
  end
end

% calc. a suitable scaling for showing modeshape
%-------------------------------------------------
in.sp = myfun_calsp(loc_xyz)*in.stretch;

% calc. common axis, no plotting yet
%------------------------------------
tmpin = in;
% tmpin.loc_num = loc_num;
for i=1:ns
  % calc. common axis
  phi = nan(n,1);
  tmpin = myfun_axis(phi(loc_dof),loc_xyz,tmpin);  
  Ich = find(~isnan(dof_setup(i,:))); % valid ch.# in setup i
  for j=Ich(:).'
    tmpin.ch = nan(n,1);
    tmpin.ch(dof_setup(i,j)) = j;
    tmpin.ch = tmpin.ch(loc_dof);    
    phi = nan(n,1); phi(dof_setup(i,j)) = 1;
    phi = myfun_setzero(phi,loc_dof);
    phi = myfun_local2global(phi,loc_dof,loc_rot);       
    tmpin = myfun_axis(phi(loc_dof),loc_xyz,tmpin);
  end  
end

% keyboard

% plot results
%---------------
if in.flag_movie
  fprintf('Generating movies ... ');
end
for i=1:ns 
  
  if i==1
    tmpin.hfig = figure;
    pause(3);
  else
    if ~in.flag_holdon & ~in.flag_movie
      tmpin.hfig = figure;
      pause(3);
    end
  end

  
  % basic mesh 
  if ~(in.flag_movie & i>1)
    tmpin.flag_line = 1; tmpin.flag_undeform = 1; tmpin.flag_arrow = 0;
    tmpin.flag_node_txt = in.flag_node_txt;
    tmpin.flag_setup_txt = 0;
    %
    tmpin = myfun_plt(nan(size(loc_xyz)),loc_xyz,loc_num,tmpin);
    tmpin.h = []; % reset so that it will be erased later 

    if ~isempty(in.aspect)
      aspect(in.aspect{:});
    end
    
    M(1:2) = getframe(gcf);
  
  end

  hold on;
  
  h = title(['Setup ',num2str(in.setup_num(in.isetup(i)))]);
  set(h,'fontweight','bold');
  
  Ich = find(~isnan(dof_setup(i,:))); % valid ch.# in setup i
 
  p = []; xyz = []; num = []; 
  tmpin.ch = []; tmpin.arrow_color = [];
  for j=Ich(:).'
    [iloc,idof] = find(loc_dof==dof_setup(i,j)); % node # of ch.
    for k=1:length(iloc(:)) % a ch. may be shared by multiple dofs
      xyz = [xyz;loc_xyz(iloc(k),:)]; % [??,3]
      tmp = nan(1,3); tmp(idof(k)) = j;    
      tmpin.ch = [tmpin.ch; tmp];

      % p
      ptmp = zeros(1,3); ptmp(idof(k)) = 1;          
      ptmp = myfun_local2global(ptmp(:),[1:3],loc_rot(iloc(k),:));
      p = [p; ptmp(:).']; % (??,3)      
      
      % num
      num = [num;loc_num(iloc(k))];
      
      tmpin.arrow_color = [tmpin.arrow_color,tmpin.color_ch(j)];
      
    end        
  end 
  tmpin.flag_line = 0; tmpin.flag_undeform = 0; tmpin.flag_arrow = 1;
  tmpin.flag_node_txt = 0; 
  tmpin.flag_setup_txt = in.flag_setup_txt;
  tmpin.setupnum = in.setup_num(in.isetup(i));
  %
  tmpin = myfun_plt(p,xyz,num,tmpin);    

  if in.flag_movie
    M(end+1) = getframe(gcf);
    if ~in.flag_holdon
      delete(tmpin.h); tmpin.h = [];
    end
  end

  
%   rotate3d on;
  
end

if in.flag_movie
  fprintf('done\n');
end

if in.flag_holdon
  title(''); % remove title since it is no longer current
end

% play movie
%---------------
if in.flag_movie
  fprintf('Movie being played.\n');
  movie(gcf,M,largenumber,in.fps);
end


% produce avi file
%------------------
if in.flag_avi
  fprintf('Generating avi file ... ');
  movie2avi(M,in.avi,'fps',in.fps,'quality',100);
  fprintf('done => %s\n',strrep(in.avi,'\','/'));
end


% set output
if in.flag_movie && ns>0
  o.M = M;
else
  o.M = [];
end

%========================================================================
% private function
%========================================================================
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
  if all(~isnan(x(loc_dof(i,:),:)))
    a1 = loc_rot(i,1)*pi/180; % strike angle in rad
    a2 = loc_rot(i,2)*pi/180; % dip angle in rad
    R = [cos(a1) -sin(a1) 0;sin(a1) cos(a1) 0;0 0 1]; % (3,3), rot. about z
    R = R*[1 0 0; 0 cos(a2) -sin(a2);0 sin(a2) cos(a2)]; % rot. about x'
    y(loc_dof(i,:),:) = R*x(loc_dof(i,:),:); % global orient.
  end
end
%-------------------------------------------------------------------------
function y = myfun_setzero(x,loc_dof)
% set remaining dofs to zero if at least one dof is measured
y = x;
if isempty(x)
  return;
end

nloc = size(loc_dof);
for i=1:nloc
  if any(~isnan(x(loc_dof(i,:),1)))
    J = find(isnan(x(loc_dof(i,:),1)));    
    y(loc_dof(i,J),:) = 0;
  end
end

%-------------------------------------------------------------------------
function o = myfun_axis(uvw,xyz,in)
% calc. axis

uvw = uvw*in.sp;
r = 1.2; % factor for margin

% calc. axis limits
if ~isfield(in,'axisv')
  in.axisv = zeros(1,6);

  in.axisv(1) = min(xyz(:,1));
  in.axisv(2) = max(xyz(:,1));  
  in.axisv(3) = min(xyz(:,2));
  in.axisv(4) = max(xyz(:,2));  
  in.axisv(5) = min(xyz(:,3));
  in.axisv(6) = max(xyz(:,3));  
end

  in.axisv(1) = min(in.axisv(1),min(xyz(:,1)+r*uvw(:,1))); 
  in.axisv(2) = max(in.axisv(2),max(xyz(:,1)+r*uvw(:,1)));   
  in.axisv(3) = min(in.axisv(3),min(xyz(:,2)+r*uvw(:,2))); 
  in.axisv(4) = max(in.axisv(4),max(xyz(:,2)+r*uvw(:,2)));   
  in.axisv(5) = min(in.axisv(5),min(xyz(:,3)+r*uvw(:,3))); 
  in.axisv(6) = max(in.axisv(6),max(xyz(:,3)+r*uvw(:,3)));   
  
  in.axisv(1) = min(in.axisv(1),min(xyz(:,1)-r*uvw(:,1))); 
  in.axisv(2) = max(in.axisv(2),max(xyz(:,1)-r*uvw(:,1)));   
  in.axisv(3) = min(in.axisv(3),min(xyz(:,2)-r*uvw(:,2))); 
  in.axisv(4) = max(in.axisv(4),max(xyz(:,2)-r*uvw(:,2)));   
  in.axisv(5) = min(in.axisv(5),min(xyz(:,3)-r*uvw(:,3))); 
  in.axisv(6) = max(in.axisv(6),max(xyz(:,3)-r*uvw(:,3)));   

% make sure axis lim is not too short
lmax = in.axisv(2)-in.axisv(1);
lmax = max(lmax,in.axisv(4)-in.axisv(3));
lmax = max(lmax,in.axisv(6)-in.axisv(5));
%
in.axisv(2) = max(in.axisv(2),in.axisv(1)+lmax/2);
in.axisv(4) = max(in.axisv(4),in.axisv(3)+lmax/2);
in.axisv(6) = max(in.axisv(6),in.axisv(5)+lmax/2);

o = in;
%------------------------------------------------------------------------
function sp = myfun_calsp(loc_xyz)
% calc. a suitable scaling for mode shape to show

sp = zeros(1,3);

for i=1:3
  tmp = kron(loc_xyz(:,i),ones(size(loc_xyz(:,i)))) ...
        - kron(ones(size(loc_xyz(:,i))),loc_xyz(:,i));
  sp(i) = max(abs(tmp(:)));
end
sp = geomean(sp(sp~=0))/10; % this is different from mvmode

%------------------------------------------------------------------------
function o = myfun_plt(uvw,xyz,num,in)
% this function is different from that used in MVMODE

o = in;
if ~isfield(o,'h')
  o.h = []; % will accumulate
end

figure(in.hfig);

uvw = uvw*in.sp;

    % plot original loc.
    if in.flag_undeform
      h = plot3(xyz(:,1),xyz(:,2),xyz(:,3),[in.color_undeform,'.']);
      set(h,'markersize',in.markersize_node);
      o.h = [o.h(:);h(:)];
    end

    % text for nodes
    if in.flag_node_txt
      h = text(xyz(:,1),xyz(:,2),xyz(:,3),num2str(num(:)));
      set(h,'fontsize',in.fontsize_node);
      set(h,'horizontalalignment','left','verticalalignment','bottom');
      o.h = [o.h(:);h(:)];
    end
    
    % plot arrow
    if in.flag_arrow
      % add corners to have better scaling of arrows
      xyz1 = [nan(1,3);in.axisv([1 3 5;2 4 6])];
      uvw1 = [nan(1,3);nan(2,3)];
      for k=1:size(xyz,1)
        xyz1(1,:) = xyz(k,:);
        uvw1(1,:) = uvw(k,:);
        h = quiver3(xyz1(:,1),xyz1(:,2),xyz1(:,3),uvw1(:,1),uvw1(:,2),uvw1(:,3),0);
        set(h,'linewidth',in.arrow_width,'color',in.arrow_color{k});
        o.h = [o.h(:);h(:)];
        
        % marker of arrow
        h = plot3(xyz1(1,1),xyz1(1,2),xyz1(1,3),'.');
        set(h,'markersize',in.markersize_arrow,'color',in.arrow_color{k});
        o.h = [o.h(:);h(:)];
      end
    end
    
    % channel number text
    if in.flag_ch_txt
      xyz1 = xyz + 1.1*uvw;
      for k=find(~isnan(in.ch(:).'))
        [I,dum] = ind2sub(size(xyz),k); % row # of k
        h = text(xyz1(I,1),xyz1(I,2),xyz1(I,3),num2str(in.ch(k)));
        set(h,'fontweight','bold','fontsize',in.fontsize_ch);
      o.h = [o.h(:);h(:)];
      end
    end
    
    % lines
    if in.flag_line
      linex = reshape(xyz(in.lineI,1),size(in.lineI));
      liney = reshape(xyz(in.lineI,2),size(in.lineI));
      linez = reshape(xyz(in.lineI,3),size(in.lineI));     
      
      h = line(linex,liney,linez);
      set(h,'color',in.color_undeform,'linestyle','-','linewidth',in.linewidth);
      o.h = [o.h(:);h(:)];
    end

    % setup number txt
    if in.flag_setup_txt
      [dum,I]=unique(num(:));
      for i=I(:).'
        h = text(xyz(i,1),xyz(i,2),xyz(i,3),num2str(in.setupnum));
        set(h,'fontweight','bold','fontsize',in.fontsize_setup,...
              'horizontalalignment','right','verticalalignment','top',...
              'edgecolor',in.arrow_color{i},'color','black');
        o.h = [o.h(:);h(:)];
      end
    end
    
    % polish figure
    if in.flag_undeform
      axis equal;
      set(gca,'box','off');
      if ~in.flag_tick
        set(gca,'xtick',[],'ytick',[],'ztick',[]);
      end
      set(gcf,'color',[1 1 1]); % set background to white
      xlabel('x'); ylabel('y'); zlabel('z');
      axis(in.axisv);
      if in.flag_visible==0
        set(gca,'visible','off');
      end
      view(in.view);
    end  


