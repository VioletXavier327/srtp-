% Mode shape assembly and plotting toolbox
% Version 120424, 100626
%==========================================
% Main functions
%----------------
% gluemode      - glue mode shapes by least square
% glueref       - glue mode shapes using conventional method with a given reference setup
% viewsetup	    - assembly and plot mode shapes, interface provided for file input
% rptsetup      - summarize results in differnt setups
% chksetup      - graphical checking of channel/dof in setups
%
% Other supporting functions
%---------------------------
% minquadc	    - constrained eigenvalue problem
% mvmode       - show mode shape movie; called by gui_pltmode3
% refdof        - determine reference dofs 
% xlstag        - read range of xls file identified by a tag
%
% Overloaded methods
%--------------------
% gui_pltmode3, pltmode3, mvmode3_solo
%
% Notes on setup (mat) file
%===========================
% setup mat file should contain
% - mat = cell of mode id result file names corresponding to setup 1,2,...
% - xls = char of excel file name containing loc, dof, line info.
% - cellref = structure with fields loc, dof, line giving cell range of 
%             location, dof and line info
%   dof = char giving the range of dof/channel data in different setups
%         (referred as DOF data here)
%   loc = char giving the range of location data (see notes on xls info)
%         (referred as LOC data here)
%   line = char giving the range of lines (start, end)
%         (referred as LINE data here)
%
% Notes on setup (xls) file:
%============================ 
% Quick notes
%-------------
% 1. Cell address can contain sheet information, e.g., 'mysheet!E12:G21'.
%    Region covered can include blank columns and blank rows.
% 2. dofs and node numbers should be integers (+ or - ok)
%    They need not be consecutive.
%    (in principle any real number can be used but this is not recommended
%    to avoid inaccuracy in comparison of numbering during algorithm
% 3. In LOC data, nodes no. should be uniquely defined; dofs can be
%    repeated, e.g., different nodes can share some dofs in common. 
% 4. dofs defined in DOF data but not in LOC will be ignored in plotting 
%   (but their mode shape will still be assembled)
% 5. Unmeasured dofs of a measured node are set to zero in plotting 
%    mode shape
% 6. Zero '0' should not be used for dof numbering. It is reserved to
%    indicate that a dof is constrained to zero. This is needed to
%    distinguish a node that is not measured but not necessarily fixed (and
%    so will not be plotted in VIEWSETUP) and a node that is not measured
%    but is fixed (which will be plotted in VIEWSETUP). 
% 6. Lines involving nodes not defined in LOC data will be ignored
% 
% Measured DOF information
%--------------------------
% 1. The cell address should be given for the rectangular region covering
%    the measured dofs in different setup. Exclude the column for the 
%    setup numbers and the row for the channel numbers. These are assumed
%    to be numbered consecutively as 1,2,...
% 2. In each row, give the measured dof number under each channel. 
% 3. Normally the measurd dofs are defined in location data. If they are
%    not defined, the mode shape information will be ignored in VIEWSETUP
% 4. Different setups can have different no. of channels. 
%    If a channel is absent in a setup, leave the entry blank (or any char)
%    Correspondingly, the mode id result file for the setup should not 
%    contain the dof (e.g., by setting 'ichan' in the input of MODEIDFFT02)
%
% Location information 
%-----------------------
% 1. To define a node the following information must be given:
%      node no., x-, y-, z-coord.
%    Any row without full information of the above will be ignored.
%    The cell address needs to cover at least these columns.
% 2. To define dofs for the node, provide the next 3 columns:
%      node no., x-, y-, z-coord., x-, y-, z-dof, 
% 3. To define local orientation of the node (sensor), provide the next 2
%    columns:
%      node no., x-, y-, z-coord., x-, y-, z-dof, alpha, beta
%    The (local) orientation of the senor at the node is obtained by first
%    rotating the sensor placed in the global axis orientation
%    anti-clockwise about the z-axis by an angle alpha, then
%    anti-clockwise about the sensor (local) x-axis by an angle beta. Both
%    alpha and beta are measured in degree.
% 4. If all dofs of a node are not defined (blank entries) then no mode 
%    mode shape value will be plotted at the node. 
%    If only some of the dofs are defined and measured, the remaining ones
%    that are not defined or measured will be set to zero.
% 
% Line information
%-------------------
% 1. To define a line, provide in a row the starting and ending node no.
% 2. Blank rows (or char) will be ignored
% 3. Lines with node number not defined in LOC will be ignored


