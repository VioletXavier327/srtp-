% STRU2FIELD breaks the structure 'stru2field_S' to its original fields 
% (if exist), leaving them in the current workspace. 'stru2field_S' is also 
% removed from the current work space. 
% Note that the following variables are reserved:
% stru2field_dummy, stru2field_i, stru2field_names, stru2field_nfield,
% stru2field_field stru2field_f
% 
% if the cell array 'stru2field_field' is given or is a non emtpy cell, 
% only those fields with names
% in stru2field_field will be extracted.
%  
% Example
% 
% Suppose the current work space is empty
%   S.x=1; S.y=2;
%   stru2field_S = S;
%   stru2field
% will leave S, x (=1), y (=2) in the current work space
% If in addition, 
%   stru2field_field = {'x'}
% Then y will not be in the current work space

% Siu-Kui Au, Caltech, 10-28-98
% Revision 2.0

if ~exist('stru2field_field')|isempty('stru2field_field')
  stru2field_field=fieldnames(stru2field_S);
end
stru2field_nfield=length(stru2field_field);

for stru2field_i=1:stru2field_nfield
  stru2field_f=stru2field_field{stru2field_i};
  if isfield(stru2field_S,stru2field_f)
    stru2field_tmpstr=[stru2field_f,'=getfield(stru2field_S,''',stru2field_f,''');'];
    eval(stru2field_tmpstr);
  end
end

% clean up temporary variables
clear stru2field_dummy stru2field_i stru2field_names stru2field_nfield
clear stru2field_S stru2field_field stru2field_f stru2field_tmpstr
