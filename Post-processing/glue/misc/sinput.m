function in = sinput(varargin)
% in = sinput(varargin{:})
%
% Putting the above in the beginning of a matlab function, say, FUN, has
% the following function:
%
% If FUN is called with FUN(A) when A is a structure, then IN = A
% If FUN is called with FUN(fieldname1,value1,fieldname2,value2,...) then
% IN is the structure with all the fields input
% 
% Note: FUN must be declared as 
%   function FUN(varargin)

% Written by SK Au, University of Liverpool
% 130519 - original code

nin = nargin;

if nin==0
  in = [];
elseif nin==1 % a structure is given
  in = varargin{1}; % structure
  if ~isstruct(in)
    error('A single input argument must be a structure');
  end
else % fields and values are given
  if rem(nargin,2)
    error('Field name and field values must be given');
  end
  nfield = nargin/2;
  in = [];
  for i=1:nfield
    in = setfield(in,varargin{2*i-1},varargin{2*i});
  end
end
