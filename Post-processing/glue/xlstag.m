function [num,txt,raw]=xlstag(xls,tag1,tag2,varargin)
% [NUM,TXT,RAW] = XLSTAG(XLS,TAG1,TAG2) opens the excel file XLS and get
% get the contents in the cell range enclosed (but excluding) by the cells
% with char (tag) in TAG1 (upper left corner) and TAG2 (lower right corner)
%
% When TAG2 = [nrow,ncol], it gives the size of the cell range to the lower
% right corner of the cell of TAG1
%
% ... XLSTAG(xls,tag1,tag2,...) provides additional options
%  'cleanrow' to remove rows with all nan in NUM
%  'cleancol' to remove columns with all nan in NUM
%
% Notes:
%  1. XLS can contain sheet info, e.g., 'mysheet!B15:C23'
%  2. If there is more than one cell with the same tag, the one with the
%     smallest cell address will be used 
%  3. An error will be returned if the tag cannot be found
%
% See also xlsget, xlsput, xlsput2

% Written by SK Au, TCU
% 120427 - original

% check inputs
%---------------
if nargin<3
  error('At least 3 inputs are required.');
end

% xls - convert to cell to separate sheet info (if any)
%------------------------------------------------------
if ~ischar(xls)
  error('XLS must be a char of excel file name.');
end
I = strfind(xls,'!');
if length(I(:))>1
  error('cell reference name has more than one ''!''.');
end
if isempty(I)
  xls = {xls};
else
  xls = {xls(1:I-1),xls(I+1:end)}; % cell
end

% tag1 & tag2
%-------------
if ~ischar(tag1)
  error('TAG1 must be a char');
end
if ~ischar(tag2)
  if length(tag2(:))~=2
    error('When TAG2 is numeric it should have two entries');
  end
end

varargin = lower(varargin);

% read data
[dum,dum,raw] = xlsread(xls{:});

[I1,J1] = myfun_find(raw,tag1);

if ischar(tag2)
  [I2,J2] = myfun_find(raw,tag2);
else % numeric, 
  I2 = I1+tag2(1)+1;
  J2 = J1+tag2(2)+1;
end

if isempty(I1) | isempty(I2) % can't find tags
  [num,txt,raw] = deal([]);
  return;
end

nrow = I2-I1-1;
ncol = J2-J1-1;

if nrow<=0 | ncol<=0
  error('Tags do not enclose a valid cell range.');
end
raw = raw(I1+[1:nrow],J1+[1:ncol]); % (nrow,ncol)

% put entries with range into NUM and TXT
num = nan(nrow,ncol);
txt = cell(nrow,ncol);
[txt{:}]=deal('');
for i=1:nrow
  for j=1:ncol
    x = raw{i,j};
    if isnumeric(x)
      num(i,j) = x;
    else
      txt{i,j} = x;
    end
  end
end

% clean up output according to options
if ~isempty(intersect(varargin,{'cleanrow'}))
  num(all(isnan(num),2),:) = []; % remove blank rows
end
if ~isempty(intersect(varargin,{'cleancol'}))
  num(:,all(isnan(num),1)) = []; % remove blank rows
end


%=======================================================================
% private function below
%=======================================================================
function [I,J] = myfun_find(raw,tag)
% raw = cell, each entry can be anything
% tag = string
I = []; J = [];
[nrow,ncol]=size(raw);
for i=1:nrow
  for j=1:ncol
    if isequal(raw{i,j},tag)
      I = i;
      J = j;
      return;
    end
  end
end


