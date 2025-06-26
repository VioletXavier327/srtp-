function [fullnames,fnames,paths] = uigetfile2(varargin)
% fullfnames = uigetfile2(...)
% does a similar job as UIGETFILE except that the output fullfnames is always a
% cell array of full file names. The file names are also sorted.
%
% [fullfnames,fnames,p] = ... also returns the names only and the common
% path


[fnames,paths]=uigetfile(varargin{:});

if isequal(fnames,0)
  disp('Action cancelled.');
  fullnames = {}; 
  fnames = {};
  paths = {};
  return;
end

if ~isa(fnames,'cell')
    fnames = {fnames};
end
fnames = sort(fnames);
nfiles = length(fnames);

% convert fnames to fullfile name
fullnames = cell(1,nfiles);
for i=1:nfiles
  fullnames{i} = fullfile(paths,fnames{i});
end
