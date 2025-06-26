function savestruct(savestruct_fname,savestruct_S,savestruct_str)
% savestruct(FNAME,S) saves all the fields in the structure array S to 
% the file FNAME.mat. 
% 
% savestruct(FNAME,S,STR) provides additional string to the command SAVE

% Written by S. K. Au, Caltech
% 110303 - added savestruct_str
% 090923 - change names to avoid clash
% 011003 - original code

if nargin<2
  error('Both file name and structure must be given.');
end

stru2field_S = savestruct_S;
stru2field;

if nargin<3
  savestruct_str = '';
end

savestruct_names = fieldnames(savestruct_S);

save(savestruct_fname,savestruct_names{:},savestruct_str);


