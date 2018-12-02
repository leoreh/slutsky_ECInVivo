function CatStruct = catstruct(varargin)

% catstruct = catstruct(varargin)
%
% goes through all folders and subfolders in parentdir, loads
% *structname*.mat if exists and concatenates the struct fields.
% for example, can be used to concatenate spikes from multiple sessions.
%
% INPUT:
%   parentdir       path to recording folder {pwd}.
%   structname      basename of saved struct {*spikes.cellinfo*}.
% 
% OUTPUT:
%   catstruct       concatenated structure 
%
% 01 dec 18 LH
% 
% to do list:
% add field indicating origins (i.e. unit x came from file y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'parentdir', pwd, @ischar);
addOptional(p, 'structname', 'spikes.cellinfo', @ischar);

parse(p, varargin{:});
parentdir = p.Results.parentdir;
structname = p.Results.structname;

% get struct files
subd = regexp(genpath(parentdir), ['[^;]*'], 'match');

j = 1;
for i = 1 : length(subd)
    ls = dir(subd{i});
    ls = {ls.name};
    sidx = contains(ls, structname);
    if any(sidx)
        sname{j} = char(fullfile(subd{i}, ls(sidx)));
        j = j + 1;
    end
end

if length(sname) <= 1
    error('there are no two structs named %s in %s', structname, parentdir);
end

% load structures
for i = 1 : length(sname)
    s{i} = load(sname{i});
end

% get structure fields
sparent = char(fieldnames(s{1}));
sfields = fieldnames(getfield(s{1}, sparent));

CatStruct = [];
for i = 1 : length(sfields)
    for j = 1 : length(sname) - 1
        CatStruct = setfield(CatStruct, sfields{i}, [getfield(s{j}, sparent, sfields{i})...
            getfield(s{j + 1}, sparent, sfields{i})]);
    end
end

end

% EOF