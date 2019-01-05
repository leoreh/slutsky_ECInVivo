function catstruct = catStruct(varargin)

% catstruct = catstruct(varargin)
%
% loads *structname*.mat from all folders and subfolders in parentdir and
% concatenates the struct fields. for example, can be used to concatenate
% spikes from multiple sessions.
%
% INPUT:
%   parentdir       path to recording folder {pwd}.
%   structname      basename of saved struct {*spikes.cellinfo*}.
% 
% OUTPUT:
%   catstruct       concatenated structure 
%
% 01 dec 18 LH.     Updates:
% 08 dec 17 -       handling different field classes 
% 
% TO DO LIST:
%   fix cat of mats. currently creates cell within cell (e.g. spindices)
%   add field indicating origins (i.e. unit x came from file y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'parentdir', pwd, @ischar);
addOptional(p, 'structname', 'spikes.cellinfo', @ischar);

parse(p, varargin{:});
parentdir = p.Results.parentdir;
structname = p.Results.structname;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get struct files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subd = regexp(genpath(parentdir), ['[^;]*'], 'match');

j = 1;
for i = 1 : length(subd)
    ls = dir(subd{i});
    ls = {ls.name};
    sidx = contains(ls, structname);
    for k = 1 : sum(sidx)
        sname{j} = char(fullfile(subd{i}, ls(sidx)));
        j = j + 1;
    end
end

% validate there are structures to concatenate
nfiles = length(sname);
if nfiles <= 1
    error('there are no two structs named %s in %s', structname, parentdir);
end

% load structures
for i = 1 : nfiles
    s{i} = load(sname{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concatenate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% get structure fields and initialize
sparent = char(fieldnames(s{1}));
sfields = fieldnames(getfield(s{1}, sparent));
catstruct = s{1}.(sparent);

% validate all structers have same fields and order
for i = 2 : nfiles
    s{i}.(sparent) = orderfields(s{i}.(sparent), catstruct);
end

% concatenate
for i = 1 : length(sfields)
    for j = 2 : nfiles
        if all(size(s{j}.(sparent).(sfields{i})) > 1) || ischar(s{j}.(sparent).(sfields{i}))
            if j == 2
                catstruct.(sfields{i}) = {catstruct.(sfields{i})};
            end
            catstruct.(sfields{i}) = cat(1, catstruct.(sfields{i}), {s{j}.(sparent).(sfields{i})});
        elseif isrow(s{j}.(sparent).(sfields{i}))
            catstruct.(sfields{i}) = horzcat(catstruct.(sfields{i}), s{j}.(sparent).(sfields{i}));
        else
            catstruct.(sfields{i}) = vertcat(catstruct.(sfields{i}), s{j}.(sparent).(sfields{i}));
        end
    end 
end

end

% EOF