function tbl = v2tbl(varargin)
% V2TBL Converts a struct array from basepaths2vars into a table (e.g., for fitlme).
%
% SUMMARY:
% This function creates a "tidy" data table from a struct array where each
% element typically corresponds to a recording session. It is designed to
% prepare data for linear mixed-effects modeling (fitlme).
%
% INPUT (Required):
%   v           - Struct array, output of basepaths2vars. Each element v(i) 
%                 contains data for one session/mouse.
%   varMap      - Struct defining the mapping from desired table column names 
%                 to their location in the struct v.
%
% INPUT (Name-Value):
%   idxCol      - Column index or indices to select from arranged data (default: all columns).
%   tagFiles    - Struct with fields corresponding to column names. Each field should contain
%                 a cell array of values, one per file in v. Example: tagFiles.Name = {'mouse1', 'mouse2'}.
%   tagAll      - Struct with fields corresponding to column names. Each field should contain
%                 a single value to be applied to all files. Example: tagAll.Group = 'Control'.
%
% OUTPUT:
%   tbl         - Table with all data, organized for LME analysis.
%
% EXAMPLE:
%   varMap.RecTime = 'frr.recovTime';
%   varMap.BurstMiz = 'st.mizuseki';
%   tagFiles.Name = {'mouse1', 'mouse2', 'mouse3'};
%   tagAll.Group = 'Control';
%   tagAll.Day = 'Baseline';
%   tbl = v2tbl(v, varMap, 'tagFiles', tagFiles, 'tagAll', tagAll, 'idxCol', 1:5);
%
% DEPENDENCIES:
%   None

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'v', @isstruct);
addParameter(p, 'varMap', @isstruct);
addParameter(p, 'idxCol', [], @(x) isnumeric(x) && all(x > 0));
addParameter(p, 'tagFiles', struct(), @isstruct);
addParameter(p, 'tagAll', struct(), @isstruct);

parse(p, varargin{:});

v = p.Results.v;
varMap = p.Results.varMap;
idxCol = p.Results.idxCol;
tagFiles = p.Results.tagFiles;
tagAll = p.Results.tagAll;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cell array to hold the table from each session
tblCell = {};
mapFields = fieldnames(varMap);

% Global unit counter to ensure unique UnitIDs across all potential groups
persistent uOffset;
if isempty(uOffset)
    uOffset = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS EACH SESSION/MOUSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iFile = 1:length(v)
    
    % Get all field sizes recursively from v(iFile) and set nUnits to the
    % most common size that isn't 1
    fldSz = get_fldSz(v(iFile));
    fldSz = fldSz(fldSz ~= 1);
    nUnits = mode(fldSz);
    unitIdx = (1:nUnits)';

    % Create a table for the current session and add metadata, including
    % unique unitID based on mouse and unit number within mouse
    sessionTable = table();    
    
    % Add tagAll metadata (applied to all files)
    tagAllFields = fieldnames(tagAll);
    for iTag = 1:length(tagAllFields)
        tagName = tagAllFields{iTag};
        tagValue = tagAll.(tagName);
        sessionTable.(tagName) = repmat(categorical({tagValue}), nUnits, 1);
    end
    
    % Add tagFiles metadata (per-file specific)
    tagFilesFields = fieldnames(tagFiles);
    for iTag = 1:length(tagFilesFields)
        tagName = tagFilesFields{iTag};
        if iFile <= length(tagFiles.(tagName))
            tagValue = tagFiles.(tagName){iFile};
            sessionTable.(tagName) = repmat(categorical({tagValue}), nUnits, 1);
        end
    end
    
    sessionTable.UnitID = uOffset + (iFile * 1000) + unitIdx;

    % Extract all variables based on the map
    for iVar = 1:length(mapFields)
        colName = mapFields{iVar};
        structPath = varMap.(colName);    
        pathParts = strsplit(structPath, '.');
        data = getfield(v(iFile), pathParts{:});
        
        % Ensure data has nUnits as the first dimension (rows)
        if size(data, 1) ~= nUnits
            if size(data, 2) == nUnits
                data = data'; % Transpose to make nUnits the first dimension
            else
                % If neither dimension matches, try to reshape
                data = reshape(data, nUnits, []);
            end
        end
        
        % Select specific columns if idxCol is specified
        if ~isempty(idxCol) && size(data, 2) > 1
            data = data(:, idxCol);
        end
        
        sessionTable.(colName) = data;
    end

    % Append the table for this session to our list
    tblCell{end+1} = sessionTable;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMBINE ALL TABLES AND UPDATE GLOBAL OFFSET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vertically concatenate all session tables into one big table
tbl = vertcat(tblCell{:});

% Update the global offset for the next time this function is called
[~, maxID] = max(tbl.UnitID);
uOffset = tbl.UnitID(maxID);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: GET_FLDSZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fldSz = get_fldSz(s)
% Recursively gets sizes of all fields in a struct, including nested structs
    
    fldSz = [];
    fields = fieldnames(s);
    
    for iFld = 1:length(fields)
        fieldVal = s.(fields{iFld});
        
        if isstruct(fieldVal)
            % Recursively get sizes from nested struct
            nestedSizes = get_fldSz(fieldVal);
            fldSz = [fldSz, nestedSizes];
        else
            % Get size of non-struct field
            fieldSize = size(fieldVal);
            fldSz = [fldSz, fieldSize];
        end
    end
end

% EOF