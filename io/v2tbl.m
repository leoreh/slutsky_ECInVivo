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
%   groupName   - String or char array for the 'Group' label.
%   mouseNames  - Cell array of strings with mouse names, corresponding to 
%                 each element of v.
%
% INPUT (Name-Value):
%   colIdx      - Column index or indices to select from arranged data (default: all columns).
%
% OUTPUT:
%   tbl         - Table with all data, organized for LME analysis.
%
% EXAMPLE:
%   varMap.RecTime = 'frr.recovTime';
%   varMap.BurstMiz = 'st.mizuseki';
%   tbl = v2tbl(v, varMap, 'Control', mouseNames, 'colIdx', 1:5);
%
% DEPENDENCIES:
%   None

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'v', @isstruct);
addOptional(p, 'varMap', @isstruct);
addOptional(p, 'groupName', @(x) ischar(x) || isstring(x));
addOptional(p, 'mouseNames', @iscell);
addParameter(p, 'colIdx', [], @(x) isnumeric(x) && all(x > 0));

parse(p, varargin{:});

v = p.Results.v;
varMap = p.Results.varMap;
groupName = p.Results.groupName;
mouseNames = p.Results.mouseNames;
colIdx = p.Results.colIdx;

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
    sessionTable.Group = repmat(string(groupName), nUnits, 1);
    sessionTable.Mouse = repmat(string(mouseNames{iFile}), nUnits, 1);
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
        
        % Select specific columns if colIdx is specified
        if ~isempty(colIdx) && size(data, 2) > 1
            data = data(:, colIdx);
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