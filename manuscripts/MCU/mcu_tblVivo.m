function [uTbl, basepaths, v] = mcu_tblVivo(varargin)

% MCU_TBLVIVO Loads and processes unit data for the MCU project.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                loads defaults.
%   v            (struct) Pre-loaded data structure.
%   varMap       (struct) Variable mapping for table creation.
%   flgClean     (logical) Remove bad units and bac on / off.
%
% OUTPUT:
%   uTbl         (table) Unit table with metadata.
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'v', [], @isstruct);
addOptional(p, 'varMap', [], @isstruct);
addOptional(p, 'flgClean', false, @islogical);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
v = p.Results.v;
varMap = p.Results.varMap;
flgClean = p.Results.flgClean;

if isempty(basepaths)
    basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
end

%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Metadata
tagFiles = struct();
tagFiles.Name = get_mname(basepaths);
[~, fileNames] = fileparts(basepaths);
tagFiles.File = fileNames;

% Load
cfg = mcu_cfg;
if isempty(v)
    v = basepaths2vars('basepaths', basepaths, 'vars', cfg.vars);
end

% Set varMap
if isempty(varMap)
    varMap = cfg.varMap;
end

% Table
uTbl = v2tbl('v', v, 'varMap', varMap, 'tagAll',...
    struct(), 'tagFiles', tagFiles, 'idxCol', []);

%% ========================================================================
%  PROCESS METADATA
%  ========================================================================

% Group metadata
uTbl.Group = ones(height(uTbl), 1) * 1;
uTbl.Group(ismember(uTbl.Name, cfg.miceMCU), :) = 2;
uTbl.Group = categorical(uTbl.Group, [1, 2], cfg.lbl.grp);

% Day metadata
fileTbl = unique(uTbl(:, {'Name', 'File'}), 'rows');
fileGrp = findgroups(fileTbl.Name);
dayCell = splitapply(@(x) {(1:numel(x))'}, fileTbl.File, fileGrp);
fileTbl.Day = vertcat(dayCell{:});
uTbl = join(uTbl, fileTbl, 'Keys', {'Name', 'File'});
uTbl.Day = categorical(uTbl.Day, [1 : 7], cfg.lbl.day);

% Reorder columns
varOrder = {'Group', 'Name', 'File', 'Day', 'UnitID', 'UnitType'};
uTbl = movevars(uTbl, varOrder, 'Before', 1);

% Assert category order
uTbl.Group = reordercats(uTbl.Group, cfg.lbl.grp);
uTbl.UnitType = reordercats(uTbl.UnitType, cfg.lbl.unit);
uTbl.Day = reordercats(uTbl.Day, cfg.lbl.day);

if flgClean
    % Remove bad units
    uTbl(uTbl.UnitType == 'Other', :) = [];
    uTbl.UnitType = removecats(uTbl.UnitType, 'Other');

    % Remove bac on / off
    uTbl(uTbl.Day == 'BAC_ON', :) = [];
    uTbl(uTbl.Day == 'BAC_OFF', :) = [];
    uTbl.Day = removecats(uTbl.Day, {'BAC_ON', 'BAC_OFF'});
end

end     % EOF
