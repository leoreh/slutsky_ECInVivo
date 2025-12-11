function uTbl = mcu_unitData()

% MCU_UNITDATA Loads and processes unit data for the MCU project.

% Load Data
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];

% Metadata
tagFiles = struct();
tagFiles.Name = get_mname(basepaths);
[~, fileNames] = fileparts(basepaths);
tagFiles.File = fileNames;

% Load
cfg = mcu_cfg;
v = basepaths2vars('basepaths', basepaths, 'vars', cfg.vars);

% Table
uTbl = v2tbl('v', v, 'varMap', cfg.varMap, 'tagAll',...
    struct(), 'tagFiles', tagFiles, 'idxCol', []);

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

% Clean up
% uTbl = rmmissing(uTbl);
% uTbl(uTbl.UnitType == 'Other', :) = [];
% uTbl.UnitType = removecats(uTbl.UnitType, 'Other');

end
