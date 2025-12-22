function utypes_push(basepaths, fetTbl)
% UTYPES_PUSH Saves classification results to units.mat files.
%
%   utypes_push(basepaths, fetTbl)
%
%   Moves existing units.mat to 'bkup' folder and creates a new units.mat
%   containing only classification fields:
%       - clean: (2 x nUnits) logical array (Row 1: RS, Row 2: FS)
%       - type:  (nUnits x 1) double array (0: Other, 1: RS, 2: FS)
%       - date:  datetime of creation
%
%   INPUT:
%       basepaths - Cell array of recording folder paths
%       fetTbl    - Table containing 'UnitType' and 'File' columns
%

if ~iscell(basepaths), basepaths = {basepaths}; end

for iPath = 1 : length(basepaths)
    basepath = basepaths{iPath};
    [~, basename] = fileparts(basepath);
    uFile = fullfile(basepath, [basename, '.units.mat']);

    % Find units belonging to this file
    % Assumes fetTbl has 'File' column matching basename
    fFiles = string(fetTbl.File);
    uIdx = (fFiles == string(basename));

    % Get subset of unitTypes
    subTypes = fetTbl.UnitType(uIdx);
    nUnits = length(subTypes);

    % Construct 'clean' (Row 1: RS, Row 2: FS)
    clean = false(2, nUnits);
    clean(1, :) = (subTypes == 'RS');
    clean(2, :) = (subTypes == 'FS');

    % Construct 'type' (0=Other, 1=RS, 2=FS)
    type = zeros(nUnits, 1);
    type(subTypes == 'RS') = 1;
    type(subTypes == 'FS') = 2;
    type = categorical(type, [0, 1, 2], {'Other', 'RS', 'FS'});

    % Prepare structure to save
    units = struct();
    units.clean = clean;
    units.type = type;
    units.date = datetime('now');

    % Handle Backup
    if exist(uFile, 'file')
        bkupDir = fullfile(basepath, 'bkup');
        if ~exist(bkupDir, 'dir')
            mkdir(bkupDir);
        end
        % Backup with timestamp: units_YYMMDD_hhmmss.mat
        ts = datestr(now, 'yymmdd_HHMMSS');
        movefile(uFile, fullfile(bkupDir, ['units_', ts, '.mat']));
    end

    % Save new file
    save(uFile, 'units');
    fprintf('Saved units for %s\n', basename);
end

end
