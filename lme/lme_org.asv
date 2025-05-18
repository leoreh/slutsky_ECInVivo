function [lmeData, lmeCfg] = lme_org(varargin)
% LME_ORG Organizes data from multiple sessions for linear mixed effects analysis.
%
% SUMMARY:
% This function creates a table for linear mixed effects modeling by loading
% and organizing data using lme_load. The data structure depends on the formula.
%
% INPUT:
%   grppaths    (Required) Cell (per group) of string arrays [mouse x session].
%   frml        (Required) Formula for LME analysis. Determines how the table is organized.
%   flgEmg      (Optional) Logical. Load vars as states (false) or EMG states (true).
%   varField    (Optional) Sub-field specifier for the variable.
%   vCell       (Optional) Pre-loaded data as cell array of structs.
%
% OUTPUT:
%   lmeData     Table organized for LME with fields based on formula.
%   lmeCfg      Struct with metadata and analysis parameters.
%
% DEPENDENCIES:
%   lme_load, cell2padmat
%
% HISTORY:
%   06 Jan 24 - Initial version
%   23 Jul 24 - Changed unit handling
%   05 May 25 - Restructured to work directly with lme_load
%   08 May 25 - Simplified data processing loop, standardized formatting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define input parameters and parse arguments
p = inputParser;
addOptional(p, 'grppaths', @iscell);
addOptional(p, 'frml', @(x) ischar(x) || isstring(x));
addOptional(p, 'flgEmg', false, @islogical);
addOptional(p, 'varField', '');
addOptional(p, 'vCell', {}, @iscell);
parse(p, varargin{:});

% Extract parameters
grppaths = p.Results.grppaths;
frml = p.Results.frml;
flgEmg = p.Results.flgEmg;
varField = p.Results.varField;
vCell = p.Results.vCell;

% Ensure frml is a char
if isstring(frml)
    frml = char(frml);
end

% Set labels according to specific experiment and formula
ngrps = length(grppaths);
[strGrp, strDay, strState] = org_strLbls(grppaths, frml, flgEmg);
strUnit = {'pPYR', 'pPV', 'Other'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE VARIABLES NEEDED BASED ON FORMULA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables based on formula
uCell = {};
bDurCell = {};
yName = '';
flgBDur = false;
flgUnits = false;

% Parse formula to determine needed variables and response variable name
if contains(frml, 'FR ~')
    if flgEmg
        varName = 'frEmg';
    else
        varName = 'fr';
    end
    flgUnits = true;
    yName = 'FR';
    
    if contains(frml, 'State')
        varField = 'states';
    end

elseif contains(frml, 'Burst ~')
    varName = 'st_metrics';
    varName = 'st_brst';
    flgUnits = true;
    if isempty(varField)
        varField = 'lidor';
    end
    yName = 'Burst';
    
elseif contains(frml, 'BDur ~') 
    if flgEmg
        varName = 'frEmg';
        varName = 'sleep_statesEmg';
        % varName = 'psdEmg';
    else
        varName = 'fr';
    end
    yName = 'BDur';
    varField = 'bDur';
    flgBDur = true;

elseif contains(frml, 'Band ~')
    if flgEmg
        varName = 'psdEmg';
    else
        varName = 'psd';
    end
    yName = 'Band';
    if isempty(varField)
        varField = 1;
    end
    
elseif contains(frml, 'FOOOF ~')    
    if flgEmg
        varName = 'psdEmg_1of';
        varName = 'f1f_eeg';
        % varName = 'f1f';
    else
        varName = 'psd_1of';
    end
    yName = 'FOOOF';
    if isempty(varField)
        varField = 'bands';
    end
    
elseif contains(frml, 'RippSpks ~')
    varName = 'ripp';
    flgUnits = true;
    yName = 'RippSpks';
    if isempty(varField)
        varField = 'frGain';
    end

elseif contains(frml, 'RippSpkLfp ~')
    varName = 'ripp';
    flgUnits = true;
    yName = 'RippSpkLfp';
    if isempty(varField)
        varField = 'theta';
    end

elseif contains(frml, 'Ripp ~')
    varName = 'ripp';
    yName = 'Ripp';
    if isempty(varField)
        varField = 'peakAmp';
    end
end

% Check if we need bout length as an additional variable
if contains(frml, 'BoutDur')
    varField = 'bouts';
    flgBDur = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA USING LME_LOAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load main variable data
dataCell = lme_load('grppaths', grppaths, 'varName', varName,...
    'varField', varField, 'vCell', vCell);

% Load units data if needed
if flgUnits
    uCell = lme_load('grppaths', grppaths, 'varName', 'units', 'vCell', vCell);
end

% Load bout duration data if needed, in addition to the response data
if flgBDur
    bDurCell = lme_load('grppaths', grppaths, 'varName', varName, 'varField', 'bDur', 'vCell', vCell);   
end

% Get mouse names for labeling
mnames = cell(1, ngrps);
for igrp = 1:ngrps
    mnames{igrp} = org_mnames(grppaths{igrp});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE TABLE SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrows = 0;
for igrp = 1:ngrps
    nrows = nrows + sum(~isnan(dataCell{igrp}(:))); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREALLOCATE ARRAYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vecData = nan(nrows, 1);
vecBDur = nan(nrows, 1);
lblGrp = strings(nrows, 1);
lblState = strings(nrows, 1);
lblMouse = strings(nrows, 1);
lblDay = strings(nrows, 1);
lblUnit = strings(nrows, 1);
idUnit = nan(nrows, 1);

% Initialize row counter
idxCurr = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL ARRAYS WITH DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for igrp = 1:ngrps
    [nmice, ndays, ~, nstates, nbouts] = size(dataCell{igrp});
    
    for imouse = 1:nmice
        for iday = 1:ndays
            
            % Determine valid units and their types
            if ~isempty(uCell)
                unitsType = squeeze(uCell{igrp}(imouse, iday, :, 1, 1));
            else
                unitsType = 1;  % Default to unit type 1
            end
            
            for istate = 1:nstates
                for ibout = 1:nbouts
                    
                    % Get data of response variable 
                    dataTmp = squeeze(dataCell{igrp}(imouse, iday, :, istate, ibout));
                    
                    % Identify valid (non-NaN) values
                    validMask = ~isnan(dataTmp);
                    nValid = sum(validMask);
                    
                    % index to current placement of data
                    idxRange = idxCurr:(idxCurr + nValid - 1);

                    % Store data points
                    vecData(idxRange) = dataTmp(validMask);

                    % Create corresponding unit indices and types
                    if length(unitsType) > 1
                        % For multiple units, determine which unit each value belongs to
                        [unitsIdx, ~] = ind2sub(size(dataTmp), find(validMask));
                        unitTypeVal = unitsType(unitsIdx); % Renamed to avoid conflict with outer unitsType
                    else
                        % For single unit or no units
                        unitTypeVal = repmat(unitsType, nValid, 1); % Renamed
                    end                    

                    % Get bout length if needed
                    if flgBDur
                        currBDur = bDurCell{igrp}(imouse, iday, 1, istate, ibout);
                        vecBDur(idxRange) = repmat(currBDur, nValid, 1);
                    end

                    % Assign common labels to all data points
                    lblGrp(idxRange) = repmat(strGrp{igrp}, nValid, 1);
                    lblState(idxRange) = repmat(strState{istate}, nValid, 1);
                    lblMouse(idxRange) = repmat(mnames{igrp}{imouse}, nValid, 1);
                    lblDay(idxRange) = repmat(strDay{iday}, nValid, 1);

                    % Handle unit information
                    if flgUnits
                        for i = 1:nValid
                            idxValid = find(validMask);
                            currIdx = idxValid(i);

                            lblUnit(idxCurr + i - 1) = strUnit{unitTypeVal(currIdx)}; % Used renamed unitTypeVal
                            idUnit(idxCurr + i - 1) = unitsIdx(currIdx) + ...
                                1000 * (imouse - 1) + ...
                                10000 * (iday - 1) + ...
                                100000 * (igrp - 1);
                        end
                    else
                        lblUnit(idxRange) = repmat(string(strUnit{1}), nValid, 1);
                        idUnit(idxRange) = nan(nValid, 1);
                    end

                    idxCurr = idxCurr + nValid;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE OUTPUT TABLE AND CONFIG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create table with appropriate columns
if flgBDur && ~strcmp(yName, 'BDur')
    lmeData = table(vecData, lblGrp, lblState, lblMouse, lblUnit, vecBDur, idUnit, lblDay,...
        'VariableNames', {yName, 'Group', 'State', 'Mouse', 'UnitType', 'BoutDur', 'UnitID', 'Day'});
else
    lmeData = table(vecData, lblGrp, lblState, lblMouse, lblUnit, idUnit, lblDay,...
        'VariableNames', {yName, 'Group', 'State', 'Mouse', 'UnitType', 'UnitID', 'Day'});
end

% Clean Table
% Find rows with NaN or Inf values
idxRemove = isnan(lmeData.(yName)) | isinf(lmeData.(yName));

% Find rows with 'Other' as unit type
idxRemove = idxRemove | lmeData.UnitType == 'Other';

% Remove these rows from the table and print warning
nRemoved = sum(idxRemove);
if nRemoved > 0
    lmeData(idxRemove, :) = [];
    warning('Removing %d rows containing NaN or Inf values from %s', nRemoved, yName)
end

% Convert columns to categorical
lmeData.State = categorical(lmeData.State);
lmeData.Mouse = categorical(lmeData.Mouse);
lmeData.Day = categorical(lmeData.Day);
lmeData.UnitType = categorical(lmeData.UnitType);
lmeData.Group = categorical(lmeData.Group);

% Set reference categories
% Makes Control the reference
if any(unique(lmeData.Group) == categorical({'Control'}))
    cats = categories(lmeData.Group);
    otherCat = cats(~strcmp(cats, 'Control'));
    lmeData.Group = reordercats(lmeData.Group, ['Control', otherCat]);
end

% Set specific Day order if BSL is the first day from org_strLbls
% Makes BSL the reference and assert the order in strDay
if any(unique(lmeData.Day) == categorical({'BSL'}))
    tblCats = categories(lmeData.Day);
    orderedCats = strDay(ismember(strDay, tblCats));
    lmeData.Day = reordercats(lmeData.Day, orderedCats);
end

% Makes pPYR the reference
if any(unique(lmeData.UnitType) == categorical({'pPYR'}))
    cats = categories(lmeData.UnitType);
    otherCat = cats(~strcmp(cats, 'pPYR'));
    lmeData.UnitType = reordercats(lmeData.UnitType, ['pPYR'; otherCat]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE CONFIG OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store metadata in config structure
lmeCfg.grps = grppaths;
lmeCfg.frml = frml;
lmeCfg.flgEmg = flgEmg;
lmeCfg.varField = varField;
lmeCfg.nRemoved = nRemoved;  

end % end function lme_org

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mnames = org_mnames(basepaths)
% Extracts mice names from first day's paths (assumes consistent across days)

[nmice, ~] = size(basepaths);
mnames = cell(nmice, 1);
for imouse = 1:nmice
    [~, mname] = fileparts(basepaths(imouse, 1));
    mnames{imouse} = regexp(mname, '^[^_]*', 'match', 'once');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [strGrp, strDay, strState] = org_strLbls(grppaths, frml, flgEmg)
% Creates string labels for groups, days, and states

% Group labels from paths (assumes 1st cell is Control and 2nd MCU-KO)
ngrps = length(grppaths);
if ngrps == 2
    strGrp = {'Control'; 'MCU-KO'};
else
    strGrp = split(num2str(1:ngrps));
end

% State labels
if flgEmg
    strState = {'High EMG'; 'Low EMG'};
else
    strState = {'AW'; 'NREM'; 'REM'};
end

% Day labels. Assumes 2nd dim of grppaths is session (day)
ndays = unique(cellfun(@(x) size(x, 2), grppaths, 'uni', true));
if contains(frml, ' Day')
    if ndays == 7
        strDay = {'BSL'; 'BAC On'; 'BAC1'; 'BAC2'; 'BAC3'; 'BAC Off'; 'WASH'};
    elseif ndays == 5
        strDay = {'BSL'; 'BAC1'; 'BAC2'; 'BAC3'; 'WASH'};
    end
else
    if ndays == 1
        strDay = {'BSL'};
    else
        strDay = split(num2str(1:ndays));
    end
end
end