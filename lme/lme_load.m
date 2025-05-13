function [dataCell] = lme_load(varargin)
% LME_LOAD Organizes data from multiple sessions for linear mixed effects analysis.
%
% SUMMARY:
% This function loads and organizes data from multiple recording sessions into a format
% compatible with linear mixed effects (LME) analysis. It processes data based on the
% specified variable type (e.g., firing rate, burstiness, power spectral density, ripples).
%
% INPUT:
%   grppaths    (Required) Cell (per group) of string arrays [mouse x session].
%   varName    (Required) Variable name to load ('fr', 'psd', 'ripp', etc.).
%   varField   (Optional) Sub-field specifier for the variable.
%   vCell       (Optional) Pre-loaded data as cell array of structs.
%
% OUTPUT:
%   dataCell    5D numeric matrix [mouse x day x unit x state x bout].
%
% DEPENDENCIES:
%   basepaths2vars, catfields, cell2padmat
%
% HISTORY:
%   06 Jan 24 - Initial version
%   23 Jul 24 - Changed unit handling
%   31 Jul 24 - Refactored to return 5D matrix based on 'var' input
%   01 Aug 24 - Simplified: Removed error handling, checks, direct var loading
%   05 May 25 - Refactored: Improved organization and documentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define input parameters and parse arguments
p = inputParser;
addOptional(p, 'grppaths', @iscell);
addOptional(p, 'varName', @(x) ischar(x) || isstring(x) || iscell(x));
addOptional(p, 'varField', '');
addOptional(p, 'vCell', {}, @iscell);
parse(p, varargin{:});

% Extract parameters
grppaths = p.Results.grppaths;
varField = p.Results.varField;
vCell = p.Results.vCell;

% Ensure varName is a char
varName = p.Results.varName;
if iscell(varName)
    varName = varName{1};
elseif isstring(varName)
    varName = char(varName);
end

ngrps = length(grppaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(vCell)
    vCell = cell(1, ngrps);
    for igrp = 1:ngrps
        [~, ndays] = size(grppaths{igrp});
        vCell{igrp} = cell(1, ndays);
        
        for iday = 1:ndays
            basepaths = grppaths{igrp}(:, iday);
            vCell{igrp}{iday} = basepaths2vars('basepaths', basepaths, 'vars', {varName});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataCell = cell(1, ngrps);

for igrp = 1:ngrps
    [~, ndays] = size(grppaths{igrp});
    dataDay = cell(1, ndays);

    for iday = 1:ndays
        v = vCell{igrp}{iday};
        
        % Process data based on variable type
        if contains(varName, 'fr')
            dataDay{iday} = org_fr(v, varField);
        
        elseif contains(varField, 'bDur')
            dataDay{iday} = org_bDur(v);

        elseif contains(varName, 'st_')
            dataDay{iday} = org_brst(v, varField);
        
        elseif contains(varName, 'psd')
            dataDay{iday} = org_band(v, varField);
        
        elseif contains(varName, 'f1f')
            dataDay{iday} = org_f1f(v, varField);

        elseif contains(varName, 'units')
            dataDay{iday} = org_units(v);

        elseif contains(varName, 'ripp')
            if contains(varField, {'fr', 'Rates'})
                dataDay{iday} = org_rippSpks(v, varName, varField);
            else
                dataDay{iday} = org_ripp(v, varName, varField);
            end
        end
    end

    % Concatenate days for current group
    dataCell{igrp} = cell2padmat(dataDay, 2, NaN);
end

end % end function lme_load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESSING FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function frData = org_fr(v, varField)
% Organizes firing rate data: [mouse x day x unit x state x bout]

fr = catfields([v(:).fr], 'addim', true);
nmice = length(v);

if contains(varField, 'states')
    % FR in states averaged across bouts
    frData = permute(fr.states.mfr, [3, 4, 1, 2, 5]);
    return
end

if contains(varField, 'bouts')
    % FR in states per bout
    frTmp = cell(nmice, 1);
    for imouse = 1:nmice
        frTmp{imouse} = cell2padmat(fr.states.fr(1, :, imouse), 3);
    end
    frData = permute(cell2padmat(frTmp, 4), [4, 5, 1, 3, 2]);

else
    % FR irrespective of states
    frTmp = cell(nmice, 1);
    for imouse = 1:nmice
        frTmp{imouse} = fr.mfr(:, :, imouse);
    end
    frData = permute(cell2padmat(frTmp, 2), [2, 3, 1, 4, 5]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function brstData = org_brst(v, brstField)
% Organizes burstiness data: [mouse x day x unit x 1 x 1]

varNameLocal = fieldnames(v(1)); % Renamed to avoid conflict with outer varName
brst = catfields([v(:).(varNameLocal{1})], 'addim', true);
brstData = permute(brst.(brstField), [3, 4, 2, 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fooofData = org_f1f(v, fooofParam)
% Organizes FOOOF parameters: [mouse x day x 1 x state x bout]

iBand = 1;
apParam = 'exp';
bandsParam = 'pow';

flds = fieldnames(v(1));
f1f = catfields([v(:).(flds{1})], 'addim', true, []);
% assumes struct fields organized as    [state x boutGrp x band]
% after concatenation will be           [state x boutGrp x band x mouse]
if strcmp(fooofParam, 'bands')
    fooofData = squeeze(f1f.(fooofParam).(bandsParam)(:, :, iBand, :));
elseif strcmp(fooofParam, 'ap')
    fooofData = squeeze(f1f.(fooofParam).(apParam)(:, :, :));
end

fooofData = permute(fooofData, [3, 4, 5, 1, 2]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bandData = org_band(v, bandIdx)
% Organizes band power data: [mouse x day x 1 x state x bout]

flgNorm = true;
psd = catfields([v(:).psd], 'addim', true);
nmice = length(v);
bandTmp = cell(nmice, 1);

for imouse = 1:nmice
    bandMat = cell2padmat(psd.bands.bouts(1, :, imouse), 3);
    if flgNorm
        bandTmp{imouse} = bandMat(bandIdx, :, :) ./ bandMat(1, :, :);
    else
        bandTmp{imouse} = bandMat(bandIdx, :, :);
    end
end

bandData = permute(cell2padmat(bandTmp, 4), [4, 5, 1, 3, 2]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bDurData = org_bDur(v)
% Organizes bout length data: [mouse x day x 1 x state x bout]

% psd = catfields([v(:).psd], 'addim', true, [], true);
% nmice = length(v);
% bDurTmp = cell(nmice, 1);
% for imouse = 1:nmice
% 
%     for istate = 1 : 2
%         btimes = psd.bouts.times{:, istate, imouse};
%         bb{istate} = btimes(:, 2) - btimes(:, 1);
%     end
% 
%     bDurTmp{imouse} = cell2padmat(bb, 2);
% end
% tmp = cell2padmat(bDurTmp, 3);
% bDurData = permute(tmp, [3, 4, 5, 2, 1]); % [mouse x day x unit x state x bout]
% size(bDurData)



ss = catfields([v(:).ssEmg], 'addim', true, [], true);
nmice = length(v);
bDurTmp = cell(nmice, 1);
for imouse = 1:nmice
    bDurTmp{imouse} = cell2padmat(ss.bouts.boutLen(:, :, imouse), 2);
end
tmp = cell2padmat(bDurTmp, 3);
bDurData = permute(tmp, [3, 4, 5, 2, 1]); % [mouse x day x unit x state x bout]
size(bDurData)

% nmice = length(v);
% bDurTmp = cell(nmice, 1);
% for imouse = 1:nmice
%     bDurTmp{imouse} = cell2padmat(boutDur(:, imouse), 1);
% end
% 
% 
% fr = catfields([v(:).fr], 'addim', true);
% binedges = squeeze(fr.states.binedges);
% boutDur = cellfun(@(x) cellfun(@(y) diff(y), x, 'uni', true), binedges, 'uni', false);
% 
% nmice = length(v);
% bDurTmp = cell(nmice, 1);
% for imouse = 1:nmice
%     bDurTmp{imouse} = cell2padmat(boutDur(:, imouse), 1);
% end
% 
% padded_data = cell2padmat(bDurTmp, 3); % [bout x state x mouse]
% bDurData = permute(padded_data, [3, 4, 5, 1, 2]); % [mouse x day x unit x state x bout]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unitData = org_units(v)
% Organizes unit type data: [mouse x day=1 x unit x 1 x 1]

units = catfields([v(:).units], 'addim', true);
nmice = length([v(:).units]);

unitMouse = cell(nmice, 1);
for imouse = 1 : nmice
    cleanMouse = units.clean(:, :, imouse);        % [2 x Nunits_mouse]
    nunitsMouse = size(cleanMouse, 2);
    
    unitType = ones(1, nunitsMouse) * 3;          % Row vector for types
    unitType(cleanMouse(1,:)) = 1;                % Type 1 (PYR)
    unitType(cleanMouse(2,:)) = 2;                % Type 2 (PV)
    
    unitMouse{imouse} = reshape(unitType, [1, 1, nunitsMouse, 1, 1]);
end

% Pad across mice and units to create a consistent matrix
unitData = cell2padmat(unitMouse, 1); % Pads along mouse dim (dim 1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rippData = org_rippSpks(v, varName, varField)
% Organizes ripple spike modulation: [mouse x day x unit x state=1 x bout=1]

normType = 'modulation';        % modulation, zscore

nMice = length(v);
spks = cell(nMice, 1);
for iMouse = 1 : nMice
    spksAvg = mean(v(iMouse).(varName).spks.su.rippRates, 2, 'omitnan');
    ctrlAvg = mean(v(iMouse).(varName).spks.su.ctrlRates, 2, 'omitnan');
    ctrlVar = std(v(iMouse).(varName).spks.su.ctrlRates, [], 2, 'omitnan');
    if contains(varField, 'normRates')
        switch normType
            case 'zscore'
                spksNorm = (spksAvg - ctrlAvg) ./ ctrlVar;
            case 'modulation'
                spksNorm = (spksAvg - ctrlAvg) ./ (spksAvg + ctrlAvg);
        end
        spks{iMouse} = spksNorm;
    elseif contains(varField, 'ctrlRates')
        spks{iMouse} = ctrlAvg;
    elseif contains(varField, 'rippRates')
        spks{iMouse} = spksAvg;
    else
        spks{iMouse} = mean(v(iMouse).(varName).spks.su.(varField), 2, 'omitnan');
    end
end
rippData = cell2padmat(spks, 2);
rippData = permute(rippData, [2, 3, 1, 4, 5]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rippData = org_ripp(v, varName, varField)
% Organizes ripple event data: [mouse x day=1 x unit=1 x state=1 x bout(=ripple)]

nmice = length(v);
nRipp = 500;
stateIdx = 4;

if strcmp(varField, 'rate') ||  strcmp(varField, 'density')
    rippTmp = cell(nmice, 1);
    for imouse = 1:nmice
        rippTmp{imouse} = v(imouse).(varName).states.(varField){stateIdx};
    end
    rippMat = cell2padmat(rippTmp, 2);
    rippData = permute(rippMat, [2, 3, 4, 5, 1]);

else

    rippMat = nan(nRipp, nmice);
    for imouse = 1:nmice
        rippTmp = v(imouse).(varName).(varField);
        idxState = v(imouse).(varName).states.idx(:, stateIdx);
        rippTmp = rippTmp(idxState);

        maxRipp = length(rippTmp);
        nSlct = min(nRipp, maxRipp);
        idxRipp = randperm(maxRipp, nSlct);
        rippMat(1:nSlct, imouse) = rippTmp(idxRipp);
    end
end
rippData = permute(rippMat, [2, 3, 4, 5, 1]);

end