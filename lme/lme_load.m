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
%   var_name    (Required) Variable name to load ('fr', 'psd', 'ripp', etc.).
%   var_field   (Optional) Sub-field specifier for the variable.
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
addOptional(p, 'var_name', @(x) ischar(x) || isstring(x) || iscell(x));
addOptional(p, 'var_field', '');
addOptional(p, 'vCell', {}, @iscell);
parse(p, varargin{:});

% Extract parameters
grppaths = p.Results.grppaths;
var_field = p.Results.var_field;
vCell = p.Results.vCell;

% Ensure var_name is a char
var_name = p.Results.var_name;
if iscell(var_name)
    var_name = var_name{1};
elseif isstring(var_name)
    var_name = char(var_name);
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
            vCell{igrp}{iday} = basepaths2vars('basepaths', basepaths, 'vars', {var_name});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataCell = cell(1, ngrps);

for igrp = 1:ngrps
    [~, ndays] = size(grppaths{igrp});
    data_day = cell(1, ndays);

    for iday = 1:ndays
        v = vCell{igrp}{iday};
        
        % Process data based on variable type
        if contains(var_name, 'fr')
            if contains(var_field, 'bouts')
                data_day{iday} = org_fr(v, var_field);
            else
                data_day{iday} = org_fr(v, var_field);
            end
        
        elseif contains(var_field, 'bDur')
            data_day{iday} = org_bDur(v);

        elseif contains(var_name, 'st_')
            data_day{iday} = org_brst(v, var_field);
        
        elseif contains(var_name, 'psd')
            data_day{iday} = org_band(v, var_field);
        
        elseif contains(var_name, 'f1f')
            data_day{iday} = org_f1f(v, var_field);

        elseif contains(var_name, 'units')
            data_day{iday} = org_units(v);

        elseif contains(var_name, 'ripp')
            if contains(var_field, {'fr', 'Rates'})
                data_day{iday} = org_rippSpks(v, var_name, var_field);
            else
                data_day{iday} = org_ripp(v, var_name, var_field);
            end
        end
    end

    % Concatenate days for current group
    dataCell{igrp} = cell2padmat(data_day, 2, NaN);
end

end % end function lme_load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESSING FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fr_data = org_fr(v, var_field)
% Organizes firing rate data: [mouse x day x unit x state x bout]

fr = catfields([v(:).fr], 'addim', true);
nmice = length(v);

if contains(var_field, 'states')
    % FR in states averaged across bouts
    fr_data = permute(fr.states.mfr, [3, 4, 1, 2, 5]);
    return
end

if contains(var_field, 'bouts')
    % FR in states per bout
    fr_tmp = cell(nmice, 1);
    for imouse = 1:nmice
        fr_tmp{imouse} = cell2padmat(fr.states.fr(1, :, imouse), 3);
    end
    fr_data = permute(cell2padmat(fr_tmp, 4), [4, 5, 1, 3, 2]);

else
    % FR irrespective of states
    fr_tmp = cell(nmice, 1);
    for imouse = 1:nmice
        fr_tmp{imouse} = fr.mfr(:, :, imouse);
    end
    fr_data = permute(cell2padmat(fr_tmp, 2), [2, 3, 1, 4, 5]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function brst_data = org_brst(v, brstField)
% Organizes burstiness data: [mouse x day x unit x 1 x 1]

var_name = fieldnames(v(1));
brst = catfields([v(:).(var_name{1})], 'addim', true);
brst_data = permute(brst.(brstField), [3, 4, 2, 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fooof_data = org_f1f(v, fooof_param)
% Organizes FOOOF parameters: [mouse x day x 1 x state x bout]

iBand = 1;
ap_param = 'exp';
bands_param = 'pow';

flds = fieldnames(v(1));
f1f = catfields([v(:).(flds{1})], 'addim', true, []);
% assumes struct fields organized as    [state x boutGrp x band]
% after concatenation will be           [state x boutGrp x band x mouse]
if strcmp(fooof_param, 'bands')
    fooof_data = squeeze(f1f.(fooof_param).(bands_param)(:, :, iBand, :));
elseif strcmp(fooof_param, 'ap')
    fooof_data = squeeze(f1f.(fooof_param).(ap_param)(:, :, :));
end

fooof_data = permute(fooof_data, [3, 4, 5, 1, 2]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function band_data = org_band(v, band_idx)
% Organizes band power data: [mouse x day x 1 x state x bout]

flg_norm = true;
psd = catfields([v(:).psd], 'addim', true);
nmice = length(v);
band_tmp = cell(nmice, 1);

for imouse = 1:nmice
    band_mat = cell2padmat(psd.bands.bouts(1, :, imouse), 3);
    if flg_norm
        band_tmp{imouse} = band_mat(band_idx, :, :) ./ band_mat(1, :, :);
    else
        band_tmp{imouse} = band_mat(band_idx, :, :);
    end
end

band_data = permute(cell2padmat(band_tmp, 4), [4, 5, 1, 3, 2]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bDur_data = org_bDur(v)
% Organizes bout length data: [mouse x day x 1 x state x bout]

% psd = catfields([v(:).psd], 'addim', true, [], true);
% nmice = length(v);
% bDur_tmp = cell(nmice, 1);
% for imouse = 1:nmice
% 
%     for istate = 1 : 2
%         btimes = psd.bouts.times{:, istate, imouse};
%         bb{istate} = btimes(:, 2) - btimes(:, 1);
%     end
% 
%     bDur_tmp{imouse} = cell2padmat(bb, 2);
% end
% tmp = cell2padmat(bDur_tmp, 3);
% bDur_data = permute(tmp, [3, 4, 5, 2, 1]); % [mouse x day x unit x state x bout]
% size(bDur_data)



ss = catfields([v(:).ssEmg], 'addim', true, [], true);
nmice = length(v);
bDur_tmp = cell(nmice, 1);
for imouse = 1:nmice
    bDur_tmp{imouse} = cell2padmat(ss.bouts.boutLen(:, :, imouse), 2);
end
tmp = cell2padmat(bDur_tmp, 3);
bDur_data = permute(tmp, [3, 4, 5, 2, 1]); % [mouse x day x unit x state x bout]
size(bDur_data)

% nmice = length(v);
% bDur_tmp = cell(nmice, 1);
% for imouse = 1:nmice
%     bDur_tmp{imouse} = cell2padmat(boutDur(:, imouse), 1);
% end
% 
% 
% fr = catfields([v(:).fr], 'addim', true);
% binedges = squeeze(fr.states.binedges);
% boutDur = cellfun(@(x) cellfun(@(y) diff(y), x, 'uni', true), binedges, 'uni', false);
% 
% nmice = length(v);
% bDur_tmp = cell(nmice, 1);
% for imouse = 1:nmice
%     bDur_tmp{imouse} = cell2padmat(boutDur(:, imouse), 1);
% end
% 
% padded_data = cell2padmat(bDur_tmp, 3); % [bout x state x mouse]
% bDur_data = permute(padded_data, [3, 4, 5, 1, 2]); % [mouse x day x unit x state x bout]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unit_data = org_units(v)
% Organizes unit type data: [mouse x day=1 x unit x 1 x 1]

units = catfields([v(:).units], 'addim', true);
nmice = length([v(:).units]);

unit_mouse = cell(nmice, 1);
for imouse = 1 : nmice
    clean_mouse = units.clean(:, :, imouse);        % [2 x Nunits_mouse]
    nunits_mouse = size(clean_mouse, 2);
    
    unit_type = ones(1, nunits_mouse) * 3;          % Row vector for types
    unit_type(clean_mouse(1,:)) = 1;                % Type 1 (PYR)
    unit_type(clean_mouse(2,:)) = 2;                % Type 2 (PV)
    
    unit_mouse{imouse} = reshape(unit_type, [1, 1, nunits_mouse, 1, 1]);
end

% Pad across mice and units to create a consistent matrix
unit_data = cell2padmat(unit_mouse, 1); % Pads along mouse dim (dim 1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ripp_data = org_rippSpks(v, var, var_field)
% Organizes ripple spike modulation: [mouse x day x unit x state=1 x bout=1]

normType = 'modulation';        % modulation, zscore

nMice = length(v);
spks = cell(nMice, 1);
for iMouse = 1 : nMice
    spksAvg = mean(v(iMouse).(var).spks.su.rippRates, 2, 'omitnan');
    ctrlAvg = mean(v(iMouse).(var).spks.su.ctrlRates, 2, 'omitnan');
    ctrlVar = std(v(iMouse).(var).spks.su.ctrlRates, [], 2, 'omitnan');
    if contains(var_field, 'normRates')
        switch normType
            case 'zscore'
                spksNorm = (spksAvg - ctrlAvg) ./ ctrlVar;
            case 'modulation'
                spksNorm = (spksAvg - ctrlAvg) ./ (spksAvg + ctrlAvg);
        end
        spks{iMouse} = spksNorm;
    elseif contains(var_field, 'ctrlRates')
        spks{iMouse} = ctrlAvg;
    elseif contains(var_field, 'rippRates')
        spks{iMouse} = spksAvg;
    else
        spks{iMouse} = mean(v(iMouse).(var).spks.su.(var_field), 2, 'omitnan');
    end
end
ripp_data = cell2padmat(spks, 2);
ripp_data = permute(ripp_data, [2, 3, 1, 4, 5]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ripp_data = org_ripp(v, var, var_field)
% Organizes ripple event data: [mouse x day=1 x unit=1 x state=1 x bout(=ripple)]

nmice = length(v);
nRipp = 500;
stateIdx = 4;

if strcmp(var_field, 'rate') ||  strcmp(var_field, 'density')
    rippTmp = cell(nmice, 1);
    for imouse = 1:nmice
        rippTmp{imouse} = v(imouse).(var).states.(var_field){stateIdx};
    end
    rippMat = cell2padmat(rippTmp, 2);
    ripp_data = permute(rippMat, [2, 3, 4, 5, 1]);

else

    rippMat = nan(nRipp, nmice);
    for imouse = 1:nmice
        rippTmp = v(imouse).(var).(var_field);
        idxState = v(imouse).(var).states.idx(:, stateIdx);
        rippTmp = rippTmp(idxState);

        maxRipp = length(rippTmp);
        nSlct = min(nRipp, maxRipp);
        idxRipp = randperm(maxRipp, nSlct);
        rippMat(1:nSlct, imouse) = rippTmp(idxRipp);
    end
end
ripp_data = permute(rippMat, [2, 3, 4, 5, 1]);

end