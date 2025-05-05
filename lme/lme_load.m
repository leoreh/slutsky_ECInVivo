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
            elseif contains(var_field, 'bDur')
                data_day{iday} = org_bDur(v);
            else
                data_day{iday} = org_fr(v, var_field);
            end
        elseif contains(var_name, 'st_')
            data_day{iday} = org_brst(v, var_field);
        elseif contains(var_name, 'psd')
            if contains(var_name, '1of')
                data_day{iday} = org_1of(v, var_field);
            else
                data_day{iday} = org_band(v, var_field);
            end
        elseif contains(var_name, 'units')
            data_day{iday} = org_units(v);

        elseif contains(var_name, 'ripp')
            if contains(var_field, 'fr')
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
function fooof_data = org_1of(v, fooof_param, iband)
% Organizes FOOOF parameters: [mouse x day x 1 x state x bout]

psd_1of = catfields([v(:).psd_1of], 'addim', true);
nmice = length(v);
fooof_tmp = cell(nmice, 1);

for imouse = 1:nmice
    if contains(fooof_param, 'ap')
        curr_param = psd_1of.(fooof_param)(:, :, imouse)';
    else
        curr_param = psd_1of.(fooof_param)(:, :, iband, imouse);
    end
    fooof_tmp{imouse} = permute(curr_param, [3, 1, 2]);
end

fooof_data = permute(cell2padmat(fooof_tmp, 4), [4, 5, 1, 2, 3]);
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

fr = catfields([v(:).fr], 'addim', true);
binedges = squeeze(fr.states.binedges);
boutDur = cellfun(@(x) cellfun(@(y) diff(y), x, 'uni', true), binedges, 'uni', false);

nmice = length(v);
bDur_tmp = cell(nmice, 1);
for imouse = 1:nmice
    bDur_tmp{imouse} = cell2padmat(boutDur(:, imouse), 1);
end

padded_data = cell2padmat(bDur_tmp, 3); % [bout x state x mouse]
bDur_data = permute(padded_data, [3, 4, 5, 1, 2]); % [mouse x day x unit x state x bout]
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

nmice = length(v);
gain_tmp = cell(nmice, 1);
max_nunits = 0;

for imouse = 1:nmice
    gain_vec = v(imouse).(var).spks.su.(var_field)(:);
    nunits_mouse = length(gain_vec);
    gain_tmp{imouse} = gain_vec';
    max_nunits = max(max_nunits, nunits_mouse);
end

ripp_data = nan(nmice, 1, max_nunits, 1, 1);
for imouse = 1:nmice
    nunits_mouse = length(gain_tmp{imouse});
    ripp_data(imouse, 1, 1:nunits_mouse, 1, 1) = gain_tmp{imouse};
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ripp_data = org_ripp(v, var, var_field)
% Organizes ripple event data: [mouse x day=1 x unit=1 x state=1 x bout(=ripple)]

nmice = length(v);
nRipp = 200;
rippMat = nan(nmice, nRipp);

for imouse = 1:nmice
    rippTmp = v(imouse).(var).(var_field);
    idxState = v(imouse).(var).states.idx(:, 4);
    rippInState = rippTmp(idxState);
    
    nAvailable = length(rippInState);
    nSamples = min(nRipp, nAvailable);
    idxSampled = randperm(nAvailable, nSamples);
    rippMat(imouse, 1:nSamples) = rippInState(idxSampled);
end

ripp_data = permute(rippMat, [1, 3, 4, 5, 2]);
end