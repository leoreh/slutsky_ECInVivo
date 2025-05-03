function [lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, flg_emg, var_field, var_idx, vCell)

% organizes data from multiple sessions into format compatible with lme
% analysis. loads data according to mcu_sessions and arranges in a table.
% data can be firing rate or burstiness, determined by the formula.
%
% INPUT
%   grppaths    cell (per group) of string arrays [mouse x session]
%   frml        char of fomrula for lme. determines how the table is
%               organized
%   flg_emg     logical. load vars by as states (false) or emg states {true}
%   var_field   string. variable to organize as fixed effects
%   var_idx     integer. index for band power or frequency band
%   vCell       cell (per group) of cell arrays (per day) of structs
%
% OUTPUT
%   lme_tbl    table organized for lme with fields FR, Group, State, Mouse, UnitType
%   lme_cfg    struct with metadata and analysis parameters
%
% CALLS
%   mcu_sessions
%   basepaths2vars
%   catfields
%
% 06 Jan 24
% 23 Jul 24 - Changed unit handling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    flg_emg = false;
end
if nargin < 4
    var_field = [];
end
if nargin < 5
    var_idx = [];
end
if nargin < 6
    vCell = {};
end

% set labels according to specific experiment (e.g., mcu) and formula
ngrps = length(grppaths);
[str_grp, str_day, str_state] = org_strLbls(grppaths, frml, flg_emg);
str_unit = {'pPYR', 'pPV'};
nunits = length(str_unit);

% set which variable to organize as fixed effects, according to formula
vars = {};
if contains(frml, 'FR ~')
    if flg_emg
        vars = {'frEmg'; 'units'};
    else
        vars = {'fr'; 'units'};
    end
    yName = 'FR';

elseif contains(frml, 'Burst ~')
    % vars = {'st_brst'; 'units'};
    % brstField = 'spkprct';
    % brstField = 'rateNorm';

    vars = {'st_metrics'; 'units'};
    if isempty(var_field)
        var_field = 'lidor';
    end

    yName = 'Burst';

elseif contains(frml, 'BLen ~') || contains(frml, ' BoutLength')
    % bout lengths extracted from fr struct
    if flg_emg
        vars = {'frEmg'};
    else
        vars = {'fr'};
    end
    yName = 'BLen';
end

if contains(frml, 'Band ~')
    if flg_emg
        vars = [vars, {'psdEmg'}];
    else
        vars = [vars, {'psd'}];
    end
    yName = 'Band';

    % bandNames = ["broad", "swa", "delta", "theta", "beta", "gamma"];
    if isempty(var_idx)
        var_idx = 1;
    end

elseif contains(frml, 'FOOOF ~')
    if flg_emg
        vars = [vars, {'psdEmg_1of'}];
    else
        vars = [vars, {'psd_1of'}];
    end
    yName = 'FOOOF';

    % parameter and frequency band to analyze
    if isempty(var_field)
        var_field = 'ap_exp';
    end
    if isempty(var_idx)
        var_idx = 1;
    end

elseif contains(frml, 'RippSpks ~')
    vars = {'ripp'; 'units'}; 
    yName = 'RippSpks';
    if isempty(var_field)
        var_field = 'frModulation';
    end

elseif contains(frml, 'Ripp ~')
    vars = {'ripp'}; 
    yName = 'Ripp';
    if isempty(var_field)
        var_field = 'peakAmp';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize cells for raw data
dataCell = cell(1, ngrps);
uCell = cell(1, ngrps);
blenCell = cell(1, ngrps);
mnames = cell(1,  ngrps);

for igrp = 1 : ngrps

    % extract mice names for labeling
    mnames{igrp} = org_mnames(grppaths{igrp});

    % initialize
    [~, ndays] = size(grppaths{igrp});
    data_day = cell(1, ndays);
    u_day = cell(1, ndays);
    blen_day = cell(1, ndays);

    % Loop through days and load data
    for iday = 1 : ndays

        % Load vars for each day
        if ~isempty(vCell)
            v = vCell{igrp}{iday}; % Use preloaded data
        else
            basepaths = grppaths{igrp}(:, iday);
            v = basepaths2vars('basepaths', basepaths, 'vars', vars); % Load data
        end

        % Extract data
        u_day{iday} = org_units(v);                                         % units 
        blen_day{iday} = org_blen(v, frml);                                 % bout length (fixed / random effect)

        if contains(frml, 'FR ~')
            data_day{iday} = org_fr(v, frml);                               % firing rate

        elseif contains(frml, 'Burst ~')
            data_day{iday} = org_brst(v, var_field);                        % burstiness

        elseif contains(frml, 'Band ~')
            data_day{iday} = org_band(v, var_idx);                          % band power

        elseif contains(frml, 'BLen ~')
            data_day{iday} = org_blen(v, 'BoutLength');                     % bout length (response variable)

        elseif contains(frml, 'FOOOF ~')
            data_day{iday} = org_1of(v, var_field, var_idx);                % psd fooof parameter
        
        elseif contains(frml, 'RippSpks ~') 
            data_day{iday} = org_rippSpks(v, var_field);                      % ripple FR increase
        
        elseif contains(frml, 'Ripp ~')
            data_day{iday} = org_ripp(v, var_field);                        % specified ripple property
            u_day = cell(1, ndays);
        end
    end

    % Concatenate across days
    dataCell{igrp} = cell2padmat(data_day, 2);                              % [mouse x day x unit x state x bout]
    uCell{igrp} = cell2padmat(u_day, 2);                                    % [mouse x day x unit x 1 x 1], value is type (1 or 2)
    blenCell{igrp} = cell2padmat(blen_day, 2);                              % [mouse x day x 1 x state x bout]
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize arrays for lme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate table size
nrows = 0;
for igrp = 1 : ngrps
    [nmice, ndays, max_nunits_grp, nstates, nbouts] = size(dataCell{igrp});

    if ~isempty(uCell{igrp})
        for imouse = 1 : nmice
            for iday = 1 : ndays
                % Get unit types for this mouse/day [max_nunits_grp x 1]
                unit_types_curr = squeeze(uCell{igrp}(imouse, iday, 1:max_nunits_grp, 1, 1)); 
                valid_unit_indices = find(unit_types_curr > 0); % Indices of units with assigned type

                if ~isempty(valid_unit_indices)
                    for istate = 1 : nstates
                        for ibout = 1 : nbouts
                            % Get data for valid units [n_valid_units x 1]
                            data_curr = squeeze(dataCell{igrp}(imouse, iday, valid_unit_indices, istate, ibout));
                            % Count non-NaN data points among valid units
                            nrows = nrows + sum(~isnan(data_curr));
                        end
                    end
                end
            end
        end
    else % Case for LFP data or if uCell is empty
        nrows = nrows + sum(~isnan(dataCell{igrp}(:)));
    end
end

% Preallocate arrays
vec_data = nan(nrows, 1);
vec_blen = nan(nrows, 1);
lbl_grp = strings(nrows, 1);
lbl_state = strings(nrows, 1);
lbl_mouse = strings(nrows, 1);
lbl_day = strings(nrows, 1);
lbl_unit = strings(nrows, 1);
id_unit = nan(nrows, 1);

% Initialize counter
idx_curr = 1;

% Organize data into table format
for igrp = 1 : ngrps
    [nmice, ndays, ~, nstates, nbouts] = size(dataCell{igrp});

    for imouse = 1 : nmice
        for istate = 1 : nstates
            for iday = 1 : ndays
                if ~isempty(uCell{igrp})
                    % Get unit info for this mouse/day
                    unit_types_all = squeeze(uCell{igrp}(imouse, iday, :, 1, 1)); % [max_nunits x 1]
                    unit_indices_exist = find(unit_types_all > 0); % Indices of units that exist
                    actual_unit_types = unit_types_all(unit_indices_exist); % Types (1 or 2) of existing units

                    if isempty(unit_indices_exist) % Skip if no units found for this mouse/day
                        continue;
                    end

                    for ibout = 1 : nbouts
                        % Get data for existing units [n_exist_units x 1]
                        data_exist_units = squeeze(dataCell{igrp}(imouse, iday, unit_indices_exist, istate, ibout));
                        
                        idx_valid_data = ~isnan(data_exist_units); % Mask for non-NaN data among existing units
                        nvalid = sum(idx_valid_data);

                        if nvalid > 0
                            idx_range = idx_curr : (idx_curr + nvalid - 1);

                            % Store data
                            vec_data(idx_range) = data_exist_units(idx_valid_data);

                            % Get bout length (if requested)
                            if contains(frml, ' BoutLength')
                                curr_blen = blenCell{igrp}(imouse, iday, 1, istate, ibout);
                                vec_blen(idx_range) = curr_blen;
                            end

                            % Assign labels
                            lbl_grp(idx_range) = str_grp{igrp};
                            lbl_state(idx_range) = str_state{istate};
                            lbl_mouse(idx_range) = mnames{igrp}{imouse};
                            lbl_day(idx_range) = str_day{iday};
                            
                            % Assign unit type labels based on the type value (1 or 2)
                            unit_types_valid = actual_unit_types(idx_valid_data);
                            lbl_unit(idx_range) = str_unit(unit_types_valid);

                            % Assign unique unit ID
                            unit_indices_valid = unit_indices_exist(idx_valid_data);
                            id_unit(idx_range) = unit_indices_valid + ...
                                1000 * (imouse - 1) + ...
                                10000 * (iday - 1) + ...
                                100000 * (igrp - 1);
                            
                            % Update counter
                            idx_curr = idx_curr + nvalid;
                        end
                    end % end ibout loop
                else % Handle non-unit specific data (LFP, etc.) or case where flg_hasUnits is false
                    nunits_loop = 1; % Treat as single "unit"
                    for iunit_dummy = 1 : nunits_loop % Loop once
                       for ibout = 1 : nbouts
                            % Get data (assuming unit dim is singleton or irrelevant)
                            curr_data = squeeze(dataCell{igrp}(imouse, iday, 1, istate, ibout)); % Access first 'unit' index
                            idx_valid = ~isnan(curr_data);
                            nvalid = sum(idx_valid);
                            
                            if nvalid > 0
                                idx_range = idx_curr : (idx_curr + nvalid - 1);
                                vec_data(idx_range) = curr_data(idx_valid);

                                % get bout length (if requested)
                                if contains(frml, ' BoutLength')
                                    curr_blen = blenCell{igrp}(imouse, iday, 1, istate, ibout);
                                    vec_blen(idx_range) = curr_blen;
                                end

                                % Append labels (UnitType will be 'NA' or similar)
                                lbl_grp(idx_range) = repmat(str_grp{igrp}, nvalid, 1);
                                lbl_state(idx_range) = repmat(str_state{istate}, nvalid, 1);
                                lbl_mouse(idx_range) = repmat(mnames{igrp}{imouse}, nvalid, 1);
                                lbl_unit(idx_range) = repmat("NA", nvalid, 1); % No specific unit type
                                lbl_day(idx_range) = repmat(str_day{iday}, nvalid, 1);

                                % Update counter
                                idx_curr = idx_curr + nvalid;
                            end
                       end % end ibout
                    end % end iunit_dummy
                end % end if flg_hasUnits
            end % end iday loop
        end % end istate loop
    end % end imouse loop
end % end igrp loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create output table and config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create table
if contains(frml, ' BoutLength')    % here bout length is a random variable
    lme_tbl = table(vec_data, lbl_grp, lbl_state, lbl_mouse, lbl_unit, vec_blen, id_unit, lbl_day,...
        'VariableNames', {yName, 'Group', 'State', 'Mouse', 'UnitType', 'BoutLength', 'UnitID', 'Day'});
else
    lme_tbl = table(vec_data, lbl_grp, lbl_state, lbl_mouse, lbl_unit, id_unit, lbl_day,...
        'VariableNames', {yName, 'Group', 'State', 'Mouse', 'UnitType', 'UnitID', 'Day'});
end

% Convert to categorical
lme_tbl.State = categorical(lme_tbl.State);
lme_tbl.Mouse = categorical(lme_tbl.Mouse);
lme_tbl.Day = categorical(lme_tbl.Day);
lme_tbl.UnitType = categorical(lme_tbl.UnitType);
lme_tbl.Group = categorical(lme_tbl.Group);

% Makes WT the reference
if any(unique(lme_tbl.Group) == categorical({'WT'}))
    cats = categories(lme_tbl.Group);
    other_cat = cats(~strcmp(cats, 'WT'));
    lme_tbl.Group = reordercats(lme_tbl.Group, ['WT', other_cat]);
end

% Makes BSL the reference
if any(unique(lme_tbl.Day) == categorical({'BSL'}))
    cats = categories(lme_tbl.Day);
    other_cat = cats(~strcmp(cats, 'BSL'));
    lme_tbl.Day = reordercats(lme_tbl.Day, ['BSL'; other_cat]);
end

% Makes pPYR the reference
if any(unique(lme_tbl.UnitType) == categorical({'pPYR'}))
    cats = categories(lme_tbl.UnitType);
    other_cat = cats(~strcmp(cats, 'pPYR'));
    lme_tbl.UnitType = reordercats(lme_tbl.UnitType, ['pPYR'; other_cat]);
end

% Store metadata
lme_cfg.grps = grppaths;
lme_cfg.vars = vars;
lme_cfg.frml = frml;
lme_cfg.flg_emg = flg_emg;

end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fr_data = org_fr(v, frml)

% extracts fr data from concatenated struct and returns a 5D array even
% when some dimensions are singleton [mouse x day x unit x state x bout]

fr = catfields([v(:).fr], 'addim', true);
nmice = length(v);

% fr irrespective of states [unit x 1 x mouse]
if ~contains(frml, 'State')
    fr_data = permute(fr.mfr, [3, 4, 1, 2, 5]);
    return
end

% fr in states averaged across bouts [unit x state x mouse]
if ~contains(frml, 'BoutLength')
    fr_data = permute(fr.states.mfr, [3, 4, 1, 2, 5]);
    return
end

% fr in states per bout {1 x state x mouse}[unit x bout]
fr_tmp = cell(nmice, 1);
for imouse = 1 : nmice
    % [unit x state x bout]
    fr_tmp{imouse} = cell2padmat(fr.states.fr(1, :, imouse), 3);
end
fr_data = permute(cell2padmat(fr_tmp, 4), [4, 5, 1, 3, 2]);

end


% burstiness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function brst_data = org_brst(v, brstField)

% extracts burst data from concatenated struct and returns a 3D array of
% [mouse x day x unit]

fldnames = fieldnames(v(1));
flg_mea = contains(fldnames{1}, 'brst');
if flg_mea
    var_brst = 'brst';
else
    var_brst = 'st';
end

brst = catfields([v(:).(var_brst)], 'addim', true);
brst_data = permute(brst.(brstField), [3, 4, 2, 1]);

end


% FOOOF parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fooof_data = org_1of(v, fooof_param, iband)

% extracts FOOOF parameters from psd_1of struct and returns a 5D array of
% [mouse x day x 1 x state x bout]. Organizes data according to specified
% parameter (cf, pow, bw, etc.) and specified frequency band

% extract and organize data
psd_1of = catfields([v(:).psd_1of], 'addim', true);
nmice = length(v);

% organize data as {1 x state x mouse}[1 x bout]
fooof_tmp = cell(nmice, 1);
for imouse = 1 : nmice

    % handle special cases of aperiodic params
    if contains(fooof_param, 'ap')
        curr_param = psd_1of.(fooof_param)(:, :, imouse)';
    else
        % extract parameter for specific band [state x bout]
        curr_param = psd_1of.(fooof_param)(:, :, iband, imouse);
    end

    fooof_tmp{imouse} = permute(curr_param, [3, 1, 2]);     % add singleton dim [1 x state x bout]
end

% organize as [mouse x day x 1 x state x bout]
fooof_data = permute(cell2padmat(fooof_tmp, 4), [4, 5, 1, 2, 3]);

end


% band power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function band_data = org_band(v, band_idx)

% extracts band power from psd struct and organizes as [mouse x day x 1 x state x bout]

% index to band for analysis. typically:
% bandNames = ["broad", "swa", "delta", "theta", "beta", "gamma"];
flg_norm = true;

% get band power of interest
psd = catfields([v(:).psd], 'addim', true);
nmice = length(v);

% band in states per bout {1 x state x mouse}[unit x bout]
band_tmp = cell(nmice, 1);
for imouse = 1 : nmice
    % [1 x state x bout]
    band_mat = cell2padmat(psd.bands.bouts(1, :, imouse), 3);
    if flg_norm
        band_tmp{imouse} = band_mat(band_idx, :, :) ./ band_mat(1, :, :);
    else
        band_tmp{imouse} = band_mat(band_idx, :, :);
    end
end
band_data = permute(cell2padmat(band_tmp, 4), [4, 5, 1, 3, 2]);

end

% bout length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function blen_data = org_blen(v, frml)

% extracts bout time edges from fr struct, calculates bout length, and
% organizes as [mouse x day x 1 x state x bout]

blen_data = [];
if ~contains(frml, 'BoutLength')
    return
end

% calculate bout length
fr = catfields([v(:).fr], 'addim', true);
binedges = squeeze(fr.states.binedges);
boutLen = cellfun(@(x) cellfun(@(y) diff(y), x, 'uni', true),...
    binedges, 'uni', false);

% ogranize
nmice = length(v);
blen_tmp = cell(nmice, 1);
for imouse = 1:nmice
    % [state x bout]
    blen_tmp{imouse} = cell2padmat(boutLen(:, imouse), 1);
end
blen_data = permute(cell2padmat(blen_tmp, 3), [3, 4, 5, 1, 2]);

end

% units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unit_data = org_units(v)

% extracts unit type data from loaded struct and returns a 5D array
% [mouse x 1 x unit x 1 x 1] where the value indicates unit type (1=PYR, 2=PV).
% Returns empty if 'units' field is not present.

units = catfields([v(:).units], 'addim', true);
nmice = length([v(:).units]);

unit_data_mouse = cell(nmice, 1);
for imouse = 1 : nmice
    clean_mouse = units.clean(:, :, imouse); % [2 x Nunits_mouse]
    nunits_mouse = size(clean_mouse, 2);
    
    unit_type_vec_mouse = zeros(1, nunits_mouse); % Row vector for types
    unit_type_vec_mouse(clean_mouse(1,:)) = 1;    % Type 1 (PYR)
    unit_type_vec_mouse(clean_mouse(2,:)) = 2;    % Type 2 (PV)
    
    % Reshape to [1 x 1 x nunits_mouse x 1 x 1] for consistent padding later
    unit_data_mouse{imouse} = reshape(unit_type_vec_mouse, [1, 1, nunits_mouse, 1, 1]);
end

% Pad across mice and units to create a consistent matrix
% Pad value 0 indicates no unit or type not assigned.
unit_data = cell2padmat(unit_data_mouse, 1, 0); % Pads along mouse dim (dim 1)

% Final output shape: [nmice x 1 x max_nunits_all x 1 x 1]

end

% mice names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mnames = org_mnames(basepaths)
% extracts mice names from first day's paths (assumes consistent across days)

[nmice, ~] = size(basepaths);
mnames = cell(nmice, 1);
for imouse = 1 : nmice
    [~, mname] = fileparts(basepaths(imouse, 1));
    mnames{imouse} = regexp(mname, '^[^_]*', 'match', 'once');
end

end


% string labels for vars according to specific experiments (e.g., mcu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [str_grp, str_day, str_state] = org_strLbls(grppaths, frml, flg_emg)


% group labels from paths (assumes 1st cell is WT and 2nd MCU-KO)
ngrps = length(grppaths);
if ngrps == 2
    str_grp = {'WT'; 'MCU-KO'};
else
    str_grp = split(num2str(1 : ngrps));
end

% state labels
if flg_emg
    str_state = {'High EMG'; 'Low EMG'};
else
    str_state = {'AW'; 'NREM'; 'REM'};
end


% day labels. Assumes 2nd dim of grppaths is session (day)
ndays = unique(cellfun(@(x) size(x ,2), grppaths, 'uni', true));
if contains (frml, ' Day')
    if ndays == 7
        str_day = {'BSL'; 'BAC On'; 'BAC1'; 'BAC2'; 'BAC3'; 'BAC Off'; 'WASH'};
    elseif ndays == 5
        str_day = {'BSL'; 'BAC1'; 'BAC2'; 'BAC3'; 'WASH'};
    end
else
    if ndays == 1
        str_day = {'BSL'};
    else
        str_day = split(num2str(1 : ndays));
    end
end


end

% ripple Gain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ripp_data = org_rippSpks(v, var_field)
% extracts ripple Gain data from ripp struct and returns a 5D array of
% [mouse x day x unit x state x bout]. Assumes Gain is calculated per unit.

nmice = length(v);
gain_tmp = cell(nmice, 1);
max_nunits = 0;

% Extract Gain vector for each mouse and find max number of units
for imouse = 1 : nmice
    gain_vec = v(imouse).ripp.spks.su.(var_field)(:); % Ensure column vector
    nunits_mouse = length(gain_vec);
    gain_tmp{imouse} = gain_vec;
    max_nunits = max(max_nunits, nunits_mouse);
end

% Preallocate and fill data matrix, padding with NaNs
% Output format: [mouse x day x unit x state x bout]
ripp_data = nan(nmice, 1, max_nunits, 1, 1);
for imouse = 1 : nmice
    if ~isempty(gain_tmp{imouse})
        nunits_mouse = length(gain_tmp{imouse});
        ripp_data(imouse, 1, 1:nunits_mouse, 1, 1) = gain_tmp{imouse};
    end
end

end

% specified ripple data field (e.g., peakAmp, peakPower)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ripp_data = org_ripp(v, var_field)
% extracts specified ripple data field from ripp struct and returns a 5D array
% [mouse x day x 1 x 1 x ripple], randomly sampling 100 ripples per mouse.

nmice = length(v);
nRipp = 200; % Max ripples to sample per mouse
rippMat = nan(nmice, nRipp);

% Randomly sample without replacement, specifically from ripples during
% nrem
for imouse = 1 : nmice
        rippTmp = v(imouse).ripp.(var_field);
        idxState = find(v(imouse).ripp.states.idx(:, 4));
        idxTmp = randperm(length(idxState), nRipp);
        idxRipp = idxState(idxTmp);
        rippMat(imouse, :) = rippTmp(idxRipp);
end

% Permute to [mouse x 1 x 1 x 1 x ripple]
ripp_data = permute(rippMat, [1, 3, 4, 5, 2]); % Add singleton dimensions for day, unit, state

end