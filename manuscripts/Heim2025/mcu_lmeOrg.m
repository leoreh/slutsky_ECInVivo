function [lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, flg_emg)

% organizes data from multiple sessions into format compatible with lme
% analysis. loads data according to mcu_sessions and arranges in a table.
% data can be firing rate or burstiness, determined by the formula. 
%
% INPUT
%   grppaths    cell (per group) of string arrays [mouse x session]
%   frml        char of fomrula for lme. determines how the table is
%               organized
%   flg_emg     logical. load fr by as states (false) or emg states {true}
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    flg_emg = false;
end

% set labels according to specific experiment (e.g., mcu) and formula
ngrps = length(grppaths);
[str_grp, str_day, str_state] = org_strLbls(grppaths, frml, flg_emg);
str_unit = {'pPYR', 'pPV'};

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
    vars = {'st_brst'; 'units'};
    % vars = {'st_metrics'; 'units'};
    yName = 'Burst';

elseif contains(frml, 'BLen ~') || contains(frml, ' BoutLength')
    % bout lengths extracted from fr struct
    vars = {'fr'};       
    yName = 'BLen';
end

if contains(frml, 'Band ~')
    if flg_emg
        vars = [vars, {'psdEmg'}];
    else
        vars = [vars, {'psd'}];
    end
    yName = 'Band';
    
elseif contains(frml, 'FOOOF ~')
   if flg_emg
        vars = [vars, {'psdEmg_1of'}];
    else
        vars = [vars, {'psd_1of'}];
   end
   yName = 'FOOOF';
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
        basepaths = grppaths{igrp}(:, iday);
        v = basepaths2vars('basepaths', basepaths, 'vars', vars);

        % Extract data
        u_day{iday} = org_units(v);                                         % units
        blen_day{iday} = org_blen(v, frml);                                 % bout length (fixed / random effect)

        if contains(frml, 'FR ~')
            data_day{iday} = org_fr(v, frml);                               % firing rate

        elseif contains(frml, 'Burst ~')
            data_day{iday} = org_brst(v);                                   % burstiness

        elseif contains(frml, 'Band ~')
            data_day{iday} = org_band(v);                                   % band power

        elseif contains(frml, 'BLen ~')
            data_day{iday} = org_blen(v, 'BoutLength');                     % bout length (response variable)

        elseif contains(frml, 'FOOOF ~')
            data_day{iday} = org_1of(v);                                    % psd foood parameter
        end
    end

    % Concatenate across days
    dataCell{igrp} = cell2padmat(data_day, 2);                              % [mouse x day x unit x state x bout]
    uCell{igrp} = cell2padmat(u_day, 2);                                    % [mouse x day x unit x type]
    blenCell{igrp} = cell2padmat(blen_day, 2);                              % [mouse x day x 1 x state x bout]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize for lme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize empty arrays
fr_vec = [];
blen_vec = [];
lbl_grp = [];
lbl_state = [];
lbl_mouse = [];
lbl_day = [];
lbl_unit = [];
id_unit = [];

% Organize data into table format
for igrp = 1 : ngrps
    [nmice, ndays, ~, nstates, nbouts] = size(dataCell{igrp});

    for imouse = 1 : nmice
        for iday = 1 : ndays
            for istate = 1 : nstates
                for iunit = 1 : 2

                    % Get indices for current unit type
                    if ~isempty(uCell{igrp})
                        unit_idx = squeeze(uCell{igrp}(imouse, iday, :, iunit));
                        curr_unit_ids = find(unit_idx)' + ...
                            1000 * (imouse - 1) + ...
                            10000 * (iday - 1) + ...
                            100000 * (igrp - 1);
                    else
                        unit_idx = 1;
                        curr_unit_ids = 1;
                    end

                    for ibout = 1 : nbouts

                        % Get firing rate and bout length
                        curr_data = squeeze(dataCell{igrp}(imouse, iday, unit_idx, istate, ibout));
                        valid_idx = ~isnan(curr_data);

                        if contains(frml, 'BoutLength')
                            curr_blen = blenCell{igrp}(imouse, iday, 1, istate, ibout);
                            valid_idx = ~isnan(curr_data) & valid_idx;
                            blen_vec = [blen_vec; repmat(curr_blen, sum(valid_idx), 1)];
                        end

                        % Append to vectors
                        fr_vec = [fr_vec; curr_data(valid_idx)];
                        lbl_grp = [lbl_grp; repmat(str_grp(igrp), sum(valid_idx), 1)];
                        lbl_state = [lbl_state; repmat(str_state(istate), sum(valid_idx), 1)];
                        lbl_mouse = [lbl_mouse; repmat(mnames{igrp}{imouse}, sum(valid_idx), 1)];
                        lbl_unit = [lbl_unit; repmat(str_unit(iunit), sum(valid_idx), 1)];
                        id_unit = [id_unit; curr_unit_ids(valid_idx)'];
                        lbl_day = [lbl_day; repmat(str_day(iday), sum(valid_idx), 1)];
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create output table and config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create table
if contains(frml, ' BoutLength')    % here bout length is a random variable
    lme_tbl = table(fr_vec, lbl_grp, lbl_state, lbl_mouse, lbl_unit, blen_vec, id_unit, lbl_day,...
        'VariableNames', {yName, 'Group', 'State', 'Mouse', 'UnitType', 'BoutLength', 'UnitID', 'Day'});
else
    lme_tbl = table(fr_vec, lbl_grp, lbl_state, lbl_mouse, lbl_unit, lbl_day,...
        'VariableNames', {yName, 'Group', 'State', 'Mouse', 'UnitType', 'Day'});
end

% Convert to categorical
lme_tbl.Group = categorical(lme_tbl.Group);
lme_tbl.State = categorical(lme_tbl.State);
lme_tbl.Mouse = categorical(lme_tbl.Mouse);
lme_tbl.Day = categorical(lme_tbl.Day);

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
function brst_data = org_brst(v)

% extracts burst data from concatenated struct and returns a 3D array of
% [mouse x day x unit]

% brstVar = 'lidor';
% st = catfields([v(:).st], 'addim', true);
% brst_data = permute(st.(brstVar), [3, 4, 2, 1]);

brstVar = 'rateNorm';
brst = catfields([v(:).brst], 'addim', true);
brst_data = permute(brst.(brstVar), [3, 4, 2, 1]);

end


% FOOOF parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fooof_data = org_1of(v)

% extracts FOOOF parameters from psd_1of struct and returns a 5D array of
% [mouse x day x 1 x state x bout]. Organizes data according to specified
% parameter (cf, pow, bw, etc.) and specified frequency band

% parameter and frequency band to analyze
fooof_param = 'ap_exp';        
iband = 1;                

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
function band_data = org_band(v)

% extracts band power from psd struct and organizes as [mouse x day x 1 x state x bout]

% index to band for analysis. typically: 
% bandNames = ["broad", "swa", "delta", "theta", "beta", "gamma"];
% see calc_bands.m for updated bands
band_idx = 6;
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

% extracts unit data from loaded struct and returns a 3D [mouse x day x unit x type]
unit_data = [];
if ~isfield(v, 'units')
    return
end
units = catfields([v(:).units], 'addim', true);
unit_data = permute(units.clean, [3, 4, 2, 1]);

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
if contains (frml, ' Group')
    if ngrps == 2
        str_grp = {'WT'; 'MCU-KO'};
    end
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
    end
else
    str_day = split(num2str(1 : ndays));
end


end