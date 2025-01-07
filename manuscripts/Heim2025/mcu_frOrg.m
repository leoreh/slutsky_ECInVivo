function [fr_tbl, fr_cfg] = mcu_frOrg(grppaths, frml, flg_emg)

% organizes firing rate data from multiple sessions into format compatible with
% lme analysis. loads data according to mcu_sessions and arranges in a table
%
% INPUT
%   grppaths    cell (per group) of string arrays [mouse x session]
%   frml        char of fomrula for lme. determines how the table is
%               organized
%   flg_emg     logical. load fr by as states (false) or emg states {true}
%
% OUTPUT
%   fr_data    table organized for lme with fields FR, Group, State, Mouse, UnitType
%   fr_cfg     struct with metadata and analysis parameters
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

% set which variable to organize, according to formula
if contains(frml, 'FR ~')
    if flg_emg
        vars = {'frEmg'; 'units'};
    else
        vars = {'fr'; 'units'};
    end
elseif contains(frml, 'Burst ~')
    vars = {'st_brst'; 'units'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize cells for raw data
ngrps = length(grppaths);
frCell = cell(1, ngrps);
uCell = cell(1, ngrps);
blenCell = cell(1, ngrps);
brstCell = cell(1, ngrps);
mnames = cell(1,  ngrps);

for igrp = 1 : ngrps
    
    % extract mice names for labeling
    mnames{igrp} = org_mnames(grppaths{igrp});

    % initialize
    [~, ndays] = size(grppaths{igrp});
    fr_day = cell(1, ndays);
    u_day = cell(1, ndays);
    blen_day = cell(1, ndays);
    brst_day = cell(1, ndays);

    % Loop through days and load data
    for iday = 1 : ndays

        % Load vars for each day
        basepaths = grppaths{igrp}(:, iday);
        v = basepaths2vars('basepaths', basepaths, 'vars', vars);

        % Extract data using standardized format functions
        u_day{iday} = org_units(v);
        fr_day{iday} = org_fr(v, frml);
        brst_day{iday} = org_brst(v);
        blen_day{iday} = org_blen(v, frml);
    end

    % Concatenate across days
    frCell{igrp} = cell2padmat(fr_day, 2);                                  % [mouse x day x unit x state x bout]
    uCell{igrp} = cell2padmat(u_day, 2);                                    % [mouse x day x unit x type]
    blenCell{igrp} = cell2padmat(blen_day, 2);                              % [mouse x day x state x bout]
    brstCell{igrp} = cell2padmat(blen_day, 2);                              % [mouse x day x unit]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize for lme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize empty arrays
fr_vec = [];
grp_label = [];
state_label = [];
mouse_label = [];
unit_label = [];
blen_vec = [];
unit_id = [];
day_label = [];

% Organize data into table format
for igrp = 1 : ngrps
    grp_idx = igrp - 1;
    [nmice, ndays, ~, nstates, nbouts] = size(frCell{igrp});

    for imouse = 1 : nmice
        for iday = 1 : ndays
            for istate = 1 : nstates
                for iunit = 1 : 2

                    % Get indices for current unit type
                    unit_idx = squeeze(uCell{igrp}(imouse, iday, :, iunit));
                    curr_unit_ids = find(unit_idx)' + ...
                        1000 * (imouse - 1) + ...
                        10000 * (iday - 1) + ...
                        100000 * (igrp - 1);

                    for ibout = 1 : nbouts

                        % Get firing rate and bout length
                        curr_fr = squeeze(frCell{igrp}(imouse, iday, unit_idx, istate, ibout));
                        valid_idx = ~isnan(curr_fr);

                        if contains(frml, 'BoutLength')
                            curr_blen = blenCell{igrp}(imouse, iday, istate, ibout);
                            valid_idx = ~isnan(curr_fr) & valid_idx;
                            blen_vec = [blen_vec; repmat(curr_blen, sum(valid_idx), 1)];
                        end

                        % Append to vectors
                        fr_vec = [fr_vec; curr_fr(valid_idx)];
                        grp_label = [grp_label; repmat(grp_idx, sum(valid_idx), 1)];
                        state_label = [state_label; repmat(istate, sum(valid_idx), 1)];
                        mouse_label = [mouse_label; repmat(mnames{igrp}{imouse}, sum(valid_idx), 1)];
                        unit_label = [unit_label; repmat(iunit, sum(valid_idx), 1)];
                        unit_id = [unit_id; curr_unit_ids(valid_idx)'];
                        day_label = [day_label; repmat(iday, sum(valid_idx), 1)];
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
if contains(frml, 'BoutLength')
    fr_tbl = table(fr_vec, grp_label, state_label, mouse_label, unit_label, blen_vec, unit_id, day_label,...
        'VariableNames', {'FR', 'Group', 'State', 'Mouse', 'UnitType', 'BoutLength', 'UnitID', 'Day'});
else
    fr_tbl = table(fr_vec, grp_label, state_label, mouse_label, unit_label, day_label,...
        'VariableNames', {'FR', 'Group', 'State', 'Mouse', 'UnitType', 'Day'});
end

% Convert to categorical
fr_tbl.Group = categorical(fr_tbl.Group);
fr_tbl.State = categorical(fr_tbl.State);
fr_tbl.Mouse = categorical(fr_tbl.Mouse);
fr_tbl.Day = categorical(fr_tbl.Day);

% Store metadata
fr_cfg.grps = grppaths;
fr_cfg.vars = vars;
fr_cfg.frml = frml;
fr_cfg.flg_emg = flg_emg;

end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fr_data = org_fr(v, frml)
% extracts fr data from concatenated struct and returns a 5D array even
% when some dimensions are singleton [mouse x day x unit x state x bout]

fr_data = [];
if ~isfield(v, 'fr')
    return
end

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
for imouse = 1:nmice
    % [unit x state x bout]
    fr_tmp{imouse} = cell2padmat(fr.states.fr(1, :, imouse), 3);
end
fr_data = permute(cell2padmat(fr_tmp, 4), [4, 5, 1, 3, 2]);

end


function brst_data = org_brst(v)
% extracts burst data from concatenated struct and returns a 3D array of
% [mouse x day x unit]

brst_data = [];
if ~isfield(v, 'brst')
    return
end

brst = catfields([v(:).brst], 'addim', true);
brst_data = permute(brst.spkprct, [3, 4, 2, 1]);

end


function blen_data = org_blen(v, frml)
% extracts bout time edges from fr struct, calculates bout length, and
% organizes as [mouse x day x state x bout]

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
blen_data = permute(cell2padmat(blen_tmp, 3), [3, 4, 1, 2]);

end


function unit_data = org_units(v)
% extracts unit data from loaded struct and returns a 3D [mouse x day x unit x type]

units = catfields([v(:).units], 'addim', true);
unit_data = permute(units.clean, [3, 4, 2, 1]);

end


function mnames = org_mnames(basepaths)
% extracts mice names from first day's paths (assumes consistent across days)
        
[nmice, ndays] = size(basepaths);
for imouse = 1:nmice
    [~, mname] = fileparts(basepaths(imouse, 1));
    mnames{imouse} = regexp(mname, '^[^_]*', 'match', 'once');
end

end

