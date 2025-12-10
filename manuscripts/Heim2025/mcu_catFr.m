function [fr_rs, fr_fs, t_axis, id_rs, id_fs] = mcu_catFr(fetTbl)
% MCU_CATFR - Concatenate and align firing rates from all mice
%
% Returns:
%   fr_rs: Matrix of RS units (Rows x Time). Aligned to perturbation.
%          Rows are packed by Mouse (MaxUnits per mouse).
%   fr_fs: Matrix of FS units (Rows x Time). Aligned to perturbation.
%   t_axis: Time vector relative to perturbation (in Hours).
%   id_rs: Mouse ID for each row in fr_rs (e.g. 96 for 'lh96').
%   id_fs: Mouse ID for each row in fr_fs.

if nargin < 1, fetTbl = []; end

grp = 'wt';
mice = unique(get_mname(mcu_basepaths(grp)));
vars = {'fr'; 'units'};

% -------------------------------------------------------------------------
% Pass 1: Stitch, Find Perturbation, Analyze Dimensions
% -------------------------------------------------------------------------

data_store = struct();
global_min = inf;
global_max = -inf;

flgPlot = true;

for iMouse = 1 : length(mice)

    % Get Mouse ID Number (e.g. 'lh96' -> 96)
    mStr = mice{iMouse};
    mNum = str2double(regexp(mStr, '\d+', 'match', 'once'));
    if isnan(mNum), mNum = iMouse; end % Fallback

    % 1. Load All Data for Mouse
    basepaths = mcu_basepaths(mice{iMouse});
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);

    % 2. Stitch and Find Perturbation
    % Extract all STRD traces for stitching (Unit agnostic for now)
    all_strd = {};
    for i = 1:length(v)
        all_strd{i} = v(i).fr.strd;
    end
    frt_stitched = cell2padmat(all_strd, 2);

    % Mean FR for perturbation detection
    mfr_mouse = mean(frt_stitched, 1, 'omitnan');
    t_dummy = 1:length(mfr_mouse);
    mfr_smooth = mea_frDenoise(mfr_mouse, t_dummy, 'flgPlot', false, 'frameLenSec', 300);

    % Find Perturbation
    pertWin = 1 : min(3000, length(mfr_smooth));
    [pertOnset, ~] = findchangepts(mfr_smooth(pertWin), ...
        'Statistic', 'mean', 'minDistance', 5, ...
        'MaxNumChanges', 1);

    if flgPlot
        figure; plot(mfr_smooth); hold on; xline(pertOnset, 'r');
        title(['Mouse ' mice{iMouse}]);
    end

    % 3. Analyze Unit Counts and Time
    curr_time = 0;
    max_rs = 0;
    max_fs = 0;
    day_meta = struct();

    for iDay = 1 : length(v)
        % Identify Unit Types for this day
        if ~isempty(fetTbl)
            [~, fName] = fileparts(basepaths{iDay});

            % Find rows for this file
            isMatch = string(fetTbl.File) == string(fName);
            subTbl = fetTbl(isMatch, :);

            % Extract types. Assume order matches v(iDay)
            idx_rs = (subTbl.unitType == 1)';
            idx_fs = (subTbl.unitType == 2)';

        else
            u_clean = v(iDay).units.clean;
            idx_rs = u_clean(1, :) == 1; % RS: Row 1 == 1
            idx_fs = u_clean(2, :) == 1; % FS: Row 2 == 1
        end

        max_rs = max(max_rs, sum(idx_rs));
        max_fs = max(max_fs, sum(idx_fs));

        % Time Alignment
        traces = v(iDay).fr.strd;
        nSamples = size(traces, 2);

        rel_start = curr_time + 1 - pertOnset;
        rel_end   = curr_time + nSamples - pertOnset;

        global_min = min(global_min, rel_start);
        global_max = max(global_max, rel_end);

        % Store info for Pass 2
        day_meta(iDay).fr = traces;
        day_meta(iDay).idx_rs = idx_rs;
        day_meta(iDay).idx_fs = idx_fs;
        day_meta(iDay).t_start = rel_start;

        curr_time = curr_time + nSamples;
    end

    % Store Mouse Data
    data_store(iMouse).days = day_meta;
    data_store(iMouse).max_rs = max_rs;
    data_store(iMouse).max_fs = max_fs;
    data_store(iMouse).id = mNum;
end

% -------------------------------------------------------------------------
% Pass 2: Construct Final Matrices
% -------------------------------------------------------------------------

t_axis = (global_min : global_max) / 60; % Hours
n_total_points = length(t_axis);

% Pre-allocate based on MAX rows per mouse (summed)
total_rs_rows = sum([data_store.max_rs]);
total_fs_rows = sum([data_store.max_fs]);

fr_rs = nan(total_rs_rows, n_total_points);
fr_fs = nan(total_fs_rows, n_total_points);

id_rs = nan(total_rs_rows, 1);
id_fs = nan(total_fs_rows, 1);

curr_row_rs = 1;
curr_row_fs = 1;

for iMouse = 1 : length(mice)
    d = data_store(iMouse);

    % Calculate Row Ranges for this Mouse
    if d.max_rs > 0
        range_rs = curr_row_rs : (curr_row_rs + d.max_rs - 1);
        id_rs(range_rs) = d.id;
    else
        range_rs = [];
    end

    if d.max_fs > 0
        range_fs = curr_row_fs : (curr_row_fs + d.max_fs - 1);
        id_fs(range_fs) = d.id;
    else
        range_fs = [];
    end

    % Fill Days
    for iDay = 1 : length(d.days)
        day = d.days(iDay);

        % Column Indices
        offset = day.t_start - global_min + 1;
        col_idx = offset : (offset + size(day.fr, 2) - 1);

        % Extract and Fill RS
        if any(day.idx_rs)
            traces = day.fr(day.idx_rs, :);
            n = size(traces, 1);
            if ~isempty(range_rs)
                fr_rs(range_rs(1:n), col_idx) = traces;
            end
        end

        % Extract and Fill FS
        if any(day.idx_fs)
            traces = day.fr(day.idx_fs, :);
            n = size(traces, 1);
            if ~isempty(range_fs)
                fr_fs(range_fs(1:n), col_idx) = traces;
            end
        end
    end

    % Advance Row Counters
    curr_row_rs = curr_row_rs + d.max_rs;
    curr_row_fs = curr_row_fs + d.max_fs;
end

% -------------------------------------------------------------------------
% Visualization
% -------------------------------------------------------------------------
if flgPlot
    figure;
    subplot(2,1,1);
    plot(t_axis, mean(fr_rs, 1, 'omitnan'), 'k');
    xline(0, 'r--');
    title(['Grand Average RS (Rows=' num2str(size(fr_rs,1)) ')']);

    subplot(2,1,2);
    plot(t_axis, mean(fr_fs, 1, 'omitnan'), 'r');
    xline(0, 'r--');
    title(['Grand Average FS (Rows=' num2str(size(fr_fs,1)) ')']);
end

end
