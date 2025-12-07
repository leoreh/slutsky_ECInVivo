function [fr_rs_mat, fr_fs_mat, t_axis, id_rs, id_fs] = mcu_catFr()
% MCU_CATFR - Concatenate and align firing rates from all mice
%
% Returns:
%   fr_rs_mat: Matrix of RS units (Rows x Time). Aligned to perturbation.
%   fr_fs_mat: Matrix of FS units (Rows x Time). Aligned to perturbation.
%   t_axis:    Time vector relative to perturbation (0 = onset).
%   id_rs:     Vector tracking Mouse ID for each row in fr_rs_mat.
%   id_fs:     Vector tracking Mouse ID for each row in fr_fs_mat.

grp = 'wt';
mice = mcu_sessions(grp);
vars = {'fr'; 'units'};

% Initialize collections
% We will store traces in a struct array first, then pad/matrixify
raw_rs = struct('trace', {}, 't_start', {}, 'mouse', {});
raw_fs = struct('trace', {}, 't_start', {}, 'mouse', {});

% Global time bounds (relative to perturbation)
t_global_min = inf;
t_global_max = -inf;

flg_plot = true;

for iMouse = 1 : length(mice)

    % 1. Load All Data for Mouse
    basepaths = mcu_sessions(mice{iMouse});
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);

    % 2. Find Perturbation Onset for the Mouse
    % We reconstruct the "stitched" timeline for the mouse to find the global
    % perturbation event.
    % Extract all STRD traces (Units x Time) for each day
    all_strd = {};
    for i = 1:length(v)
        % Ensure consistent orientation (Units x Time)
        all_strd{i} = v(i).fr.strd;
    end

    % Use cell2padmat to create a "chimera" matrix just for finding the mean
    % frt_stitched: (MaxUnits x TotalTime)
    frt_stitched = cell2padmat(all_strd, 2);

    % Calculate Mean FR of the mouse (averaging all units)
    mfr_mouse = mean(frt_stitched, 1, 'omitnan');

    % Denoise / Smooth if necessary (following original code logic)
    % t_dummy is just for the function call, absolute time doesn't matter yet
    t_dummy = 1:length(mfr_mouse);
    mfr_smooth = mea_frDenoise(mfr_mouse, t_dummy, 'flgPlot', false, 'frameLenSec', 300);

    % Find Change Point
    pertWin = 1 : min(3000, length(mfr_smooth)); % Use original window or bounded
    [pertOnset, ~] = findchangepts(mfr_smooth(pertWin), ...
        'Statistic', 'mean', 'minDistance', 5, ...
        'MaxNumChanges', 1);

    if flg_plot
        figure; plot(mfr_smooth); hold on; xline(pertOnset, 'r');
        title(['Mouse ' mice{iMouse} ' Perturbation: ' num2str(pertOnset)]);
    end

    % 3. Extract and Align Individual Units
    % Iterate through days again to pull specific units
    current_time_idx = 0;

    for iFile = 1 : length(v)

        % Data for this day
        traces = v(iFile).fr.strd; % Units x Time
        unit_types = v(iFile).units.clean; % 0=RS, 1=FS

        [nUnits, nSamples] = size(traces);

        % Time range for this file relative to perturbation
        % Start: (Cumulative + 1) - PertOnset
        % End:   (Cumulative + nSamples) - PertOnset
        rel_start = (current_time_idx + 1) - pertOnset;
        rel_end   = (current_time_idx + nSamples) - pertOnset;

        % Update global bounds
        t_global_min = min(t_global_min, rel_start);
        t_global_max = max(t_global_max, rel_end);

        % Separate Units
        % Assume: 0 = pPYR (RS), 1 = pINT (FS)

        % RS Units
        idx_rs = unit_types(1, :) == 1;
        if any(idx_rs)
            rs_traces = traces(idx_rs, :);
            n_rs = size(rs_traces, 1);
            for u = 1:n_rs
                new_entry.trace = rs_traces(u, :);
                new_entry.t_start = rel_start;
                new_entry.mouse = iMouse;
                raw_rs(end+1) = new_entry;
            end
        end

        % FS Units
        idx_fs = unit_types(2, :) == 1;
        if any(idx_fs)
            fs_traces = traces(idx_fs, :);
            n_fs = size(fs_traces, 1);
            for u = 1:n_fs
                new_entry.trace = fs_traces(u, :);
                new_entry.t_start = rel_start;
                new_entry.mouse = iMouse;
                raw_fs(end+1) = new_entry;
            end
        end

        % Advance cumulative time
        current_time_idx = current_time_idx + nSamples;
    end
end

% 4. Construct Final Matrices
% Create time axis
t_axis = t_global_min : t_global_max;
n_total_points = length(t_axis);

% Helper to build matrix
    function [mat, ids] = build_matrix(raw_struct)
        n_rows = length(raw_struct);
        mat = nan(n_rows, n_total_points);
        ids = splitapply(@(x) x, [raw_struct.mouse], 1:n_rows)'; % Or simple extraction
        ids = [raw_struct.mouse]';

        for k = 1:n_rows
            % Find indices in the global matrix
            % signal starts at t_start.
            % matrix starts at t_global_min.
            % offset = t_start - t_global_min + 1
            offset = raw_struct(k).t_start - t_global_min + 1;
            len = length(raw_struct(k).trace);

            % Insert
            mat(k, offset : offset + len - 1) = raw_struct(k).trace;
        end
    end

[fr_rs_mat, id_rs] = build_matrix(raw_rs);
[fr_fs_mat, id_fs] = build_matrix(raw_fs);

% Verification Plot
if flg_plot
    figure;
    subplot(2,1,1);
    plot(t_axis, mean(fr_rs_mat, 1, 'omitnan'), 'k');
    xline(0, 'r--'); title('Grand Average RS (n=' + string(length(id_rs)) + ')');

    subplot(2,1,2);
    plot(t_axis, mean(fr_fs_mat, 1, 'omitnan'), 'r');
    xline(0, 'r--'); title('Grand Average FS (n=' + string(length(id_fs)) + ')');
end

end
