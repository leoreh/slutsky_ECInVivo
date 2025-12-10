function [tAxis, frTbl] = mcu_frTbl(fetTbl)
% MCU_FRTBL - Create a table of firing rates aligned to perturbation
%
% Usage:
%   [tAxis, frTbl] = mcu_frTbl(fetTbl)
%
% Inputs:
%   fetTbl (optional) - Table containing unit metadata. If provided, the
%                       function assumes that for each file, the rows in
%                       fetTbl correspond 1-to-1 with the unit order in
%                       the loaded firing rate data.
%                       Must contain: Mouse, File.
%                       If empty, defaults to processing all 'wt' mice.
%
% Outputs:
%   tAxis  - Time vector relative to perturbation (in Hours).
%   frTbl  - Table where each row is a unit.
%
% See also: MCU_CATFR, MEAS_FRDENOISE

%% ========================================================================
%  SETUP
%  ========================================================================
if nargin < 1, fetTbl = []; end

if isempty(fetTbl)
    grp = 'wt';
    mice = mcu_sessions(grp);
    useFetTbl = false;
else
    % Use mice present in the table
    mice = cellstr(unique(fetTbl.Mouse));
    useFetTbl = true;
end

vars = {'fr'};
rows = {}; % Collector for table rows
global_min = inf;
global_max = -inf;

%% ========================================================================
%  ITERATE MICE
%  ========================================================================
for iMouse = 1 : length(mice)

    mouseName = mice{iMouse};

    % Load All Data for Mouse
    basepaths = mcu_sessions(mouseName);
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);

    % Stitch and Find Perturbation (Per Mouse)
    % Extract all STRD traces for stitching (Unit agnostic)
    all_strd = {};
    for i = 1:length(v)
        if isfield(v(i).fr, 'strd')
            all_strd{i} = v(i).fr.strd;
        else
            all_strd{i} = [];
        end
    end
    frt_stitched = cell2padmat(all_strd, 2);


    % Mean FR for perturbation detection
    mfr_mouse = mean(frt_stitched, 1, 'omitnan');
    t_dummy = 1:length(mfr_mouse);
    % Denoise mean to find pert
    mfr_smooth = mea_frDenoise(mfr_mouse, t_dummy, 'flgPlot', false, 'frameLenSec', 300);

    % Find Perturbation
    % Limit search to first 3000 pts (arbitrary as per mcu_catFr) or length
    pertWin = 1 : min(3000, length(mfr_smooth));
    [pertOnset, ~] = findchangepts(mfr_smooth(pertWin), ...
        'Statistic', 'mean', 'minDistance', 5, ...
        'MaxNumChanges', 1);

    % Iterate Days/Files
    curr_time = 0;

    for iDay = 1 : length(v)
        [~, fileName] = fileparts(basepaths{iDay});

        % Data for this day
        traces = v(iDay).fr.strd;
        traceLen = size(traces, 2);

        % Determine Rows to process
        % Find rows in fetTbl corresponding to this File
        isMatch = string(fetTbl.File) == string(fileName);
        dayRows = fetTbl(isMatch, :);
        nToProcess = height(dayRows);

        % Process each unit
        for iUnit = 1:nToProcess
            % Assumption: Table row k corresponds to Trace k
            uIdx = iUnit;
            if uIdx > size(traces, 1), continue; end

            % Get Metadata Row
            uMeta = dayRows(iUnit, :);

            % Get FR & Denoise
            raw_fr = traces(uIdx, :);
            dn_fr = mea_frDenoise(raw_fr, 1:length(raw_fr), ...
                'flgPlot', false, 'frameLenSec', 300);

            % Calculate Timing
            rel_start_idx = (curr_time + 1) - pertOnset;
            rel_end_idx   = (curr_time + traceLen) - pertOnset;

            % Update Global Bounds (in indices)
            global_min = min(global_min, rel_start_idx);
            global_max = max(global_max, rel_end_idx);

            % pertOnset (relative to unit's firing rate)
            unit_pert_idx = pertOnset - curr_time;

            % Store in row
            uMeta.pertOnset = unit_pert_idx / 60; % In Hours
            uMeta.FR = {dn_fr};

            % Store placement info for later alignment
            uMeta.tmp_rel_start_idx = rel_start_idx;
            uMeta.tmp_len = traceLen;

            rows{end+1} = uMeta;
        end

        curr_time = curr_time + traceLen;
    end
end

%% ========================================================================
%  CONSTRUCT TABLE
%  ========================================================================
if isempty(rows)
    frTbl = table();
    tAxis = [];
    return;
end

% Vertcat all one-row tables
frTbl = vertcat(rows{:});

%% ========================================================================
%  CONSTRUCT FINAL TIME AXIS & ALIGN
%  ========================================================================
if isinf(global_min)
    tAxis = [];
    frTbl.FR_aligned = [];
    return;
end

tAxis = (global_min : global_max) / 60; % Hours
nPoints = length(tAxis);

% Preallocate aligned FR as a Matrix (Rows x Time)
% nUnits x nPoints
fr_aligned = nan(height(frTbl), nPoints);

for i = 1:height(frTbl)
    % Get data
    fr = frTbl.FR{i};
    rel_start = frTbl.tmp_rel_start_idx(i);

    % Identify bounds in tAxis
    start_idx = rel_start - global_min + 1;
    end_idx = start_idx + length(fr) - 1;

    % Check bounds (sanity check)
    if start_idx < 1 || end_idx > nPoints
        warning('Alignment out of bounds for unit %d', i);
        target_idx = max(1, start_idx) : min(nPoints, end_idx);
        src_idx = (target_idx - start_idx) + 1;
        fr_aligned(i, target_idx) = fr(src_idx);
    else
        fr_aligned(i, start_idx : end_idx) = fr;
    end
end

frTbl.FR_aligned = fr_aligned;

% Cleanup temporary fields
frTbl.tmp_rel_start_idx = [];
frTbl.tmp_len = [];

end
