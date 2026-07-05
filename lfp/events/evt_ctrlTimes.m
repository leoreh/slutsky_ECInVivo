function ctrlTimes = evt_ctrlTimes(evtTimes, varargin)
% EVT_CTRLTIMES Generates matched-duration control intervals.
%
%   ctrlTimes = EVT_CTRLTIMES(evtTimes, varargin)
%
%   SUMMARY:
%       Identifies non-event intervals of identical duration to the detected events.
%       Ensures control events do not overlap with any event (plus padding) and
%       are allocated linearly to fill available "safe" gaps in the recording.
%
%   ALGORITHM:
%       1. Define Exclusion Zones: Union of all Event times +/- padding.
%       2. Merge Overlaps: Consolidate overlapping exclusions into continuous blocks.
%       3. Find Gaps: Invert merged exclusions to identify 'Available Intervals'.
%       4. Linear Allocation: Iterate through events and place a matched control
%          event in the first available gap that fits and is within 'maxDist'.
%
%   INPUTS:
%       evtTimes   - (Mat) [N x 2] Start and end times of events [s].
%       varargin    - Parameter/Value pairs:
%           'vldTimes' - (Mat) [M x 2] Intervals allowed for control events.
%                             If empty, the entire recording is allowed.
%           'recDur'   - (Num) Total recording duration [s]. (Default: max(evtTimes)).
%           'maxDist'  - (Num) Max temporal distance allowed for a control event
%                             relative to its event [s]. (Default: 4 hours).
%           'padding'  - (Num) Separation from event boundaries [s]. (Default: 0.1).
%           'flgPlot'  - (Log) Validate with a visual plot? (Default: false).
%
%   OUTPUTS:
%       ctrlTimes   - (Mat) [N x 2] Start and end times of control events [s].
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Updated: 23 Jan 2026

%% ========================================================================
%  INITIALIZATION
%  ========================================================================
p = inputParser;
addRequired(p, 'evtTimes', @isnumeric);
addParameter(p, 'vldTimes', [], @isnumeric);
addParameter(p, 'recDur', [], @isnumeric);
addParameter(p, 'maxDist', 1 * 3600, @isnumeric);
addParameter(p, 'padding', 0.1, @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'lbl', 'Event', @ischar);
parse(p, evtTimes, varargin{:});

vldTimes = p.Results.vldTimes;
recDur = p.Results.recDur;
maxDist = p.Results.maxDist;
padding = p.Results.padding;
flgPlot = p.Results.flgPlot;
lbl = p.Results.lbl;

% If recDur empty, assume recording terminates with last event
if isempty(recDur)
    if ~isempty(vldTimes)
        recDur = max([evtTimes(:); vldTimes(:)]);
    else
        recDur = max(evtTimes(:));
    end
end

% Event Durations
nEvt = size(evtTimes, 1);
evtDur = evtTimes(:, 2) - evtTimes(:, 1);

%% ========================================================================
%  FIND AVAILABLE GAPS 
%  ========================================================================

% Define Exclusion Zones (Events + Padding)
excl = [evtTimes(:,1) - padding, evtTimes(:,2) + padding];

% Handle Inclusion Intervals
if ~isempty(vldTimes)
    
    % Sort and Merge Inclusions
    incOrd = sortrows(vldTimes);
    if size(incOrd, 1) > 1
        % Determine where a new block starts
        % A new block starts if the current start is > previous cumulative max end
        cumMaxEndsInc = cummax(incOrd(:, 2));
        isNewBlockInc = [true; incOrd(2:end, 1) > cumMaxEndsInc(1:end-1)];
        mergedInc = [incOrd(isNewBlockInc, 1), [cumMaxEndsInc(isNewBlockInc(2:end)); cumMaxEndsInc(end)]];
    else
        mergedInc = incOrd;
    end

    % Invert Inclusions relative to RecDur to get "Invalid Zones"
    invExcl = [];
    if mergedInc(1, 1) > 0
        invExcl = [invExcl; 0, mergedInc(1, 1)];
    end
    
    % Gaps between inclusions
    if size(mergedInc, 1) > 1
        invExcl = [invExcl; mergedInc(1:end-1, 2), mergedInc(2:end, 1)];
    end
    
    % Gap from last inclusion to end
    if mergedInc(end, 2) < recDur
        invExcl = [invExcl; mergedInc(end, 2), recDur];
    end

    % Add to Exclusions
    excl = [excl; invExcl];
end

% Crop to Recording Limits
excl(excl(:, 1) < 0, 1) = 0;
excl(excl(:, 2) > recDur, 1) = recDur;

% Fast Merge of Overlapping Exclusions
% Identifies continuous blocks of "unavailable" time.
cumMaxEnds = cummax(excl(:, 2));
isNewBlock = [true; excl(2:end, 1) > cumMaxEnds(1:end-1)];

mergedStarts = excl(isNewBlock, 1);
mergedEnds   = [cumMaxEnds(isNewBlock(2:end)); cumMaxEnds(end)];
mergedExcl   = [mergedStarts, mergedEnds];

% Invert to Find Gaps (Available Time)
% Gaps are the spaces between the end of one block and start of the next.
gapStarts = [0; mergedExcl(:,2)];
gapEnds   = [mergedExcl(:,1); recDur];
gaps      = [gapStarts, gapEnds];

% Filter Invalid Gaps
% Remove gaps with negative duration (overlaps) or too short for ANY event.
gaps(diff(gaps, [], 2) < min(evtDur), :) = [];
nGaps = size(gaps, 1);

if isempty(gaps)
    warning('No valid control intervals found (Recording too saturated).');
    ctrlTimes = nan(nEvt, 2);
    return;
end

%% ========================================================================
%  ALLOCATE CONTROL EVENTS
%  ========================================================================

ctrlTimes = nan(nEvt, 2);

currGap = 1;            % Index of current gap being filled
tCursor = gaps(1, 1);   % Current placement cursor within the gap

for iRipp = 1:nEvt

    neededDur = evtDur(iRipp);
    refTime = mean(evtTimes(iRipp, :));

    % --- Find Valid Gap ---
    while currGap <= nGaps
        gapEnd = gaps(currGap, 2);

        % Check 1: Is this gap "too old"? (maxDist constraint)
        % If the gap ended hours ago relative to this event, skip it.
        if gapEnd < (refTime - maxDist)
            currGap = currGap + 1;
            if currGap <= nGaps
                tCursor = gaps(currGap, 1);
            end
            continue;
        end

        % Check 2: Does the event fit in the remaining space?
        if tCursor + neededDur > gapEnd
            % Doesn't fit. Move to next gap.
            currGap = currGap + 1;
            if currGap <= nGaps
                tCursor = gaps(currGap, 1);
            end
            continue;
        end

        % If we are here, it fits and is valid.
        break;
    end

    % --- Assign Control Event ---
    if currGap <= nGaps
        ctrlTimes(iRipp, :) = [tCursor, tCursor + neededDur];

        % Advance cursor + padding for the NEXT control event
        tCursor = tCursor + neededDur + padding;
    end
end

%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot

    hFig = figure;

    % Events
    nEvt = size(evtTimes, 1);
    xEvt = [evtTimes'; nan(1, nEvt)] / 3600;

    % Controls
    nCtrl = size(ctrlTimes, 1);
    xCtrl = [ctrlTimes'; nan(1, nCtrl)] / 3600;

    % Y value
    yVal = repmat([1; 1; nan], nEvt, 1);

    % Plot
    figure('Color', 'w', 'Name', [lbl ' vs Control Check']);

    % Plot Events in Red
    hEvt = plot(xEvt(:), yVal, 'Color', [0.8 0.2 0.2], 'LineWidth', 15); hold on;

    % Plot Controls in Black
    hCtrl = plot(xCtrl(:), yVal - 0.1, 'k', 'LineWidth', 15);

    % Formatting
    ylim([0 2]);
    xlabel('Time (hr)');
    title(sprintf('Event Timing Check (N=%d)', nEvt));
    set(gca, 'TickDir', 'out', 'Box', 'off', 'YDir', 'normal');

    % Add legend
    legend([hEvt, hCtrl], {lbl, 'Control Events'}, 'Location', 'best');
end


end         % EOF