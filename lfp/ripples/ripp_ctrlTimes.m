function ctrlTimes = ripp_ctrlTimes(rippTimes, varargin)
% RIPP_CTRLTIMES Finds non-ripple intervals of matched duration.
%
% SUMMARY:
%   Generates control events (matched duration) in non-ripple periods.
%   Uses a fast, vectorized approach to merge exclusion zones and linear
%   allocation to place control events.
%
% ALGORITHM:
%   1. Define Exclusion Zones: Ripple times +/- padding.
%   2. Merge Overlaps: Consolidate overlapping exclusions into single blocks.
%   3. Find Gaps: Invert merged exclusions to find 'Available Intervals'.
%   4. Linear Allocation: Iterate through ripples and place a matched control
%      event in the first available gap that satisfies constraints.
%
% INPUTS:
%   rippTimes   - (N x 2) Start and end times of ripples [s].
%   recDur      - (Num) Duration of the recording [s].
%   varargin    - 'maxDist' (default 4h), 'padding' (default 0.1s).
%
% OUTPUT:
%   ctrlTimes   - (N x 2) Start and end times of control events [s].

%% ========================================================================
%  INITIALIZATION
%  ========================================================================
p = inputParser;
addRequired(p, 'rippTimes', @isnumeric);
addParameter(p, 'recDur', [], @isnumeric);
addParameter(p, 'maxDist', 4 * 60 * 60, @isnumeric); 
addParameter(p, 'padding', 0.1, @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);
parse(p, rippTimes, varargin{:});

recDur = p.Results.recDur;
maxDist = p.Results.maxDist;
padding = p.Results.padding;
flgPlot = p.Results.flgPlot;

% If recDur empty, assume recording terminates with last ripple
if isempty(recDur)
    recDur = max(rippTimes(:));
end

% Ripple Durations
nRipp = size(rippTimes, 1);
rippDur = rippTimes(:, 2) - rippTimes(:, 1);

%% ========================================================================
%  FIND AVAILABLE GAPS (Vectorized)
%  ========================================================================

% Define Exclusion Zones (Ripples + Padding)
excl = [rippTimes(:,1) - padding, rippTimes(:,2) + padding];

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
% Remove gaps with negative duration (overlaps) or too short for ANY ripple.
gaps(diff(gaps, [], 2) < min(rippDur), :) = [];
nGaps = size(gaps, 1);

if isempty(gaps)
    warning('No valid control intervals found (Recording too saturated).');
    ctrlTimes = nan(nRipp, 2);
    return;
end

%% ========================================================================
%  ALLOCATE CONTROL EVENTS
%  ========================================================================

ctrlTimes = nan(nRipp, 2);

currGap = 1;            % Index of current gap being filled
tCursor = gaps(1, 1);   % Current placement cursor within the gap

for iRipp = 1:nRipp
    
    neededDur = rippDur(iRipp);
    refTime = mean(rippTimes(iRipp, :));
    
    % --- Find Valid Gap ---
    while currGap <= nGaps
        gapEnd = gaps(currGap, 2);
        
        % Check 1: Is this gap "too old"? (maxDist constraint)
        % If the gap ended hours ago relative to this ripple, skip it.
        if gapEnd < (refTime - maxDist)
            currGap = currGap + 1;
            if currGap <= nGaps
                tCursor = gaps(currGap, 1);
            end
            continue;
        end
        
        % Check 2: Does the ripple fit in the remaining space?
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

    % Ripples
    nRipp = size(rippTimes, 1);
    xRipp = [rippTimes'; nan(1, nRipp)] / 3600; 
    
    % Controls 
    nCtrl = size(ctrlTimes, 1);
    xCtrl = [ctrlTimes'; nan(1, nCtrl)] / 3600;
    
    % Y value
    yVal = repmat([1; 1; nan], nRipp, 1);

    % Plot
    figure('Color', 'w', 'Name', 'Ripple vs Control Check');
    
    % Plot Ripples in Red
    hRipp = plot(xRipp(:), yVal, 'Color', [0.8 0.2 0.2], 'LineWidth', 15); hold on;
    
    % Plot Controls in Black
    hCtrl = plot(xCtrl(:), yVal - 0.1, 'k', 'LineWidth', 15);
    
    % Formatting
    ylim([0 2]);
    xlabel('Time (hr)');
    title(sprintf('Event Timing Check (N=%d)', nRipp));
    set(gca, 'TickDir', 'out', 'Box', 'off', 'YDir', 'normal');
    
    % Add legend
    legend([hRipp, hCtrl], {'Ripples', 'Control Events'}, 'Location', 'best');
end


end         % EOF