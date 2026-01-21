function ctrlTimes = ripp_ctrlTimes(rippTimes, recDur, varargin)
% RIPP_CTRLTIMES Finds non-ripple intervals of matched duration.
%
% SUMMARY:
%   Generates control events for ripple analysis by finding non-ripple periods
%   of identical duration to the detected ripples. It attempts to find these
%   periods as close as possible to the original ripple event (within maxDist),
%   avoiding other ripples.
%
% INPUTS:
%   rippTimes   - (N x 2) Start and end times of ripples [s].
%   recDur      - (Num) Duration of the recording [s].
%   varargin    - 'maxDist' (default 4*3600), 'padding' (default 0.1).
%
% OUTPUT:
%   ctrlTimes   - (N x 2) Start and end times of control events [s].
%
% DEPENDENCIES: SubtractIntervals
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'rippTimes', @isnumeric);
addRequired(p, 'recDur', @isnumeric);
addParameter(p, 'maxDist', 4 * 60 * 60, @isnumeric); % 4 Hours
addParameter(p, 'padding', 0.1, @isnumeric);
parse(p, rippTimes, recDur, varargin{:});

maxDist = p.Results.maxDist;
padding = p.Results.padding;

%% ========================================================================
%  FIND CONTROL INTERVALS
%  ========================================================================

nRipples = size(rippTimes, 1);
if nRipples == 0
    ctrlTimes = [];
    return;
end

rippDur = rippTimes(:, 2) - rippTimes(:, 1);
buffer = max(rippDur);

% Define available time (Rec - Ripples +/- buffer)
intervalsExclude = rippTimes + [-buffer, buffer];
intervalsInclude = [0, recDur];
intervalsCtrl = SubtractIntervals(intervalsInclude, intervalsExclude);

% Filter short intervals
% Note: Using ripple max dur as minimum requirement for control blocks
durIdx = diff(intervalsCtrl, 1, 2) < max(rippDur);
intervalsCtrl(durIdx, :) = [];

if isempty(intervalsCtrl)
    warning('No valid control intervals found.');
    ctrlTimes = nan(nRipples, 2);
    return;
end

%% ========================================================================
%  MATCH EVENTS
%  ========================================================================

ctrlTimes = nan(nRipples, 2);
nFree = size(intervalsCtrl, 1);

currentFreeIdx = 1;
lastEndTime = -inf;

for iRipp = 1:nRipples

    needed = rippDur(iRipp);

    % Save current index to avoid re-scanning from 0, but allow looking ahead
    tempIdx = currentFreeIdx;

    while tempIdx <= nFree
        blockStart = intervalsCtrl(tempIdx, 1);
        blockEnd = intervalsCtrl(tempIdx, 2);

        % Candidates must be after last assigned control event + padding
        t_start = max(blockStart, lastEndTime + padding);

        if t_start + needed <= blockEnd
            % Check distance condition
            % If the block is too far IN THE PAST from the ripple, we should skip it

            refTime = mean(rippTimes(iRipp, :));

            % If we found a valid time:
            if t_start < refTime - maxDist
                % This block is too old for this ripple (and likely future ones)
                % Move main index forward
                currentFreeIdx = currentFreeIdx + 1;
                tempIdx = currentFreeIdx;
                continue;
            end

            ctrlTimes(iRipp, :) = [t_start, t_start + needed];
            lastEndTime = t_start + needed;
            break;

        else
            % This block is not suitable (too short remaining after t_start)
            % Move to next block
            tempIdx = tempIdx + 1;
        end
    end
end

end     % EOF
