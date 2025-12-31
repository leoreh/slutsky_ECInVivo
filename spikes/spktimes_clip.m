function tOut = spktimes_clip(tIn, bouts, recDur)
% SPKTIMES_CLIP Limits spike times to specific bouts and shifts times.
%
%   SPKTIMES = SPKTIMES_CLIP(SPKTIMES, BOUTS, RECDUR) restricts the given
%   cell array of spike times (or multi-unit times) to the intervals
%   defined in BOUTS. Time segments are then concatenated to strictly match
%   the duration RECDUR.
%
%   INPUTS:
%       spktimes    - (cell) Spike times per unit (in seconds) [nUnits x 1]
%       bouts       - (mat) Start and end times of bouts [nBouts x 2]
%       recDur      - (num) Target total duration (in seconds)
%
%   OUTPUTS:
%       spktimes    - (cell) Shifted spike times [nUnits x 1]

%% ========================================================================
%  CLIP BOUTS
%  ========================================================================

% Calculate durations
durs = bouts(:, 2) - bouts(:, 1);
cumDur = cumsum(durs);

% Find how many bouts needed
idxEnd = find(cumDur >= recDur, 1);

% Clip valid intervals to recDur
validBouts = bouts(1:idxEnd, :);

% Calculate how much overflow we have
overflow = cumDur(idxEnd) - recDur;
% Shorten the last bout by the overflow amount
validBouts(end, 2) = validBouts(end, 2) - overflow;

%% ========================================================================
%  SHIFT SPIKES
%  ========================================================================

% Recalculate durations for shifting
finalDurs = validBouts(:, 2) - validBouts(:, 1);
offsets = [0; cumsum(finalDurs(1:end-1))];

% Apply to spktimes
nUnits = length(tIn);
tOut = cell(nUnits, 1); 
nBouts = size(validBouts, 1);

for iUnit = 1:nUnits
    st = tIn{iUnit};
    newSpks = cell(nBouts, 1);

    for iBout = 1:nBouts
        t1 = validBouts(iBout, 1);
        t2 = validBouts(iBout, 2);
        offset = offsets(iBout);

        % Extract
        mask = st >= t1 & st < t2;
        chunk = st(mask);

        % Shift
        chunk = chunk - t1 + offset;

        % Append
        newSpks{iBout} = chunk(:);
    end
    tOut{iUnit} = cell2mat(newSpks);
end

end     % EOF
