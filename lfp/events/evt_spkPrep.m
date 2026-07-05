function [spkTimes, muTimes, uType] = evt_spkPrep(spikes, spktimes, units, win, sigDur, fsSpk)
% EVT_SPKPREP Prepare single- and multi-unit spike times for event analysis.
%
%   [spkTimes, muTimes, uType] = EVT_SPKPREP(spikes, spktimes, units, win, sigDur, fsSpk)
%
%   SUMMARY:
%       Shared spike-times preparation for both event pipelines (ripples, ED).
%       Shifts spikes into the window-relative frame, clips to [0 sigDur], and
%       flattens the raw (sample-based) spktimes into one sorted MUA vector.
%       Every input is optional; absent data yields empty outputs so that the
%       downstream spike analyses simply skip.
%
%   INPUTS:
%       spikes   - (Struct) Sorted spikes with .times (cell, seconds), or [].
%       spktimes - (Cell)   {N x 1} Raw spike times in SAMPLES, or [].
%       units    - (Struct)  Unit classification with .type, or [].
%       win      - (Vec)     Analysis window [start end] (s).
%       sigDur   - (Num)     Window duration (s); Inf for the full recording.
%       fsSpk    - (Num)     Spike sampling rate (Hz), to convert spktimes.
%
%   OUTPUTS:
%       spkTimes - (Cell) {N_units x 1} SU spike times (s, relative). {} if absent.
%       muTimes  - (Cell) {1 x 1} Sorted pooled MUA vector (s, relative). {[]} if absent.
%       uType    - (Cat)  Unit types, or [] if absent.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 05 Jul 2026 (absorbs the duplicated wrapper cellfun chains).

% Single units (already in seconds)
if ~isempty(spikes) && isfield(spikes, 'times') && ~isempty(spikes.times)
    spkTimes = spikes.times;
    spkTimes = cellfun(@(x) x - win(1), spkTimes, 'Uni', false);
    spkTimes = cellfun(@(x) x(x >= 0 & x <= sigDur), spkTimes, 'Uni', false);
else
    spkTimes = {};
end

% Multi-unit (samples -> seconds, pooled + sorted into one vector)
if ~isempty(spktimes)
    muTimes = cellfun(@(x) x / fsSpk, spktimes, 'uni', false);
    muTimes = cellfun(@(x) x - win(1), muTimes, 'Uni', false);
    muTimes = cellfun(@(x) x(x >= 0 & x <= sigDur), muTimes, 'Uni', false);
    muTimes = {sort(vertcat(muTimes{:}))};
else
    muTimes = {[]};
end

% Unit types (optional)
uType = [];
if ~isempty(units) && isfield(units, 'type')
    uType = units.type;
end

end
