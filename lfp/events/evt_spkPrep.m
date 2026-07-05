function [spkTimes, muTimes, uType] = evt_spkPrep(v, win, sigDur, fsSpk)
% EVT_SPKPREP Prepare single- and multi-unit spike times for event analysis.
%
%   [spkTimes, muTimes, uType] = EVT_SPKPREP(v, win, sigDur, fsSpk)
%
%   SUMMARY:
%       Shared spike-times preparation for both event pipelines (ripples, ED).
%       Pulls the optional spikes / spktimes / units from the basepaths2vars
%       struct v, shifts single-unit spikes into the window-relative frame and
%       clips to [0 sigDur], and flattens the raw (sample-based) spktimes into
%       one sorted MUA vector. Absent data yields empty outputs, so the
%       downstream spike analyses simply skip.
%
%   INPUTS:
%       v      - (Struct) basepaths2vars output. Reads the optional fields
%                         .spikes.times (cell, s), .spktimes (cell, samples),
%                         and .units.type (categorical); each may be missing.
%       win    - (Vec)    Analysis window [start end] (s).
%       sigDur - (Num)    Window duration (s); Inf for the full recording.
%       fsSpk  - (Num)    Spike sampling rate (Hz), to convert spktimes.
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
%       Updated: 05 Jul 2026 (takes the loader struct v; guards the fields here).

% Single units (already in seconds)
if isfield(v, 'spikes') && isfield(v.spikes, 'times') && ~isempty(v.spikes.times)
    spkTimes = v.spikes.times;
    spkTimes = cellfun(@(x) x - win(1), spkTimes, 'Uni', false);
    spkTimes = cellfun(@(x) x(x >= 0 & x <= sigDur), spkTimes, 'Uni', false);
else
    spkTimes = {};
end

% Multi-unit (samples -> seconds, pooled + sorted into one vector)
if isfield(v, 'spktimes') && ~isempty(v.spktimes)
    muTimes = cellfun(@(x) x / fsSpk, v.spktimes, 'uni', false);
    muTimes = cellfun(@(x) x - win(1), muTimes, 'Uni', false);
    muTimes = cellfun(@(x) x(x >= 0 & x <= sigDur), muTimes, 'Uni', false);
    muTimes = {sort(vertcat(muTimes{:}))};
else
    muTimes = {[]};
end

% Unit types (optional)
uType = [];
if isfield(v, 'units') && isfield(v.units, 'type')
    uType = v.units.type;
end

end
