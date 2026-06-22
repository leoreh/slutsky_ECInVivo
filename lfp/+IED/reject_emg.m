function ied = reject_emg(ied, varargin)
% reject IED detections that coincide with high EMG (muscle contamination).
%
%   INPUT (1st):
%       ied     IED.data object after detection
%   INPUT (name-value):
%       method      'zscore' (default) or 'emg_rms'
%       emg         EMG trace [same fs/length as detection sig]. default: ied EMG channel
%       emg_rms     1-Hz log-RMS vector (e.g. sSig.emg_rms), for method 'emg_rms'
%       winSec      window around each event to measure EMG [s] {0.05}
%       baselineSec moving baseline for z-score [s] {5}
%       thrZ        reject if local EMG z-score exceeds this {3}
%       thrRms      reject if emg_rms exceeds this; default = 75th percentile
%       onlyAccepted  if true, only test currently accepted events {true}
%
%   OUTPUT:
%       ied with ied.accepted set false for contaminated events
%
%   see also IED.detect_move_z, IED.data

if ~isa(ied, 'IED.data')
    error("first input must be IED.data obj")
end

p = inputParser;
addParameter(p, 'method', 'zscore', @(x) ismember(x, {'zscore', 'emg_rms'}));
addParameter(p, 'emg', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'emg_rms', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'winSec', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'baselineSec', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'thrZ', 3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'thrRms', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'onlyAccepted', true, @islogical);
parse(p, varargin{:})

pos = ied.pos(:);
test_mask = true(size(pos));
if p.Results.onlyAccepted
    test_mask = ied.accepted(:);
end

switch p.Results.method
    case 'zscore'
        emg = p.Results.emg;
        if isempty(emg)
            if ~ismember("EMG", string(ied.data_sources.Properties.RowNames))
                error('reject_emg:no_emg', 'No EMG channel in ied; pass emg or use method ''emg_rms''.')
            end
            emg = ied.data_sources.data{"EMG"};
        end
        emg = emg(:);
        if numel(emg) ~= numel(ied.sig)
            error('reject_emg:length_mismatch', 'EMG must match detection signal length.')
        end

        win_samp = max(1, round(ied.fs * p.Results.winSec));
        base_win = max(1, round(ied.fs * p.Results.baselineSec));
        emg_abs = abs(emg);
        emg_mu = movmean(emg_abs, base_win);
        emg_sd = movstd(emg_abs, base_win);
        emg_sd(emg_sd == 0) = eps;

        peak_emg = zeros(numel(pos), 1);
        for iEv = 1:numel(pos)
            idx = pos(iEv);
            seg = emg_abs(max(1, idx - win_samp) : min(numel(emg_abs), idx + win_samp));
            peak_emg(iEv) = max(seg);
        end
        score = (peak_emg - emg_mu(pos)) ./ emg_sd(pos);
        contaminated = score > p.Results.thrZ;
        thr_msg = sprintf('thrZ=%g', p.Results.thrZ);

    case 'emg_rms'
        emg_rms = p.Results.emg_rms(:)';
        if isempty(emg_rms)
            error('reject_emg:no_emg_rms', 'Pass emg_rms (e.g. sSig.emg_rms) for method ''emg_rms''.')
        end
        if isempty(p.Results.thrRms)
            thr_rms = prctile(emg_rms, 75);
        else
            thr_rms = p.Results.thrRms;
        end
        event_sec = pos / ied.fs;
        score = interp1((1:numel(emg_rms))', emg_rms(:), event_sec, 'nearest', 'extrap');
        contaminated = score > thr_rms;
        thr_msg = sprintf('thrRms=%.3f', thr_rms);
end

reject = test_mask & contaminated;
ied.accepted(reject) = false;
fprintf('reject_emg (%s): rejected %d / %d events (%s)\n', ...
    p.Results.method, sum(reject), numel(pos), thr_msg);

end
