function ripp = ripp_stats(rippSig, ripp, fs, varargin)
% RIPP_STATS Calculates statistics and maps for detected ripples.
%
% SUMMARY:
% Computes peri-ripple maps, ACG, estimates rate, and calculates additional
% properties like peak amplitude, frequency, and energy.
%
% INPUT:
%   rippSig         Structure with fields: lfp, filt, amp, freq, z.
%   ripp            Structure with detected events (must have .times, .peakTime).
%   fs              Sampling frequency [Hz].
%   varargin        'basepath' (for saving/plots), 'flgPlot', 'flgSave'.
%                   'thr' (optional, included for compatibility/info).
%
% OUTPUT:
%   ripp            Updated structure with .maps, .acg, .rate, .info, .peakAmp etc.
%
% DEPENDENCIES: Sync, SyncMap, CCG, times2rate, ripp_plot.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'rippSig', @isstruct);
addRequired(p, 'ripp', @isstruct);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'thr', [1.5, 2.5, 2, 200, 50], @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', true, @islogical);

parse(p, rippSig, ripp, fs, varargin{:});
rippSig     = p.Results.rippSig;
ripp        = p.Results.ripp;
fs          = p.Results.fs;
basepath    = p.Results.basepath;
thr         = p.Results.thr;
flgPlot     = p.Results.flgPlot;
flgSave     = p.Results.flgSave;

%% ========================================================================
%  SETUP
%  ========================================================================
% Initialize additional fields if missing
if ~isfield(ripp, 'maps'), ripp.maps = struct('durWin', [-0.075 0.075]); end
if ~isfield(ripp, 'info'), ripp.info = struct(); end
if ~isfield(ripp, 'rate'), ripp.rate = struct(); end
if ~isfield(ripp, 'acg'), ripp.acg = struct(); end

% Define file for saving
[~, basename] = fileparts(basepath);
rippfile = fullfile(basepath, [basename, '.ripp.mat']);

timestamps = (0:length(rippSig.lfp)-1)' / fs;
ripp.rate.binsize = 60;

% If no ripples, return early but save empty implementation if requested
if isempty(ripp.times)
    warning('No ripples passed to ripp_stats.');
    if flgSave, save(rippfile, 'ripp', '-v7.3'); end
    return;
end

peakTime = ripp.peakTime;
times = ripp.times;
nEvents = length(peakTime);

% Calculate durations if not present
if ~isfield(ripp, 'dur') || isempty(ripp.dur)
    ripp.dur = (times(:,2) - times(:,1)) * 1000; % stored in ms
end


%% ========================================================================
%  PROPERTIES & STATISTICS
%  ========================================================================

% Maps configuration
nbinsMap = floor(fs * diff(ripp.maps.durWin) / 2) * 2 + 1;
centerBin = ceil(nbinsMap / 2);
mapArgs = {'durations', ripp.maps.durWin, 'nbins', nbinsMap, 'smooth', 0};

% 1. Filtered Signal Map
[r, ir] = Sync([timestamps rippSig.filt], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.raw = SyncMap(r, ir, mapArgs{:}); % Note: original code called it .raw but used .filt input
ripp.maps.filt = ripp.maps.raw;

% 2. Amplitude Envelope Map
[r, ir] = Sync([timestamps rippSig.amp], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.amp = SyncMap(r, ir, mapArgs{:});

% 3. Phase Map (if available)
if isfield(rippSig, 'sigPhase')
    [r, ir] = Sync([timestamps rippSig.sigPhase], peakTime, 'durations', ripp.maps.durWin);
    ripp.maps.phase = SyncMap(r, ir, mapArgs{:});
end

% 4. Frequency Map (if available)
if isfield(rippSig, 'freq')
    [r, ir] = Sync([timestamps rippSig.freq], peakTime, 'durations', ripp.maps.durWin);
    ripp.maps.freq = SyncMap(r, ir, mapArgs{:});
end

% Store Peak Properties
ripp.peakAmp = ripp.maps.amp(:, centerBin);
if isfield(ripp.maps, 'freq')
    ripp.peakFreq = ripp.maps.freq(:, centerBin);
end
ripp.peakFilt = ripp.maps.filt(:, centerBin);

idxWin = centerBin + (-5:5);
ripp.peakEnergy = mean(ripp.maps.filt(:, idxWin).^2, 2);

% ACG
[ripp.acg.data, ripp.acg.t] = CCG(peakTime, ones(nEvents, 1), 'binSize', 0.01, 'duration', 1);

% Rate
[ripp.rate.rate, ripp.rate.binedges, ripp.rate.timestamps] = ...
    times2rate(peakTime, 'binsize', ripp.rate.binsize, 'winCalc', [timestamps(1) timestamps(end)], 'c2r', true);


%% ========================================================================
%  PARAMETER SUGGESTIONS (Legacy Support)
%  ========================================================================
% Re-calculating some stats needed for "thrData" suggestions if not passed in
% To perfectly replicate ripp_detect's suggestions we would need the raw detection vectors.
% However, we can approximate or skip if those internal vars aren't available.
% The original code computed percentiles of 'contPowAvg' and 'peakPower'.
% We can reconstruct peakPower from `rippSig.z` at `peakTime`.

% Re-extract peak power from Z-score for stats
z_at_peak = interp1(timestamps, rippSig.z, peakTime, 'nearest');
% Note: this is an approximation since peakTime is refined.
% Better to just use the peak amplitude we found? No, detection used Z-score.
% Let's grab the Z maps to be accurate.
[r_z, ir_z] = Sync([timestamps rippSig.z], peakTime, 'durations', ripp.maps.durWin);
z_map = SyncMap(r_z, ir_z, mapArgs{:});
peakPower = z_map(:, centerBin);




%% ========================================================================
%  FINALIZE
%  ========================================================================

if flgSave
    save(rippfile, 'ripp', '-v7.3');
end

if flgPlot
    % Call the plotting function if available
    if exist('ripp_plot', 'file')
        ripp_plot(ripp, 'basepath', basepath, 'flgSaveFig', true);
    else
        warning('ripp_plot function not found. Skipping plot.');
    end
end

end
