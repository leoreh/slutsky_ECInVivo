function ripp = ripp_wrapper(varargin)
% RIPP_WRAPPER Master function to execute the complete SWR analysis pipeline.
%
%   ripp = RIPP_WRAPPER(varargin)
%
%   SUMMARY:
%       Coordinates the detection, analysis, and visualization of Sharp-Wave Ripples (SWR).
%       The pipeline proceeds in the following steps:
%       1.  Data Loading: Loads LFP, Spikes, and State data.
%       2.  Signal Prep: Filters LFP and calculates envelopes/z-scores (ripp_sigPrep).
%       3.  Detection: Identifies candidate events based on thresholds (ripp_times).
%       4.  Parameterization: Calculates stats like Amp/Freq/Energy (ripp_params).
%       5.  Maps: Generates peri-event LFP maps (ripp_maps).
%       6.  States: Classifies events by vigilance state (ripp_states).
%       7.  Spiking: Analyzes SU/MU modulation and generates PETHs (ripp_spks, ripp_spkPeth).
%       8.  Phasing: Calculates Spike-LFP coupling (spklfp_phase).
%       9.  Quality Assurance: Filters events based on spiking gain (optional).
%       10. Visualization: Runs the interactive GUI (ripp_gui) and saves plots.
%
%   INPUTS:
%       varargin - Parameter/Value pairs:
%           'basepath'   - (Char) Base directory of the session (default: pwd).
%           'rippCh'     - (Num)  Zero-indexed channel ID for ripple detection.
%                                 If empty, tries to load from session tags.
%           'thr'        - (Vec)  Detection thresholds [start, peak, cont, max, min_cont].
%           'limDur'     - (Vec)  Duration limits [min, max, inter, min_cont_dur] (ms).
%           'win'        - (Vec)  Time window to analyze [start end] (s). Default: [0 Inf].
%           'passband'   - (Vec)  Filtering frequency band [min max] (Hz). Default: [80 250].
%           'detectMet'  - (Num)  Detection method ID (see ripp_sigPrep). Default: 4.
%           'zMet'       - (Char) Z-scoring method ('adaptive', 'nrem'). Default: 'nrem'.
%           'mapDur'     - (Vec)  Window for PETH/Maps [pre post] (s). Default: [-0.05 0.05].
%           'bit2uv'     - (Num)  Conversion factor. Auto-detects Intan/OpenEphys if empty.
%           'flgPlot'    - (Log)  Generate summary plots? Default: true.
%           'flgSave'    - (Log)  Save output .mat files? Default: false.
%           'verbose'    - (Log)  Print progress steps? Default: true.
%           'steps'      - (Char) Execution mode: 'all', 'spks', 'phase'.
%                                 'all': Run full pipeline (default).
%                                 'spks': Skip sigLoad/detect. Calc Spks/PETH.
%                                 'phase': Skip detect/spks. Calc Phase.
%
%   OUTPUTS:
%       ripp         - (Struct) Combined structure containing all analysis results:
%           .times       - Event start/end times.
%           .peakTime    - Event peak times.
%           .amp,.freq   - Event parameters.
%           .state       - Vigilance state per event.
%           .spkGain     - MUA gain per event.
%
%   files saved (if flgSave=true):
%       basename.ripp.mat
%       basename.rippStates.mat
%       basename.rippSpks.mat
%       basename.rippPeth.mat
%       basename.rippSpkLfp.mat
%
%   DEPENDENCIES:
%       ripp_sigPrep, ripp_times, ripp_params, ripp_maps, ripp_states,
%       ripp_spks, ripp_spkPeth, ripp_plotSpks, spklfp_phase, ripp_gui.
%
%   HISTORY:
%       Updated: 23 Jan 2026

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'rippCh', [], @isnumeric);
addParameter(p, 'thr', [1, 3.5, 2, 200, 50], @isnumeric);
addParameter(p, 'limDur', [15, 300, 20, 10], @isnumeric);
addParameter(p, 'win', [0, Inf], @isnumeric);
addParameter(p, 'passband', [80 250], @isnumeric);
addParameter(p, 'detectMet', 3, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'bit2uv', [], @isnumeric);
addParameter(p, 'zMet', 'nrem', @ischar);
addParameter(p, 'mapDur', [-0.1 0.1], @isnumeric);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'steps', 'all', @ischar);

% Parse Inputs
parse(p, varargin{:});
basepath    = p.Results.basepath;
win         = p.Results.win;
flgPlot     = p.Results.flgPlot;
flgSave     = p.Results.flgSave;
mapDur      = p.Results.mapDur;
verbose     = p.Results.verbose;
steps       = p.Results.steps;

% Flags
doLoad   = ismember(steps, {'all', 'phase'});
doDetect = strcmpi(steps, 'all');
doSpks   = ismember(steps, {'all', 'spks'});
doPhase  = ismember(steps, {'all', 'phase'});

%% ========================================================================
%  SETUP
%  ========================================================================
cd(basepath);
[~, basename] = fileparts(basepath);

if verbose, fprintf('[RIPP]: Starting pipeline for %s...\n', basename); end
if verbose, fprintf('[RIPP]: Loading data...\n'); end

% Load Session & Data
vars = {'session', 'spikes', 'spktimes', 'sleep_states', 'units'};

if ~doDetect
    vars = [vars, {'ripp', 'rippMaps'}];
end

v = basepaths2vars('basepaths', {basepath}, 'vars', vars);

% Session Parameters
fs      = v.session.extracellular.srLfp;
fsSpk   = v.session.extracellular.sr;
uType   = v.units.type;

if isinf(win(2))
    sigDur = Inf;
else
    sigDur = win(2) - win(1);
end

% Prepare Multi-Unit (MU) Spike Times
% Flatten all units into a single sorted vector for MUA analysis
muTimes = cellfun(@(x) x / fsSpk, v.spktimes, 'uni', false);
muTimes = cellfun(@(x) x - win(1), muTimes, 'Uni', false);
muTimes = cellfun(@(x) x(x >= 0 & x <= sigDur), muTimes, 'Uni', false);
muTimes = {sort(vertcat(muTimes{:}))};

% Prepare Single-Unit (SU) Spike Times
spkTimes = v.spikes.times;
spkTimes = cellfun(@(x) x - win(1), spkTimes, 'Uni', false);
spkTimes = cellfun(@(x) x(x >= 0 & x <= sigDur), spkTimes, 'Uni', false);
nUnits = length(spkTimes);

% Prepare Bout Times
try
    boutTimes = v.ss.bouts.times;
    boutTimes = cellfun(@(x) x - win(1), boutTimes, 'UniformOutput', false);
    boutTimes = cellfun(@(x) x(x(:,2)>0 & x(:,1)<sigDur, :), boutTimes, 'UniformOutput', false);
    nremTimes = boutTimes{4}; % Assuming standard 4th cell is NREM
    vldTimes = vertcat(boutTimes{2}, boutTimes{3}, boutTimes{4});
catch
    warning('Could not load sleep states.');
    boutTimes = [];
    nremTimes = [];
    vldTimes = [];
end

% Output Files
files.ripp      = fullfile(basepath, [basename, '.ripp.mat']);
files.maps      = fullfile(basepath, [basename, '.rippMaps.mat']);
files.spks      = fullfile(basepath, [basename, '.rippSpks.mat']);
files.peth      = fullfile(basepath, [basename, '.rippPeth.mat']);
files.phase     = fullfile(basepath, [basename, '.rippSpkLfp.mat']);


%% ========================================================================
%  DETECT
%  ========================================================================
% Prepare Signals
rippSig = [];
if doLoad
    rippSig = ripp_sigLoad(p.Results, v.session, nremTimes);
end

% Detect
if doDetect
    ripp = ripp_detect(rippSig, p.Results, boutTimes, muTimes, fs, fsSpk);

    % Control Intervals
    ripp.ctrlTimes = ripp_ctrlTimes(ripp.times, ...
        'vldTimes', vldTimes, ...
        'flgPlot', false);

    % States (Final)
    [ripp.state, ~] = ripp_states(ripp.times, ...
        ripp.peakTime, ...
        boutTimes, ...
        'basepath', basepath, ...
        'flgPlot', flgPlot, ...
        'flgSave', true);

    % Params
    ripp = ripp_params(rippSig, ripp);

    % Maps
    if verbose, fprintf('[RIPP]: Generating Maps...\n'); end
    rippMaps = ripp_maps(rippSig, ripp.peakTime, fs, 'mapDur', mapDur, ...
        'flgSave', false);

else
    if verbose, fprintf('[RIPP]: Skipping Detection (loading from file)...\n'); end
    if isfield(v, 'ripp'), ripp = v.ripp; end
    if isfield(v, 'rippMaps'), rippMaps = v.rippMaps; end
end


%% ========================================================================
%  SPIKE STATS
%  ========================================================================
if doSpks
    if verbose, fprintf('[RIPP]: Calculating Spike Stats...\n'); end

    % Firing Metrics
    rippSpks = ripp_spks(spkTimes, ...
        ripp.times, ...
        ripp.ctrlTimes, ...
        ripp.peakTime, ...
        'unitType', uType, ...
        'flgSave', false);

    % Analyze single and burst FR metrics
    if isfield(v, 'brst')

        % Prepare single and burst spktimes
        brstTimes = v.brst.spktimes';
        brstTimes = cellfun(@(x) x - win(1), brstTimes, 'Uni', false);
        brstTimes = cellfun(@(x) x(x >= 0 & x <= sigDur), brstTimes, 'Uni', false);

        singleTimes = cell(1, nUnits);
        for iUnit = 1:nUnits
            singleTimes{iUnit} = setdiff(spkTimes{iUnit}, brstTimes{iUnit});
        end

        % Run ripp_spks
        rippSpks.burst = ripp_spks(brstTimes, ...
            ripp.times, ...
            ripp.ctrlTimes, ...
            ripp.peakTime, ...
            'flgSave', false);

        rippSpks.single = ripp_spks(singleTimes, ...
            ripp.times, ...
            ripp.ctrlTimes, ...
            ripp.peakTime, ...
            'flgSave', false);
    end

    % Add per ripple spike metrics to ripp struct
    ripp.spks = rippSpks.events;
    rippSpks = rmfield(rippSpks, 'events');
end


%% ========================================================================
%  PHASE COUPLING
%  ========================================================================

if doPhase

    % SPK-LFP
    % All Spikes
    if verbose, fprintf('[RIPP]: Calculating Spk-LFP Phase...\n'); end
    spkLfp = spklfp_phase(rippSig.filt, spkTimes, fs, ...
        'lfpTimes', ripp.times, ...
        'nPerms', 0);

    % % First Spikes
    % spkLfp.first = spklfp_phase(rippSig.filt, rippSpks.times.first, fs, ...
    %     'lfpTimes', ripp.times, ...
    %     'nPerms', 0);
    %
    % % Late Spikes
    % spkLfp.late = spklfp_phase(rippSig.filt, rippSpks.times.late, fs, ...
    %     'lfpTimes', ripp.times, ...
    %     'nPerms', 0);
end

% Remove times cell arrays after using them for phase coupling
if exist('rippSpks', 'var') && isfield(rippSpks, 'times')
    rippSpks = rmfield(rippSpks, 'times');
end


%% ========================================================================
%  SPIKES PETH
%  ========================================================================
if doSpks
    if verbose, fprintf('[RIPP]: Generating Spike PETHs...\n'); end

    % PETH
    rippPeth.su = ripp_spkPeth(spkTimes, ripp.peakTime, ripp.ctrlTimes, ...
        'flgSave', false, ...
        'mapDur', mapDur);
    rippPeth.mu = ripp_spkPeth(muTimes, ripp.peakTime, ripp.ctrlTimes, ...
        'flgSave', false, ...
        'mapDur', mapDur);

    tstamps = rippPeth.su.tstamps;

    % Smoothing Kernel
    dt = mode(diff(tstamps));
    sigma = 0.001;
    nSteps = ceil(3*sigma / dt);
    kRng = (-nSteps : nSteps) * dt;
    kd = normpdf(kRng, 0, sigma);
    kd = kd / sum(kd);

    % Unit PETH
    % ---------------------------------------------------------------------
    % Mean PETH across ripples (Units x Bins)
    pethSU = squeeze(mean(rippPeth.su.ripp, 2, 'omitnan'));

    % Control Statistics
    ctrlPeth = squeeze(mean(rippPeth.su.ctrl, 2, 'omitnan')); % (Units x Bins)

    % Calculate
    pethZ = peth_norm(pethSU, ctrlPeth, kd);

    % Add to rippSpks
    rippSpks.peth = pethZ;
    rippSpks.tstamps = tstamps;

    % Population PETH
    % ---------------------------------------------------------------------
    popTypes = {'RS', 'FS'};
    for iType = 1:length(popTypes)
        currType = popTypes{iType};
        idxType = uType == currType;

        if sum(idxType) == 0
            continue;
        end

        % Raw Population PETH (Ripples x Bins)
        % Sum spikes across units, divide by N units
        subRipp = rippPeth.su.ripp(idxType, :, :); % (Units x Ripples x Bins)
        popRipp = squeeze(sum(subRipp, 1)) ./ sum(idxType);

        % Control Statistics (Scalar)
        subCtrl = rippPeth.su.ctrl(idxType, :, :); % (Units x Controls x Bins)
        popCtrl = squeeze(sum(subCtrl, 1)) ./ sum(idxType); % (Controls x Bins)

        meanPopCtrl = mean(popCtrl, 1, 'omitnan'); % (1 x Bins)

        % Calculate
        pethZ = peth_norm(popRipp, meanPopCtrl, kd);

        % Store
        rippMaps.peth.(currType) = pethZ;
    end

    % MU PETH
    % ---------------------------------------------------------------------
    popRipp = squeeze(rippPeth.mu.ripp); % (Ripples x Bins)

    % Control Statistics
    popCtrl = squeeze(rippPeth.mu.ctrl); % (Controls x Bins)
    meanPopCtrl = mean(popCtrl, 1, 'omitnan'); % (1 x Bins)

    % Calculate
    pethZ = peth_norm(popRipp, meanPopCtrl, kd);

    % Store
    rippMaps.peth.MU = pethZ;
end

%% ========================================================================
%  SAVE
%  ========================================================================

% Revert times to absolute for saving/Neuroscope
ripp.times = ripp.times + win(1);
ripp.peakTime = ripp.peakTime + win(1);

if flgSave
    if verbose, fprintf('[RIPP]: Saving output files...\n'); end

    % RIPP & MAPS
    if doDetect || doSpks
        save(files.ripp, 'ripp', '-v7.3');
        save(files.maps, 'rippMaps', '-v7.3');
    end

    % SPKS
    if doSpks
        save(files.spks, 'rippSpks', '-v7.3');
        save(files.peth, 'rippPeth', '-v7.3');
    end

    % PHASE
    if doPhase
        save(files.phase, 'spkLfp', '-v7.3');
    end

end


%% ========================================================================
%  PLOT
%  ========================================================================

% Plot spikes
if flgPlot
    if verbose, fprintf('[RIPP]: Generating Summary Plot...\n'); end
    ripp_plotSpks(rippSpks, rippPeth, ...
        'basepath', basepath, ...
        'flgSaveFig', true);
end

if flgPlot

    % Per Ripple Map
    tblMap = struct2table(rmfield(rippMaps, {'tstamps', 'peth'}));

    % Add PETH
    tblMap.pethRs = rippMaps.peth.RS;
    tblMap.pethFs = rippMaps.peth.FS;
    tblMap.pethMu = rippMaps.peth.MU;

    % Plot
    tblGUI_xy(rippMaps.tstamps, tblMap, 'yVar', 'z', 'grpVar', 'states');

    % Per Unit Map
    tblUnit = table();
    tblUnit.unitType = uType;
    tblUnit.peth = rippSpks.peth;
    tblGUI_xy(rippPeth.su.tstamps, tblUnit, 'yVar', 'peth', 'grpVar', 'states');

end

if verbose, fprintf('[RIPP]: Pipeline Completed for %s.\n', basename); end


end     % EOF


%% ========================================================================
%  HELPER: RIPP_SIGLOAD
%  ========================================================================

function rippSig = ripp_sigLoad(params, session, nremTimes)
% RIPP_SIGLOAD Helper to load and prepare LFP signal for ripple analysis.
%
%   rippSig = ripp_sigLoad(params, session, nremTimes)
%
%   INPUTS:
%       params    - (Struct) Parsed inputs from ripp_wrapper (p.Results).
%       session   - (Struct) Session metadata (v.session).
%       nremTimes - (M x 2) NREM bout times (s).
%
%   OUTPUTS:
%       rippSig   - (Struct) Prepared ripple signals (raw, filt, envelope, etc).

% Unpack Parameters
basepath = params.basepath;
[~, basename] = fileparts(basepath);
win = params.win;
verbose = params.verbose;

% Setup
if verbose, fprintf('[RIPP]: Pre-processing signals...\n'); end

% Duration
if isinf(win(2))
    sigDur = Inf;
else
    sigDur = win(2) - win(1);
end

% Session Info
fs = session.extracellular.srLfp;
nchans = session.extracellular.nChannels;

% Ripple Channel
if isempty(params.rippCh)
    if isfield(session.channelTags, 'Ripple')
        rippCh = session.channelTags.Ripple;
    else
        rippCh = 1;
    end
else
    rippCh = params.rippCh;
end

% Voltage Conversion
bit2uv = params.bit2uv;
if isempty(bit2uv)
    if round(session.extracellular.sr) == 24414
        bit2uv = 1;     % TDT/Tucker
    else
        bit2uv = 0.195; % Intan
    end
end

% Load LFP
fname = fullfile(basepath, [basename, '.lfp']);
lfp = binary_load(fname, 'duration', sigDur, 'fs', fs, 'nCh', nchans,...
    'start', win(1), 'ch', rippCh, 'downsample', 1, 'bit2uv', bit2uv);

if size(lfp, 2) > 1
    lfp = mean(lfp, 2);
end

% Prepare Signals
rippSig = ripp_sigPrep(lfp, fs, ...
    'detectMet', params.detectMet, ...
    'passband', params.passband, ...
    'zMet', params.zMet, ...
    'nremTimes', nremTimes);

end


%% ========================================================================
%  HELPER: RIPP_DETECT
%  ========================================================================

function ripp = ripp_detect(rippSig, params, boutTimes, muTimes, fs, fsSpk)
% RIPP_DETECT Helper to detect ripples.
%
%   INPUTS:
%       rippSig   - (Struct) Prepared ripple signals (from ripp_sigLoad).
%       params    - (Struct) Use parameters (p.Results).
%       boutTimes - (Cell)  Sleep state bout times.
%       muTimes   - (Cell)  Multi-unit spike times.
%       fs        - (Num)   Sampling rate.
%       fsSpk     - (Num)   Spike sampling rate (for Neuroscope).
%
%   OUTPUTS:
%       ripp      - (Struct) Ripple events and stats.


% Unpack Params
basepath = params.basepath;
[~, basename] = fileparts(basepath);
win = params.win;
thr = params.thr;
limDur = params.limDur;
flgPlot = params.flgPlot;
flgSave = params.flgSave;
verbose = params.verbose;

% Load EMG
% -------------------------------------------------------------------------
if verbose, fprintf('[RIPP]: Loading EMG...\n'); end
load(fullfile(basepath, [basename, '.sleep_sig.mat']), 'emg');

% Trim EMG
s1 = round(win(1) * fs) + 1;
if isinf(win(2))
    emg = emg(s1:end);
else
    s2 = min(length(emg), round(win(2) * fs));
    emg = emg(s1:s2);
end

if length(emg) ~= length(rippSig.lfp)
    error('EMG and LFP do not fit')
end


% Detect
% -------------------------------------------------------------------------
if verbose, fprintf('[RIPP]: Thresholding candidate events...\n'); end

% Detect Candidate Events
ripp = ripp_times(rippSig, fs, ...
    'thr', thr, ...
    'limDur', limDur);

% States (Prelim)
ripp.state = ripp_states(ripp.times, ripp.peakTime, boutTimes, ...
    'basepath', basepath, ...
    'flgPlot', false, ...
    'flgSave', false);

% Quality Assurance
[idxQA, met] = ripp_qa(ripp.times, boutTimes{4}, muTimes, emg, ripp.state, fs);

% Store
ripp.emg = met.emg;
ripp.spkGain = met.gain;
ripp.idxQA = idxQA;

% Select "Good" Ripples
idxGood = idxQA.good;

if verbose
    fprintf('[RIPP]: QA Kept %d / %d events.\n', sum(idxGood), length(idxGood));
end

% Filter Struct
fnames = fieldnames(ripp);
for iField = 1:length(fnames)
    fn = fnames{iField};
    if size(ripp.(fn), 1) == length(idxGood)
        ripp.(fn) = ripp.(fn)(idxGood, :);
    end
end

% GUI
% -------------------------------------------------------------------------
if flgPlot
    if verbose, fprintf('[RIPP]: Launching Detection GUI...\n'); end
    ripp_gui(ripp.times, ripp.peakTime, rippSig, muTimes, fs, thr, emg, ...
        'basepath', basepath);
end

% Neuroscope Saving
% -------------------------------------------------------------------------
if flgSave
    if verbose, fprintf('[RIPP]: Saving Neuroscope events...\n'); end
    % Convert to samples add win(1) back to get absolute time
    absRipp = ripp.times + win(1);
    absPeak = ripp.peakTime + win(1);

    rippSamps = round(absRipp * fsSpk);
    peakSamps = round(absPeak * fsSpk);

    ripp2ns(rippSamps, peakSamps, 'basepath', basepath);
end
end


%% ========================================================================
%  HELPER: PETH_NORM
%  ========================================================================

function pethZ = peth_norm(pethRaw, refData, kd)
% peth_norm Helper function to calculate smooth, normalized PETH stats.
%
%   pethZ = peth_norm(pethRaw, refData, kd)
%
%   INPUTS:
%       pethRaw - (M x N) Raw PETH matrix (units/events x bins).
%       refData - (M x N) or (1 x N) Reference PETH data.
%       kd      - (1 x K) Smoothing kernel.
%
%   OUTPUTS:
%       pethZ   - (M x N) Z-scored PETH.

% Edge Effect Correction Vector
% Normalized kernel convolution with unity vector reveals boundary loss
nBins = size(pethRaw, 2);
kCorr = conv(ones(1, nBins), kd, 'same');

% Smooth Target
pethSmooth = conv2(pethRaw, kd, 'same') ./ kCorr;

% Reference
if size(refData, 1) == size(pethRaw, 1) && size(refData, 1) > 1
    % Unit-wise normalization (Reference is Units x Bins)
    refMean = mean(refData, 2, 'omitnan'); % (Units x 1)
    refSd = std(refData, [], 2, 'omitnan'); % (Units x 1)
else
    % Independent Reference (e.g. Population Mean Vector 1 x Bins)
    refMean = mean(refData, 'all', 'omitnan');
    refSd = std(refData, [], 'all', 'omitnan');
end

if any(refSd == 0)
    refSd(refSd == 0) = 1;
end

% Z-Score (Control Rates)
pethZ = (pethSmooth - refMean) ./ refSd;

end


%% ========================================================================
%  HELPER: RIPP_QA
%  ========================================================================

function [idxQA, met] = ripp_qa(rippTimes, nremTimes, muTimes, emg, rippState, fs)
% RIPP_QA Helper function to perform quality assurance on ripples.
%
%   [idxQA, met] = ripp_qa(rippTimes, nremTimes, muTimes, emg, rippState, fs)
%
%   INPUTS:
%       rippTimes - (N x 2) Event start/end times (s).
%       nremTimes - (M x 2) NREM bout times for baseline (s).
%       muTimes   - (Cell)  Multi-unit spike times {sorted_vec}.
%       emg       - (Vec)   EMG signal.
%       rippState - (Cell)  State labels (N x 1).
%       fs        - (Num)   Sampling rate (Hz).
%
%   OUTPUTS:
%       idxQA     - (Struct) Logical indices for filtering:
%                   .state  - Valid states (NREM, LSLEEP, QWAKE).
%                   .gain   - Positive spike gain.
%                   .emg    - Low EMG activity.
%                   .good   - Intersection of all above.
%       met       - (Struct) QA Metrics:
%                   .emg    - Z-scored EMG per ripple.
%                   .gain   - Spike gain per ripple.

% EMG Analysis
% -------------------------------------------------------------------------
emgAbs = abs(emg);
nRipp = size(rippTimes, 1);
rippEmg = nan(nRipp, 1);

for iRipp = 1:nRipp
    bStart = max(1, round(rippTimes(iRipp,1) * fs));
    bEnd = min(length(emg), round(rippTimes(iRipp,2) * fs));
    rippEmg(iRipp) = mean(emgAbs(bStart : bEnd));
end

% NREM Baseline
nremMask = false(size(emg));
if ~isempty(nremTimes)
    for iBout = 1:size(nremTimes, 1)
        bStart = max(1, round(nremTimes(iBout,1) * fs));
        bEnd = min(length(emg), round(nremTimes(iBout,2) * fs));
        nremMask(bStart : bEnd) = true;
    end
end
nremAmp = emgAbs(nremMask);

% Z-Score
if ~isempty(nremAmp)
    rippEmgZ = (rippEmg - mean(nremAmp, 'omitnan')) / std(nremAmp, 'omitnan');
else
    rippEmgZ = nan(size(rippEmg));
end

idxEmg = rippEmgZ < 2;

% Spike Gain
% -------------------------------------------------------------------------
% Generate control intervals locally
ctrlTimes = ripp_ctrlTimes(rippTimes);

% Calculate Rates
% muTimes is expected to be {vector}
rippRates = times2rate(muTimes, 'winCalc', rippTimes, 'binsize', Inf);
ctrlRates = times2rate(muTimes, 'winCalc', ctrlTimes, 'binsize', Inf);

% Gain = (Ripple Rate - Mean Control) / Std Control
muCtrl = mean(ctrlRates, 'all', 'omitnan');
sdCtrl = std(ctrlRates, [], 'all', 'omitnan');
if sdCtrl == 0, sdCtrl = 1; end

spkGain = (rippRates - muCtrl) ./ sdCtrl;
% Ensure gain is column vector for consistency
if size(spkGain, 2) > size(spkGain, 1), spkGain = spkGain'; end

idxGain = spkGain > 0;

% States
% -------------------------------------------------------------------------
% Valid states
validStates = {'QWAKE', 'LSLEEP', 'NREM'};
idxState = ismember(rippState, validStates);

% Combine
% -------------------------------------------------------------------------
idxGood = idxState & idxGain & idxEmg;

% Output
idxQA.state = idxState;
idxQA.gain = idxGain;
idxQA.emg = idxEmg;
idxQA.good = idxGood;

met.emg = rippEmgZ;
met.gain = spkGain;

end
