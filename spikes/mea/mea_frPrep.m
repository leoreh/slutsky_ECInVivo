function fr = mea_frPrep(spktimes, varargin)
% MEA_FRPREP Prepares firing rate matrix for MEA analysis.
%
%   fr = MEA_FRPREP(SPKTIMES) processes spike times to generate smoothed
%   firing rate matrices, detect network perturbations, and align the time
%   axis.
%
%   INPUTS:
%       spktimes    - (cell) Cell array of spike times for each unit.
%
%   OPTIONAL (Key-Value Pairs):
%       binSize     - (double) Bin size in seconds. {60}
%       winLim      - (vector) [start end] window in seconds. {[0 Inf]}
%       spkThr      - (double) Minimum spikes required per unit. {4}
%       sgPolyOrder - (double) Polynomial order for denoising. {3}
%       sgFrameSec  - (double) Frame length for denoising in seconds. {600}
%       flgSave     - (logical) Whether to save the output struct. {false}
%       flgPlot     - (logical) Whether to plot perturbation detection. {false}
%       basepath    - (char) Base path for saving files. {pwd}
%
%   OUTPUTS:
%       fr          - (struct) Structure containing:
%                       .fr         - Smoothed firing rate matrix.
%                       .frOrig     - Raw firing rate matrix.
%                       .t          - Corrected time vector [seconds].
%                       .uGood      - Logical vector of good units.
%                       .info       - Analysis parameters and metadata.
%
%   See also: TIMES2RATE, FR_DENOISE, MCU_DETECTPERT

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'binSize', 60, @isnumeric);
addParameter(p, 'winLim', [0 Inf], @isnumeric);
addParameter(p, 'spkThr', 4, @isnumeric);
addParameter(p, 'sgPolyOrder', 3, @isnumeric);
addParameter(p, 'sgFrameSec', 600, @isnumeric);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'basepath', pwd, @ischar);

parse(p, spktimes, varargin{:});
binSize = p.Results.binSize;
winLim = p.Results.winLim;
spkThr = p.Results.spkThr;
sgPolyOrder = p.Results.sgPolyOrder;
sgFrameSec = p.Results.sgFrameSec;
flgSave = p.Results.flgSave;
flgPlot = p.Results.flgPlot;
basepath = p.Results.basepath;


%% ========================================================================
%  FILTER AND CALCULATE RATES
%  ========================================================================

% Filter Spikes
if isinf(winLim(2))
    if isempty(spktimes)
        maxSpk = 0;
    else
        maxSpk = max(cellfun(@(x) max(x, [], 'omitnan'), spktimes));
    end
    winLim(2) = maxSpk;
end
spktimes = cellfun(@(x) x(x >= winLim(1) & x <= winLim(2)), spktimes, 'UniformOutput', false);

% Identify good units
if ~isempty(spktimes)
    nSpks = cellfun(@length, spktimes)';
else
    nSpks = [];
end
uGood = nSpks >= spkThr;

% Calculate Raw Firing Rates
[frOrig, ~, ~] = times2rate(spktimes, 'binsize', binSize, 'c2r', true);


%% ========================================================================
%  DETECT PERTURBATION
%  ========================================================================

% Use original data for detection, masking bad units
[idxPert, tAxis] = mcu_detectPert(frOrig(uGood, :), 'binSize', binSize, ...
    'srchStart', 1, 'srchEnd', 3, 'flgPlot', flgPlot, 'frameLen', 10);


%% ========================================================================
%  DENOISE
%  ========================================================================

sgFrameLen = round(sgFrameSec / binSize);
frSm = fr_denoise(frOrig, 'flgPlot', false, 'frameLen', sgFrameLen);


%% ========================================================================
%  CORRECT TIME AXIS
%  ========================================================================

% Create 0-centered index vector
nBins = size(frSm, 2);
tIdx = (1:nBins) - idxPert;

% Convert to initial physical time (seconds) assuming continuous
tSec = tIdx * binSize;

% Apply corrections based on recording protocol:
% - Baseline: 20 min saved every hour (3x scaling)
% - Post-Perturbation: 20 min saved every two hours (6x scaling)
t = zeros(size(tSec));

% Pre-perturbation
t(tIdx <= 0) = tSec(tIdx <= 0) * 3;

% Post-perturbation (including 0)
t(tIdx > 0) = tSec(tIdx > 0) * 6;


%% ========================================================================
%  OUTPUT STRUCT
%  ========================================================================

fr.fr = frSm;
fr.frOrig = frOrig;
fr.t = t;
fr.uGood = uGood;
fr.info.idxPert = idxPert;
fr.info.binSize = binSize;
fr.info.winLim = winLim;
fr.info.spkThr = spkThr;
fr.info.sgPolyOrder = sgPolyOrder;
fr.info.sgFrameSec = sgFrameSec;
fr.info.correctionNote = 'Time scaled x3 (pre) and x6 (post) relative to idxPert';

if flgSave
    [~, basename] = fileparts(basepath);
    save(fullfile(basepath, [basename, '.fr.mat']), 'fr', '-v7.3');
end

end
