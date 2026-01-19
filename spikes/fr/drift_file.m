function drft = drift_file(spktimes, varargin)
% DRIFT_FILE Wrapper for calculating population drift over time.
%
%   drft = DRIFT_FILE(SPKTIMES, ...) calculates population vector drift
%   for the given spike times, potentially split into chunks if a window
%   size is specified.
%
%   INPUTS:
%       spktimes    - (cell) Spike times per unit.
%       varargin    - (param/value) Optional parameters:
%                     'winLim'     : (num) [Start, End] time limit for analysis.
%                                    Default: [0, max(spktimes)].
%                     'winSize'    : (num) Size of chunks to split winLim into.
%                                    If empty, uses reduced winLim as one chunk.
%                     'binSize'    : (num) Bin size for drift vectors {1200} (s).
%                     'basepath'   : (char) Base path for saving {pwd}
%                     'flgSave'    : (log) Save result as drift struct {false}
%                     'flgPlot'    : (log) Plot drift for each chunk {false}
%
%   OUTPUTS:
%       drft        - (struct) Drift statistics structure with concatenated fields:
%                     .time     : (nChunks x 2) Time windows used.
%                     .dt_corr  : (nChunks x nLags) Cell/Matrix of correlations.
%                     .drate    : (nChunks x 1) Drift rates.
%                     ... etc.
%
%   See also: DRIFT_CALC, FR_NETWORK, N2CHUNKS, CATFIELDS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'winLim', [], @isnumeric);
addParameter(p, 'winSize', [], @isnumeric);
addParameter(p, 'binSize', 1200, @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'flgPlot', false, @islogical);

parse(p, spktimes, varargin{:});
winLim   = p.Results.winLim;
winSize  = p.Results.winSize;
binSize  = p.Results.binSize;
basepath = p.Results.basepath;
flgSave  = p.Results.flgSave;
flgPlot  = p.Results.flgPlot;

%% ========================================================================
%  INITIALIZE
%  ========================================================================

% Hardcoded Drift Parameters (passed to drift_calc)
params.thrFr   = 0.005;
params.thrLin  = 20;
params.limUnit = [];

% Handle winLim
maxTime = max(cellfun(@(x) max([0; x(:)]), spktimes));
if isempty(winLim)
    winLim = [0, maxTime];
end
if winLim(2) > maxTime
    winLim(2) = maxTime;
end

% Handle winSize / Chunks
if ~isempty(winSize)
    % Split into chunks, excluding last incomplete one
    chunks = n2chunks('n', winLim(2), 'chunksize', winSize, ...
        'lastChunk', 'exclude', 'clip', [0, winLim(1); winLim(2), Inf]);
else
    chunks = winLim;
end

nChunks = size(chunks, 1);
nUnits = length(spktimes);

% Info Struct
drft.info.input   = p.Results;
drft.info.params  = params;
drft.info.nUnits  = nUnits;
drft.win = chunks;

if nUnits < 2
    warning('Fewer than 2 units provided. Returning empty struct.');
    return;
end

%% ========================================================================
%  COMPUTE LOOP
%  ========================================================================

statsArr = struct('dt_corr', {}, 'm_corr', {}, 'lin_coef', {}, 'drate', {});

for iChunk = 1:nChunks

    tStart = chunks(iChunk, 1);
    tEnd   = chunks(iChunk, 2);

    % Create Firing Rate Matrix
    tVec = tStart : binSize : tEnd;
    nBins = length(tVec) - 1;
    frMat = nan(nUnits, nBins);

    for iUnit = 1:nUnits
        st = spktimes{iUnit};
        st = st(st >= tStart & st <= tEnd);
        frMat(iUnit, :) = histcounts(st, tVec, 'Normalization', 'countdensity');
    end

    % Calculate Drift
    chunkStat = drift_calc(frMat, 'flgPlot', flgPlot, ...
        'thrFr', params.thrFr, ...
        'thrLin', params.thrLin, ...
        'limUnit', params.limUnit);

    % Remove extra info field from chunkStat to avoid confusion during concatenation
    if isfield(chunkStat, 'info')
        chunkStat = rmfield(chunkStat, 'info');
    end

    statsArr{iChunk} = chunkStat;

end



%% ========================================================================
%  CONCATENATE OUTPUT
%  ========================================================================

% Concatenate using catfields with 'addim' to stack across chunks
% This will effectively make drate [nChunks x 1], etc.
statsArr = [statsArr{:}];   % Flatten from cell array of structs to struct array
catStats = catfields(statsArr, 'addim', false, [], true);

% Merge into drft
flds = fieldnames(catStats);
for iFld = 1:length(flds)
    drft.(flds{iFld}) = catStats.(flds{iFld});
end
drft.dt_corr = squeeze(drft.dt_corr);

%% ========================================================================
%  SAVE
%  ========================================================================

if flgSave
    [~, basename] = fileparts(basepath);
    fname = fullfile(basepath, [basename, '.drift.mat']);
    save(fname, 'drft');
end

end     % EOF


%% ========================================================================
%  NOTE: STATE-DEPENDENT DRIFT ANALYSIS
%  ========================================================================
%  While this function computes drift across continuous time chunks, analyzing
%  population stability often requires accounting for behavioral states
%  (e.g., REM, NREM, Wake). The relationship between drift and state can be
%  approached via several theoretical frameworks:
%
%  1. THE CONCATENATED STATE HYPOTHESIS
%     This approach assumes that representational drift is a function of the
%     "time spent in state" rather than absolute time. By concatenating all
%     epochs of a specific state (e.g., NREM) and ignoring the intervening
%     wakefulness, one tests if the manifold evolves strictly due to state-
%     specific dynamics.
%     - Implementation: Construct 'frMat' by extracting columns corresponding
%       to specific state bouts and joining them.
%     - Implication: If drift is linear in concatenated time, the state
%       itself drives the plasticity.
%
%  2. THE ELAPSED TIME HYPOTHESIS (MASKING)
%     This approach posits that drift is a continuous, background process
%     that occurs regardless of state, but we can only observe it during
%     specific stable epochs. Here, the absolute timestamps are preserved,
%     but firing rate vectors are only sampled during the state of interest.
%     - Implementation: Compute 'frMat' on the full timeline, then set columns
%       outside the desired state to NaN.
%     - Implication: Allows comparison of drift rates (slopes) between states
%       over the same absolute duration.
%
%  3. THE STATE-TRANSITION MODEL
%     This framework investigates whether the drift occurs *within* states or
%     at the *transitions* between them. By computing population vectors for
%     individual bouts and correlating them, one can assess if the
%     representation 'jumps' to a new configuration upon re-entry into a
%     state.
%     - Implementation: Treat each bout as a single data point (Population
%       Vector) rather than binning within it.
%
%  4. IMPACT OF SLEEP ARCHITECTURE
%     In datasets where sleep architecture varies (e.g., fragmented sleep),
%     the choice of method defines the result. Method 1 (Concatenation)
%     normalizes for total sleep amount, whereas Method 2 (Elapsed) normalizes
%     for total recording time.
%  ========================================================================
