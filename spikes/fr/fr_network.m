function frNet = fr_network(spktimes, varargin)
% FR_NETWORK Wrapper for network-based firing rate calculations.
%
%   frNet = FR_NETWORK(SPKTIMES, ...) calculates dimensionality and
%   pairwise correlations for the given spike times, potentially split into
%   chunks if a window size is specified.
%
%   INPUTS:
%       spktimes    - (cell) Spike times per unit (e.g., {unit1, unit2}).
%       varargin    - (param/value) Optional parameters:
%                     'winLim'     : (num) [Start, End] time limit for analysis.
%                                    Default: [0, max(spktimes)].
%                     'winSize'    : (num) Size of chunks to split winLim into.
%                                    If empty, uses reduced winLim as one chunk.
%                     'binSize'    : (num) Bin size for FR calc {0.1}
%                     'basepath'   : (char) Base path for saving {pwd}
%                     'flgSave'    : (log) Save result as frNet struct {false}
%
%   OUTPUTS:
%       frNet       - (struct) Network statistics structure:
%                     .time       : (nChunks x 2) Time windows used.
%                     .dim        : (nChunks x 1) Dimensionality.
%                     .mcc        : (nChunks x 1) Denoised Mean Correlation.
%                     .mccRaw     : (nChunks x 1) Raw Mean Correlation.
%                     .funcon     : (nChunks x nUnits) Mean Functional Connectivity.
%                     .corr       : (nChunks x 1) Struct array of detailed corr stats.
%                     .info       : Info and parameters.
%
%   See also: DIM_CALC, FR_CORR, N2CHUNKS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'winLim', [], @isnumeric);
addParameter(p, 'winSize', [], @isnumeric);
addParameter(p, 'binSize', 0.1, @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', false, @islogical);

parse(p, spktimes, varargin{:});
winLim    = p.Results.winLim;
winSize   = p.Results.winSize;
binSize   = p.Results.binSize;
basepath  = p.Results.basepath;
flgSave   = p.Results.flgSave;

%% ========================================================================
%  INITIALIZE
%  ========================================================================

% Hardcoded Parameters for DIM
params.dim.method  = 'pr';
params.dim.thrVal  = 0.8;
params.dim.flgFrac = false;

% Hardcoded Parameters for CORR
params.corr.flgPlot = true;
params.corr.nShuffles = 40;
params.corr.zMet = 'shuffle';


% Handle winLim
if isempty(winLim)
    % Determine max time from spktimes
    maxTime = max(cellfun(@(x) max([0; x(:)]), spktimes));
    winLim = [0, maxTime];
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

% Initialize Output
frNet.time    = chunks;
frNet.dim     = nan(nChunks, 1);
frNet.mcc     = nan(nChunks, 1);
frNet.mccRaw  = nan(nChunks, 1);
frNet.funcon  = nan(nChunks, nUnits);
frNet.corr    = struct([]); % Will be array of structs

% Info
frNet.info.input   = p.Results;
frNet.info.params  = params;
frNet.info.binSize = binSize;
frNet.info.nUnits  = nUnits;

if nUnits < 2
    warning('Fewer than 2 units provided. Returning empty struct.');
    return;
end

%% ========================================================================
%  COMPUTE LOOP
%  ========================================================================

for iChunk = 1:nChunks

    tStart = chunks(iChunk, 1);
    tEnd   = chunks(iChunk, 2);

    % Create Firing Rate Matrix
    t = tStart : binSize : tEnd;
    nBins = length(t) - 1;
    frMat = nan(nUnits, nBins);

    for iUnit = 1:nUnits
        st = spktimes{iUnit};
        st = st(st >= tStart & st <= tEnd);
        counts = histcounts(st, t);
        frMat(iUnit, :) = counts ./ binSize;
    end

    % Dimensionality
    frNet.dim(iChunk) = dim_calc(frMat, ...
        'method', params.dim.method, ...
        'thrVal', params.dim.thrVal, ...
        'flgFrac', params.dim.flgFrac);

    % Correlations
    resCorr = fr_corr(frMat, ...
        'nShuffles', params.corr.nShuffles, ...
        'flgPlot', params.corr.flgPlot, ...
        'zMet', params.corr.zMet);

    frNet.mcc(iChunk)    = resCorr.mcc;
    frNet.mccRaw(iChunk) = resCorr.mccRaw;
    frNet.funcon(iChunk, :) = resCorr.funcon;

    % Store full result in struct array (grow it)
    if iChunk == 1
        frNet.corr = resCorr;
    else
        frNet.corr(iChunk) = resCorr;
    end

end

% Save
if flgSave
    [~, basename] = fileparts(basepath);
    fname = fullfile(basepath, [basename, '.frNet.mat']);
    save(fname, 'frNet');
end

end     % EOF
