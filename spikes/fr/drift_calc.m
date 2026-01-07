function drft = drift_calc(frMat, varargin)
% DRIFT_CALC Calculates population vector drift over time.
%
%   drft = DRIFT_CALC(FRMAT, ...) computes the drift of the population
%   vector (PV) across the columns of the input matrix.
%
%   INPUTS:
%       frMat       - (matrix) Firing rates or Counts [nUnits x nWin].
%                     Each column represents a time window.
%       varargin    - (param/value) Optional parameters:
%                     'thrLin'   : (num) Min correlation pairs for linear fit.
%                     'thrFr'    : (num) Min FR to include unit in PV {0}.
%                     'limUnit'  : (num) Randomly subsample N units.
%                     'flgPlot'  : (log) Plot results {false}.
%
%   OUTPUTS:
%       drft        - (struct) Drift statistics:
%                     .dt_corr   : Correlation vs time lag (cell/matrix)
%                     .m_corr    : Mean correlation per lag
%                     .lin_coef  : Linear fit coefficients
%                     .drate     : Drift rate (slope)
%
%   See also: FR_NETWORK, DRIFT_PLOT

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'frMat', @isnumeric);
addParameter(p, 'thrLin', [], @isnumeric);
addParameter(p, 'thrFr', 0, @isnumeric);
addParameter(p, 'limUnit', [], @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);

parse(p, frMat, varargin{:});
thrLin   = p.Results.thrLin;
thrFr    = p.Results.thrFr;
limUnit  = p.Results.limUnit;
flgPlot  = p.Results.flgPlot;


%% ========================================================================
%  INITIALIZE
%  ========================================================================

[nUnits, nWin] = size(frMat);
nCorr = nWin - 1;

% Initialize Output
drft = struct('dt_corr', [], 'm_corr', [], 'lin_coef', [], 'drate', []);
drft.info = p.Results;
drft.info.nWin = nWin;
drft.info.nUnits = nUnits;

if nWin < 2
    warning('Not enough windows for drift calculation.');
    return;
end


%% ========================================================================
%  COMPUTE: POPULATION VECTORS
%  ========================================================================

% 1. Filter Units by FR Threshold
% -------------------------------
% Identify units with low average activity across windows (or max?)
% Usually we want units active within specific windows, but for global
% filtering:
meanUnitFr = mean(frMat, 2);
validUnits = meanUnitFr > thrFr;
pv = frMat(validUnits, :);

% 2. Subsample Units
% ------------------
currUnits = size(pv, 1);
if ~isempty(limUnit) && limUnit < currUnits
    unitIdx = randsample(currUnits, limUnit);
    pv = pv(unitIdx, :);
end

% 3. Normalize to Unit Vector (L2 Norm)
% -------------------------------------
% Each column becomes a unit vector
for iWin = 1:nWin
    v = pv(:, iWin);
    vNorm = norm(v);
    if vNorm > 0
        pv(:, iWin) = v / vNorm;
    else
        pv(:, iWin) = 0;
    end
end


%% ========================================================================
%  COMPUTE: CORRELATIONS & DRIFT
%  ========================================================================

% 1. Pairwise Correlations
% ------------------------
% pv is [Units x Win], corr(pv) gives [Win x Win]
pv_corr = corr(pv);

% 2. Organize by Time Lag (dt)
% ----------------------------
dt_corr = cell(nCorr, 1);
maxLag = nCorr;

for lag = 1:maxLag
    dt_corr{lag} = diag(pv_corr, lag);
end

% Pad for matrix form (creates [Win-1 x Lag] usually, or similar)
% Using a simple nan-pad implementation or cell approach
% Here we calculate mean immediately
m_corr = cellfun(@(x) mean(x, 'omitnan'), dt_corr);

% Store full distributions if needed
% (Skipping cell2padmat dependency if possible, or keeping simple)
drft.dt_corr = dt_corr;
drft.m_corr  = m_corr;

% 3. Linear Fit (Drift Rate)
% --------------------------
xAxis = (1:nCorr)';
yData = m_corr(:);

% Filter constraints for fit
isValid = ~isnan(yData);

if ~isempty(thrLin)
    % Count valid pairs per lag
    nPairs = cellfun(@(x) sum(~isnan(x)), dt_corr);
    isValid = isValid & (nPairs >= thrLin);
end

if sum(isValid) < 2
    drft.drate = NaN;
    drft.lin_coef = [NaN, NaN];
else
    % Polyfit
    lin_coef = polyfit(xAxis(isValid), yData(isValid), 1);
    drft.lin_coef = lin_coef;
    drft.drate = lin_coef(1);
end


%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot
    drift_plot(drft);
end

end     % EOF
