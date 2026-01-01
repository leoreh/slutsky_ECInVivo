function C = spk2cyto(spktimes, varargin)
% SPK2CYTO Converts spike times to cytosolic calcium concentration.
%
%   C = SPK2CYTO(SPKTIMES, ...) convolves a spike train with an exponential
%   decay kernel to simulate cytosolic calcium dynamics. It supports
%   non-linear ISI gain weighting to boost calcium influx for high-frequency
%   bursts.
%
%   INPUTS:
%       spktimes    - (vec) Spike times (in seconds).
%       varargin    - (param/value) Optional parameters:
%                     't'       : (vec) Time vector for simulation (s).
%                                 If empty, constructed from dt and Duration.
%                     'dt'      : (num) Time step {0.001} (s).
%                     'tauC'    : (num) Cytosolic decay time constant {0.1} (s).
%                     'isiGain' : (num) Gain for non-linear ISI weighting {0}.
%                                 Weights spikes by 1 + Gain * ((10ms - ISI)/10ms)^2.
%                                 Only applies to ISIs < 10ms.
%
%   OUTPUTS:
%       C           - (vec) Cytosolic Calcium Concentration (same length as t).
%
%   See also: SPK2CA, SPK2CA_SIM, CYTO2MITO

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spktimes', @isnumeric);
addParameter(p, 't', [], @isnumeric);
addParameter(p, 'dt', 0.001, @isnumeric);
addParameter(p, 'tauC', 0.1, @isnumeric);
addParameter(p, 'isiGain', 0, @isnumeric);

parse(p, spktimes, varargin{:});
par = p.Results;

% Ensure row vector for standard processing
spktimes = spktimes(:)';

%% ========================================================================
%  TIME VECTOR
%  ========================================================================

if isempty(par.t)
    maxTime = max(spktimes) + (par.tauC * 5);
    t = 0:par.dt:maxTime;
else
    t = par.t;
end
nPoints = length(t);

%% ========================================================================
%  DISCRETIZE COMPUTE WEIGHTS
%  ========================================================================

% Map Spikes to Time Indices

% Convert time to index (1-based)
idx = round(spktimes / par.dt) + 1;

% Filter out-of-bounds spikes
valid = idx >= 1 & idx <= nPoints;
idx = idx(valid);
spktimes = spktimes(valid); % Align spikes with valid indices

if isempty(idx)
    C = zeros(1, nPoints);
    return;
end

% Calculate Weights (ISI Gain)
w = ones(size(idx));

if par.isiGain > 0 && length(spktimes) > 1
    % Calculate ISIs (Inter-Spike Intervals)
    % Prepend infinity for the first spike
    isis = [inf, diff(spktimes)];

    % Threshold for short ISIs (11ms to capture up to 10ms comfortably)
    thresh = 0.011;
    isShort = isis < thresh;

    % Non-linear Boost Formulation:
    % Weight = 1 + Gain * ((Thresh - ISI) / Thresh)^2
    % This creates a quadratic boost as ISI approaches 0.
    boost = par.isiGain * ((thresh - isis(isShort)) ./ thresh).^2;

    w(isShort) = 1 + boost;
end

% Create Weighted Spike Vector (S)
% -------------------------------------------------------------------------
S = zeros(1, nPoints, 'single');

% Accumulate weights into bins
% Note: Multiple spikes in same bin are summed
for iBin = 1:length(idx)
    S(idx(iBin)) = S(idx(iBin)) + w(iBin);
end

%% ========================================================================
%  COMPUTE CALCIUM (FILTER)
%  ========================================================================

% Recursive Filter Implementation (Leaky Integrator)
% C(t) = C(t-1)*alpha + S(t)
alphaC = exp(-par.dt / par.tauC);

% Filter function implements: a(1)*y(n) = b(1)*x(n) + a(2)*y(n-1) ...
% We want: y(n) = S(n) + alpha*y(n-1)
% So: y(n) - alpha*y(n-1) = S(n)
% Coefficients: a = [1, -alpha], b = 1
C = filter(1, [1, -alphaC], S);

end     % EOF
