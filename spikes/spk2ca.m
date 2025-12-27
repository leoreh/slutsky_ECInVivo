function ca = spk2ca(spktimes, varargin)
% SPK2CCA Converts spike times to continuous mitochondrial calcium accumulation.
%
%   ca = SPK2CCA(SPKTIMES, ...) calculates the Virtual Mitochondrial Calcium
%   Accumulation (Psi_mCa) metric by convolving spike trains with a cytosolic
%   decay kernel and passing the result through a Hill equation.
%
%   INPUTS:
%       spktimes    - (cell) Spike times per unit (in seconds).
%       varargin    - (param/value) Optional parameters:
%                     'winCalc' : (num) Limit analysis window (s)
%                     'dt'      : (num) Time step for discretization {0.001} (s)
%                     'tauC'    : (num) Cytosolic decay time constant {0.150} (s)
%                     'tauM'    : (num) Matrix decay time constant {30} (s)
%                     'n'       : (num) Hill coefficient {2.7}
%                     'Kd'      : (num) Dissociation constant {3.0} (USE)
%                     'binSize' : (num) Output bin size {60} (s)
%                     'flgPlot' : (log) Plot results per unit {true}
%
%   OUTPUTS:
%       ca         - (struct) Mitochondiral Calcium structure.
%                     .time     : (1 x nBins) Time vector (bin centers).
%                     .cyto     : (nUnits x nBins) Cytosolic Calcium (C), averaged.
%                     .mito     : ((nUnits x nBins) Integrated Flux (Psi_mCa), summed * dt.
%                     .params   : (struct) Parameters used.
%
%   See also: BRST_DYNAMICS, N2CHUNKS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'winCalc', [], @isnumeric);
addParameter(p, 'binSize', 60, @isnumeric);
addParameter(p, 'dt', 0.001, @isnumeric);
addParameter(p, 'tauC', 0.05, @isnumeric);
addParameter(p, 'tauM', 20, @isnumeric);
addParameter(p, 'n', 10, @isnumeric);
addParameter(p, 'Kd', 20, @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'basepath', pwd, @ischar);

parse(p, spktimes, varargin{:});
winCalc = p.Results.winCalc;
binSize = p.Results.binSize;
flgPlot = p.Results.flgPlot;
flgSave = p.Results.flgSave;
basepath = p.Results.basepath;
params  = p.Results;

%% ========================================================================
%  INITIALIZE
%  ========================================================================

nUnits = length(spktimes);

% Limit analysis to window
if ~isempty(winCalc)
    spktimes = cellfun(@(x) x(x >= winCalc(1) & x <= winCalc(2)), ...
        spktimes, 'UniformOutput', false);
end

% Determine global time range
maxTime = max(cellfun(@(x) max([0; x(:)]), spktimes));
if ~isempty(winCalc)
    if winCalc(2) > maxTime
        winCalc(2) = maxTime;
    end
    maxTime = min(maxTime, winCalc(2));
end

% Add a small buffer to the end to allow decay
buffer = params.tauC * 5;
calcTime = maxTime + buffer;

% Time vector (High Resolution)
t = 0:params.dt:calcTime;

% Use n2chunks to split duration
chunks = n2chunks('n', maxTime, 'chunksize', binSize, 'lastChunk', 'exclude');
edges = [chunks(:,1)-1, chunks(:,2)];
edges = unique(edges(:))';
tBins = edges(1:end-1) + diff(edges)/2;
nBins = length(tBins);

% Initialize Matrices (nUnits x nBins)
ca.time   = tBins;
ca.cyto   = zeros(nUnits, nBins);
ca.mito   = zeros(nUnits, nBins);
ca.params = p.Results;


%% ========================================================================
%  COMPUTE LOOP
%  ========================================================================

cytoMat   = zeros(nUnits, nBins);
mitoMat   = zeros(nUnits, nBins);

if isempty(gcp('nocreate'))
    parpool('local', 6);
end
dq = parallel.pool.DataQueue;
afterEach(dq, @(msg) fprintf('%s\n', msg));

parfor iUnit = 1:nUnits

    send(dq, sprintf('Working on Unit %d', iUnit));
    st = spktimes{iUnit};

    if isempty(st)
        continue;
    end

    [c, m] = transfer_fn(st, single(t), edges, params, false, iUnit);

    cytoMat(iUnit, :)   = (c);
    mitoMat(iUnit, :)   = (m);

end

ca.cyto   = cytoMat;
ca.mito   = mitoMat;

%% ========================================================================
%  PLOT & SAVE
%  ========================================================================

if flgPlot
    tbl = struct2table(rmfield(ca, {'time', 'params'}));
    tblGUI_xy(ca.time, tbl, 'yVar', 'mito');
end

if flgSave
    [~, basename] = fileparts(basepath);
    if isempty(basename)
        basename = 'mCa_results';
    end
    saveFile = fullfile(basepath, [basename, '.ca.mat']);
    save(saveFile, 'ca');
end

end     % EOF


%% ========================================================================
%  HELPER: COMPUTE UNIT DYNAMICS
%  ========================================================================

function [cyto, mito] = transfer_fn(st, t, edges, params, flgPlot, uid)
% SPK2MCA Helper to calculate calcium dynamics for a single unit.

nTime = length(t);
dt = params.dt;
n = params.n;
Kd = single(params.Kd);
nBins = length(edges) - 1;

% Discretize
% Convert spike times to indices
stIdx = round(st / dt) + 1;
stIdx(stIdx > nTime) = [];
stIdx(stIdx < 1) = [];

if isempty(stIdx)
    cyto = zeros(1, nBins);
    mito = zeros(1, nBins);
    return;
end

S = zeros(1, nTime, 'single');
for iTime = 1:length(stIdx)
    S(stIdx(iTime)) = S(stIdx(iTime)) + 1;
end

% Cytosolic Calcium (C)
% Recursive Filter Implementation
alphaC = exp(-dt / params.tauC);
C = filter(1, [1 -alphaC], S);

% Mitochondrial Flux (J)
C_n = C .^ n;
Kd_n = Kd ^ n;
J = C_n ./ (Kd_n + C_n);

% % Mitochondrial Concentration (M)
% alphaM = exp(-dt / params.tauM);
% M = filter(1 - alphaM, [1 -alphaM], J);

% Binning
binIdx = discretize(t, edges);
valid = ~isnan(binIdx);

cyto = accumarray(binIdx(valid)', C(valid)', [nBins 1], @mean)';
mito = accumarray(binIdx(valid)', J(valid)', [nBins 1], @sum)' * dt;

% Plotting
if flgPlot
    plot_spk2ca(t, st, C, J, uid, Kd);
end

end

%% ========================================================================
%  HELPER: PLOTTING
%  ========================================================================

function plot_spk2ca(t, st, C, J, uid, Kd)
% PLOT_MCA_UNIT Plot high-res dynamics for a single unit.

figure('Name', sprintf('Unit %d Dynamics', uid), 'Color', 'w', 'NumberTitle', 'off');
hTile = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Raster
ax1 = nexttile(hTile);
plot_raster({st}, 'hAx', ax1, 'plotType', 'vertline');
title(ax1, sprintf('Unit %d Raster', uid));
ylabel(ax1, 'Spikes');
grid(ax1, 'on');

% Cyto
ax2 = nexttile(hTile);
plot(ax2, t, C, 'b-', 'LineWidth', 1.5);
yline(ax2, Kd, 'k--', 'Kd');
title(ax2, 'Cytosolic Calcium (C)');
ylabel(ax2, 'Conc. (USE)');
grid(ax2, 'on');

% Mito Flux
ax3 = nexttile(hTile);
plot(ax3, t, J, 'r-', 'LineWidth', 1.5);
title(ax3, 'Mitochondrial Influx (J)');
ylabel(ax3, 'Flux');
grid(ax3, 'on');

linkaxes([ax1, ax2, ax3], 'x');
xlim(ax1, [0, max(t)]);

end


%% ========================================================================
%  NOTE: CONCENTRATION (CYTO) VS. ACCUMULATION (MITO)
%  ========================================================================
%  Summarizing high-resolution calcium traces into 1-minute bins requires
%  different mathematical treatments for the cytosol and the mitochondria.
%  The choice between averaging and summation is dictated by whether the
%  biological compartment acts as a transient mediator or a long-term
%  homeostatic integrator.
%
%  CYTOSOLIC CALCIUM (CONCENTRATION-BASED)
%  The cytosol acts as a transient signaling mediator characterized by
%  rapid influx and nearly instantaneous efflux via PMCA and NCX pumps.
%  Because this compartment does not "store" signal for long periods,
%  the 1-minute bin represents the average prevailing concentration.
%  This value reflects the general "pressure" applied to the mitochondrial
%  gate over time.
%
%  * Measures average concentration (C_cyto).
%  * Uses arithmetic mean across the bin.
%  * Reflects the current state of the gate.
%  * Units: Unitary Spike Equivalents (USE).
%
%  MITOCHONDRIAL CALCIUM (ACCUMULATION-BASED)
%  The mitochondrion acts as a "storage tank" or a computational observer
%  that integrates activity over time to set homeostatic metabolic points.
%  Because the homeostatic error signal is driven by the total volume of
%  calcium entry, we calculate the integral (sum) of the flux. This
%  represents the cumulative "dose" or "mass" transferred to the matrix
%  within each minute.
%
%  * Measures total integrated flux (Psi_mCa).
%  * Uses summation (Area Under Curve) across the bin.
%  * Reflects the cumulative metabolic load.
%  * Units: Integrated Flux (Psi).
%
%  WHY TOTAL IS SAVED OVER RATE
%  While the "Mean Rate" (Average Flux) and "Total" (Integrated Flux)
%  share the same temporal shape in a fixed-bin recording, the Total
%  is the more fundamental unit for homeostasis. It allows for the
%  quantification of the "Cumulative Error Signal" by simply summing bins
%  across longer experimental phases, such as the period of silencing.
%
%  * Total is additive across minutes.
%  * Rate = Total / 60 seconds.
%  * Preserves total calcium "mass" data.
%
%  THE NON-LINEAR INTEGRATION REQUIREMENT
%  To maintain biological accuracy, integration must occur at the full
%  sampling rate (e.g., 1000 Hz) before downsampling to 1-minute bins.
%  Averaging spikes before processing would mask the high-frequency
%  bursts required to overcome the non-linear MCU threshold (Kd).
%  Processing at high-res and then binning ensures the "total work done"
%  by the mitochondria is accurately captured.
%
%  * Captures burst-dependent thresholding.
%  * Prevents "averaging away" signal density.
%  * Decouples spike count from flux intensity.
%  ========================================================================


%% ========================================================================
%  NOTE: FILTER (IIR) VS CONVOLUTION (FIR)
%  ========================================================================
%  The implementation of the 'leaky integrator' for both cytosolic and
%  mitochondrial compartments can be achieved via either standard
%  convolution (FIR) or recursive filtering (IIR). This script utilizes
%  the 'filter' approach to improve performance and biological accuracy.
%
%  * Complexity: Filter operates with O(N) complexity, requiring only
%    one multiplication and addition per time step, whereas convolution
%    scales with the length of the decay kernel (tau).
%  * Memory Footprint: Eliminates the need to pre-compute and store
%    long kernel vectors for the cytosolic (150ms) and mitochondrial
%    (30s) decays, which is critical for high-resolution datasets.
%  * State Continuity: Unlike convolution, which is 'blind' to preceding
%    activity, the recursive filter allows for the passing of initial
%    conditions. This ensures the 'metabolic memory' of the matrix is
%    preserved across continuous 1-minute bins.
%
%  MATHEMATICAL EQUIVALENCE AND CAUSALITY
%  For the pure exponential decays modeled here, the recursive difference
%  equation is the exact mathematical counterpart to the physical
%  processes of calcium buffering and NCLX-mediated efflux.
%
%  * Causality: The filter implementation is strictly causal, mimicking
%    the real-time response of the Mitochondrial Observer.
%  * Tail Clipping: Filter returns a vector of the same length as the
%    input. To capture the full matrix decay after the final spike of
%    a bin, the spike train is zero-padded by 3*tauM.
%  * Parameter Stability: Using 'filter' allows the return to robust
%    physiological parameters without the computational penalty of
%    extremely long FIR kernels.
%% ========================================================================

%% ========================================================================
%  NOTE: STATE VS. AVERAGE: METABOLIC MOMENTUM
%  ========================================================================
%  While mitochondrial flux (J) and matrix calcium (M) are related, they
%  represent fundamentally different biological dimensions. Flux is a
%  'transaction' - the rate of entry into the organelle - while Matrix
%  Calcium is a 'state' - the current physiological load.
%
%  * Flux (J): Quantifies 'Total Work' done.
%  * Matrix (M): Quantifies 'Current Load' state.
%  * Binning: Averaging masks temporal memory.
%  * End-State: Preserves history and momentum.
%
%  JUSTIFICATION FOR THE END-STATE REPORTING
%  Averaging Matrix Calcium over a 60-second bin creates a high
%  correlation with firing rate (FR) because the system reaches a
%  steady-state equilibrium. Reporting the value at the end of the bin
%  (t = 60s) is a superior 'Homeostatic Anchor' for the following reasons:
%
%  1. Historical Continuity: The end-state reflects the residual stress
%     carried into the next processing window.
%  2. Burst Sensitivity: It distinguishes a bin ending in a burst from
%     one ending in silence, even if both have the same total spikes.
%  3. Leaky Integration: It correctly preserves the 30-second decay (tauM)
%     as a bridge across continuous experimental phases.
%% ========================================================================
