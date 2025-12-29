function sim = spk2ca_sim(varargin)
% SPK2CA_SIM Parametric simulation of mitochondrial calcium selectivity.
%
%   sim = SPK2CA_SIM(...) performs a sweep of biophysical parameters to evaluate
%   the selectivity of mitochondrial calcium accumulation to 'Tonic' vs. 'Phasic'
%   firing patterns.
%
%   INPUTS:
%       varargin    - (param/value) Optional parameters:
%                     'sweep_rate'    : (vec) Mean firing rates (Hz) { [0.1 0.5 1 5 10] }
%                     'sweep_Kd'      : (vec) Dissociation constants { logspace(-1, log10(20), 50) }
%                     'tauC'          : (vec) Cytosolic decay time constants (s) {0.1}
%                     'burstFreq'     : (num) Intra-burst frequency (Hz) {100}
%                     'spikesPerBurst': (num) Number of spikes per burst {8}
%                     'n'             : (num) Hill coefficient {4}
%                     'tauM'          : (num) Matrix decay time constant {20} (s)
%                     'dt'            : (num) Simulation time step {0.001} (s)
%                     'dur'           : (num) Simulation duration {300} (s)
%                     'flgPlot'       : (log) Plot results {true}
%
%   OUTPUTS:
%       sim         - (struct) Simulation results structure.
%                     .T        : (table) Results table (Rate, TauC, Selectivity, Gain).
%                     .res      : (struct) Summary matrices and optimal parameters.
%                     .params   : (struct) Input parameters.
%
%   See also: SPK2CA, TBLGUI_XY, PLOT_AXSIZE

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'sweep_rate', [0.1, 0.5, 1, 5, 10], @isnumeric);
addParameter(p, 'sweep_Kd', logspace(-1, 1, 20), @isnumeric);
addParameter(p, 'tauC', logspace(-2, 0, 20), @isnumeric);
addParameter(p, 'burstFreq', 100, @isnumeric);
addParameter(p, 'spikesPerBurst', 8, @isnumeric);
addParameter(p, 'n', 4, @isnumeric);
addParameter(p, 'tauM', 20, @isnumeric);
addParameter(p, 'dt', 0.001, @isnumeric);
addParameter(p, 'dur', 300, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);

parse(p, varargin{:});
par = p.Results;

% Validate Dimensions
nRate = length(par.sweep_rate);
nKd   = length(par.sweep_Kd);
nTau  = length(par.tauC);

%% ========================================================================
%  INITIALIZE
%  ========================================================================

% Initialization for Table
nTotal = nRate * nTau;
cRate  = cell(nTotal, 1);
cTau   = cell(nTotal, 1);
cSel   = cell(nTotal, 1);
cGain  = cell(nTotal, 1);
cnt    = 0;

% Summary containers
res.peakSel    = zeros(nTau, nRate);
res.peakGain   = zeros(nTau, nRate);
res.bestParams = cell(nTau, nRate);

% Simulation Vectors
t = 0:par.dt:par.dur;
alphaM = exp(-par.dt / par.tauM);
filtBM = (1 - alphaM);

% Steady State Window
idxSS = t > (par.dur - 30);
if ~any(idxSS), idxSS = t > par.dur/2; end

fprintf('Starting Sweep [n=%d, %d Taus, %d Rates, %d Kds]...\n', ...
    par.n, nTau, nRate, nKd);

%% ========================================================================
%  COMPUTE LOOP
%  ========================================================================

for iT = 1:nTau
    currTau = par.tauC(iT);
    alphaC = exp(-par.dt / currTau);

    for iR = 1:nRate
        currRate = par.sweep_rate(iR);
        cnt = cnt + 1;

        % 1. Generate Spike Times
        if par.spikesPerBurst * currRate > par.burstFreq
            % Invalid (Duty Cycle > 1)
            cRate{cnt} = sprintf('%.3g Hz', currRate);
            cTau{cnt}  = sprintf('%.3g s', currTau);
            cSel{cnt}  = nan(1, nKd);
            cGain{cnt} = nan(1, nKd);
            continue;
        end

        [stT, stP] = gen_spktimes(currRate, par.spikesPerBurst, par.burstFreq, par.dur);

        % 2. Convert to Binary Vector
        sT = spk2vec(stT, t, par.dt);
        sP = spk2vec(stP, t, par.dt);

        % 3. Cytosolic Calcium (Linear Filter)
        cT = filter(1, [1, -alphaC], sT);
        cP = filter(1, [1, -alphaC], sP);

        % Pre-calculate Powers (Hill)
        cTn = cT .^ par.n;
        cPn = cP .^ par.n;

        % 4. Kd Sweep
        vSel  = zeros(1, nKd);
        vGain = zeros(1, nKd);

        for iK = 1:nKd
            valK = par.sweep_Kd(iK);
            Kn = valK ^ par.n;

            % Mitochondrial Flux (Hill Equation)
            jT = cTn ./ (cTn + Kn);
            jP = cPn ./ (cPn + Kn);

            % Matrix Accumulation (Integration)
            mT = filter(filtBM, [1, -alphaM], jT);
            mP = filter(filtBM, [1, -alphaM], jP);

            ssT = mean(mT(idxSS));
            ssP = mean(mP(idxSS));

            % Metrics
            vGain(iK) = ssP / (ssT + eps);
            vSel(iK)  = (ssP - ssT) / (ssP + ssT + eps);
        end

        % Handle NaNs
        vSel(isnan(vSel)) = 0;

        % 5. Store Results
        cRate{cnt} = sprintf('%.3g Hz', currRate);
        cTau{cnt}  = sprintf('%.3g s', currTau);
        cSel{cnt}  = vSel;
        cGain{cnt} = vGain;

        % Find Optimal Kd (95% Threshold)
        peakVal = max(vSel);
        thresh = 0.95 * peakVal;
        idxOpt = find(vSel >= thresh, 1, 'first');
        if isempty(idxOpt), [~, idxOpt] = max(vSel); end

        best.Kd   = par.sweep_Kd(idxOpt);
        best.sel  = vSel(idxOpt);
        best.gain = vGain(idxOpt);
        best.stT  = stT;
        best.stP  = stP;
        best.rate = currRate;
        best.tau  = currTau;

        res.bestParams{iT, iR} = best;
        res.peakSel(iT, iR)    = peakVal;
        res.peakGain(iT, iR)   = best.gain;
    end
end

% Construct Simulation Table
T = table(cRate, cTau, cSel, cGain, 'VariableNames', {'Rate', 'TauC', 'Selectivity', 'Gain'});
T.Rate = categorical(T.Rate);
T.TauC = categorical(T.TauC);
T.Rate = reordercats(T.Rate, natsort(unique(cRate)));
T.TauC = reordercats(T.TauC, unique(cTau));

sim.T      = T;
sim.res    = res;
sim.params = par;

fprintf('Done.\n');

%% ========================================================================
%  VISUALIZATION
%  ========================================================================

if par.flgPlot
    plot_results(sim, t, alphaM, filtBM);
end

end     % EOF

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function plot_results(sim, t, alphaM, filtBM)
par = sim.params;
res = sim.res;
T   = sim.T;
nTau  = length(par.tauC);
nRate = length(par.sweep_rate);
clrs  = parula(nRate);

% --- Fig 1: Parameter Sensitivity (Interactive GUI) ---
% Plots Selectivity vs Kd. Tiles by TauC, Groups by Rate.
tblGUI_xy(par.sweep_Kd, T, ...
    'yVar', 'Selectivity', ...
    'grpVar', 'Rate', ...
    'tileVar', 'none', ...
    'xLbl', 'Kd');

% --- Fig 2: Summary Analysis (Only if Multi-Tau) ---
if nTau > 1
    [f2, ~] = plot_axSize('axShape', 'wide', 'axWidth', 800, 'hFig', []);
    set(f2, 'Name', 'TauC Analysis');

    tl2 = tiledlayout(f2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Panel A: Selectivity vs Tau (Lines=Rate)
    axA = nexttile(tl2); hold(axA, 'on');
    for iR = 1:nRate
        plot(axA, par.tauC, res.peakSel(:, iR), 'o-', 'Color', clrs(iR,:), 'LineWidth', 1.5);
    end
    title(axA, 'Peak Selectivity vs TauC');
    xlabel(axA, 'TauC (s)'); ylabel(axA, 'Max Selectivity');
    grid(axA, 'on');

    % Panel B: Gain vs Tau (Lines=Rate)
    axB = nexttile(tl2); hold(axB, 'on');
    for iR = 1:nRate
        plot(axB, par.tauC, res.peakGain(:, iR), 's-', 'Color', clrs(iR,:), 'LineWidth', 1.5);
    end
    title(axB, 'Gain at Max Selectivity');
    xlabel(axB, 'TauC (s)'); ylabel(axB, 'Gain (P/T)');
    grid(axB, 'on');

    l2 = legend(axB, arrayfun(@(x) sprintf('R=%.3g', x), par.sweep_rate, 'UniformOutput', false));
    l2.Location = 'bestoutside';
end

% --- Fig 3: Dynamics Explorer ---
[f3, ~] = plot_axSize('axShape', 'wide', 'axWidth', 800, 'hFig', []);
set(f3, 'Name', 'Dynamics Explorer');
tg = uitabgroup(f3);

for iT = 1:nTau
    for iR = 1:nRate
        bs = res.bestParams{iT, iR};
        if isempty(bs), continue; end

        titleStr = sprintf('R=%.1g, T=%.2g', bs.rate, bs.tau);
        tab = uitab(tg, 'Title', titleStr);
        tl = tiledlayout(tab, 2, 2, 'Padding', 'compact');

        % Re-Simulate
        sT = spk2vec(bs.stT, t, par.dt);
        sP = spk2vec(bs.stP, t, par.dt);

        alphaC = exp(-par.dt / bs.tau);
        cT = filter(1, [1, -alphaC], sT);
        cP = filter(1, [1, -alphaC], sP);

        Kn = bs.Kd ^ par.n;
        jT = (cT.^par.n) ./ ((cT.^par.n) + Kn);
        jP = (cP.^par.n) ./ ((cP.^par.n) + Kn);

        mT = filter(filtBM, [1, -alphaM], jT);
        mP = filter(filtBM, [1, -alphaM], jP);

        % Plot Tonic
        axT1 = nexttile(tl, 1);
        plot_raster({bs.stT}, 'hAx', axT1, 'xLim', [t(1) t(end)]);
        title(axT1, 'Tonic');

        axT2 = nexttile(tl, 3);
        yyaxis(axT2, 'left'); plot(t, cT, 'b-'); ylabel('Cyto (Blue)');
        set(axT2, 'YColor', 'b');
        yyaxis(axT2, 'right'); plot(t, mT, 'r-'); ylabel('Mito (Log, Red)');
        set(axT2, 'YColor', 'r', 'YScale', 'log');
        title(axT2, sprintf('Tonic Response (Sel=%.2f)', bs.sel));
        grid(axT2, 'on');

        % Plot Phasic
        axP1 = nexttile(tl, 2);
        plot_raster({bs.stP}, 'hAx', axP1, 'xLim', [t(1) t(end)], 'clr', 'r');
        title(axP1, 'Phasic');

        axP2 = nexttile(tl, 4);
        yyaxis(axP2, 'left'); plot(t, cP, 'b-');
        set(axP2, 'YColor', 'b');
        yyaxis(axP2, 'right'); plot(t, mP, 'r-');
        set(axP2, 'YColor', 'r', 'YScale', 'log');
        title(axP2, sprintf('Phasic Response (Kd=%.1f)', bs.Kd));
        grid(axP2, 'on');

        linkaxes([axT1, axP1, axT2, axP2], 'x');
    end
end
end


function [stT, stP] = gen_spktimes(meanRate, nSpk, freq, dur)
% Tonic (Regular)
stT = 0:(1/meanRate):dur;

% Phasic (Bursts)
period = nSpk / meanRate;
bIsi   = 1/freq;
bDur   = (nSpk-1)*bIsi;
burst  = 0:bIsi:bDur;

nCycles = floor(dur/period);
stP = [];
for i = 0:(nCycles-1)
    stP = [stP, burst + (i*period)]; %#ok
end
end

function S = spk2vec(st, t, dt)
n = length(t);
idx = round(st/dt) + 1;
idx(idx>n | idx<1) = [];
S = zeros(1, n, 'single');
S(idx) = 1;
end
