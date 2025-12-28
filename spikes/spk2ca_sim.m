function sim = spk2ca_sim(varargin)
% SPK2CA_SIM Simulation framework for mitochondrial calcium dynamics.
%
%   sim = SPK2CA_SIM(...) performs a parametric sweep of the mitochondrial
%   calcium transfer function to evaluate the sensitivity of matrix accumulation
%   to 'Tonic' vs. 'Phasic' firing patterns.
%
%   PLOTS:
%   1. Parameter Sensitivity: Selectivity vs Kd (Tiled by TauC).
%   2. Summary (if vector TauC): Selectivity/Gain vs TauC (Lines = Rate).
%   3. Dynamics Explorer: Detailed Traces at Optimal Parameters.
%
%   INPUTS:
%       varargin    - (param/value) Optional parameters:
%                     'sweep_rate'    : (vec) Mean firing rates (Hz)
%                     'sweep_Kd'      : (vec) vector of Kd values
%                     'tauC'          : (num) Cytosolic Decay (Scalar or Vec)
%                     'n'             : (num) Fixed Hill coefficient {4}
%                     'burstFreq'     : (num) Intra-burst frequency (Hz)
%                     'spikesPerBurst': (num) Number of spikes per burst
%
%   EXAMPLE:
%       sim = spk2ca_sim('tauC', logspace(-2, 0, 20));
%       sim = spk2ca_sim('tauC', [0.1]);
%
%   See also: SPK2CA, PLOT_RASTER

%% Arguments

p = inputParser;
addParameter(p, 'sweep_rate', [0.1, 0.5, 1, 5, 10], @isnumeric);
addParameter(p, 'sweep_Kd', logspace(-1, log10(20), 50), @isnumeric);
addParameter(p, 'n', 4, @isnumeric);
addParameter(p, 'tauC', 0.1, @isnumeric);
addParameter(p, 'burstFreq', 100, @isnumeric);
addParameter(p, 'spikesPerBurst', 8, @isnumeric);
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

%% Initialize

res.selCurves = cell(nTau, 1);  % Each cell: [nKd x nRate]
res.peakSel   = zeros(nTau, nRate);
res.peakGain  = zeros(nTau, nRate);
res.bestParams = cell(nTau, nRate);

% Time Vector
t = 0:par.dt:par.dur;
alphaM = exp(-par.dt / par.tauM);
filtBM = (1 - alphaM);

idxSS = t > (par.dur - 30);
if ~any(idxSS), idxSS = t > par.dur/2; end

fprintf('Starting Sweep [n=%d, %d Taus, %d Rates, %d Kds]...\n', ...
    par.n, nTau, nRate, nKd);

%% Sweep Loop

for iT = 1:nTau
    curTau = par.tauC(iT);
    alphaC = exp(-par.dt / curTau);

    selMap = zeros(nKd, nRate, 'single');

    for iR = 1:nRate
        curRate = par.sweep_rate(iR);

        % Generate Trains
        if par.spikesPerBurst * curRate > par.burstFreq
            fprintf('Skipping R=%.1f: Duty Cycle > 1\n', curRate);
            continue;
        end
        [stT, stP] = gen_trains(curRate, par.spikesPerBurst, par.burstFreq, par.dur);

        % Convert to Vector
        sT = spk2vec(stT, t, par.dt);
        sP = spk2vec(stP, t, par.dt);

        % Cytosolic Ca
        cT = filter(1, [1, -alphaC], sT);
        cP = filter(1, [1, -alphaC], sP);

        % Fixed Hill Powers
        cTn = cT .^ par.n;
        cPn = cP .^ par.n;

        % Kd Sweep
        vSel  = zeros(nKd, 1);
        vGain = zeros(nKd, 1);

        for iK = 1:nKd
            valK = par.sweep_Kd(iK);
            Kn = valK ^ par.n;

            % Flux
            jT = cTn ./ (cTn + Kn);
            jP = cPn ./ (cPn + Kn);

            % Matrix Ca
            mT = filter(filtBM, [1, -alphaM], jT);
            mP = filter(filtBM, [1, -alphaM], jP);

            ssT = mean(mT(idxSS));
            ssP = mean(mP(idxSS));

            % Metrics
            vGain(iK) = ssP / (ssT + eps);
            vSel(iK)  = (ssP - ssT) / (ssP + ssT + eps);
        end
        selMap(:, iR) = vSel;

        % Optimal Kd (Min val > 95% Peak)
        peakVal = max(vSel);
        thresh = 0.95 * peakVal;
        idxOpt = find(vSel >= thresh, 1, 'first');
        if isempty(idxOpt), [~, idxOpt] = max(vSel); end

        best.Kd = par.sweep_Kd(idxOpt);
        best.sel = vSel(idxOpt);
        best.gain = vGain(idxOpt);
        best.stT = stT;
        best.stP = stP;
        best.rate = curRate;
        best.tau = curTau;

        res.bestParams{iT, iR} = best;
        res.peakSel(iT, iR) = peakVal;
        res.peakGain(iT, iR) = best.gain;
    end

    res.selCurves{iT} = selMap;
end

sim.res = res;
sim.par = par;
fprintf('Done.\n');

%% Visualization

if par.flgPlot
    plot_results(sim, t, alphaM, filtBM);
end

end

%% Helper: Plotting
function plot_results(sim, t, alphaM, filtBM)
par = sim.par;
res = sim.res;
nTau = length(par.tauC);
nRate = length(par.sweep_rate);

clrs = parula(nRate);

% --- Fig 1: Parameter Sensitivity (Always On) ---
f1 = figure('Name', 'Sensitivity', 'Color', 'w', 'Position', [50 300 300*max(1,nTau) 400]);
tl1 = tiledlayout(1, nTau, 'TileSpacing', 'compact', 'Padding', 'compact');

for iT = 1:nTau
    ax = nexttile(tl1); hold on;
    selMap = res.selCurves{iT};

    for iR = 1:nRate
        if all(selMap(:, iR)==0), continue; end
        plot(ax, par.sweep_Kd, selMap(:, iR), '-', 'Color', clrs(iR,:), 'LineWidth', 1.5);
    end
    title(sprintf('TauC = %.3g s', par.tauC(iT)));
    xlabel('Kd'); ylabel('Selectivity');
    grid on;
    if iT == 1
        lgd = legend(arrayfun(@(x) sprintf('R=%.3g', x), par.sweep_rate, 'UniformOutput', false));
        lgd.Location = 'best';
    end
end

% --- Fig 2: Summary (If Multi-Tau) ---
if nTau > 1
    f2 = figure('Name', 'TauC Analysis', 'Color', 'w', 'Position', [50 50 800 400]);
    tl2 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Panel A: Selectivity vs Tau (Lines=Rate)
    axA = nexttile(tl2); hold on;
    for iR = 1:nRate
        plot(axA, par.tauC, res.peakSel(:, iR), 'o-', 'Color', clrs(iR,:), 'LineWidth', 1.5);
    end
    title('Peak Selectivity vs TauC');
    xlabel('TauC (s)'); ylabel('Max Selectivity');
    grid on;

    % Panel B: Gain vs Tau (Lines=Rate)
    axB = nexttile(tl2); hold on;
    for iR = 1:nRate
        plot(axB, par.tauC, res.peakGain(:, iR), 's-', 'Color', clrs(iR,:), 'LineWidth', 1.5);
    end
    title('Gain at Max Selectivity');
    xlabel('TauC (s)'); ylabel('Gain (P/T)');
    grid on;

    l2 = legend(axB, arrayfun(@(x) sprintf('R=%.3g', x), par.sweep_rate, 'UniformOutput', false));
    l2.Location = 'bestoutside';
end

% --- Fig 3: Dynamics Explorer ---
f3 = figure('Name', 'Dynamics', 'Color', 'w', 'Position', [660 50 800 600]);
tg = uitabgroup(f3);

for iT = 1:nTau
    for iR = 1:nRate
        bs = res.bestParams{iT, iR};
        if isempty(bs), continue; end

        titleStr = sprintf('R=%.1g, T=%.2g', bs.rate, bs.tau);
        tab = uitab(tg, 'Title', titleStr);
        tl = tiledlayout(tab, 2, 2, 'Padding', 'compact');

        % Re-Simulate Optimal Trace
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
        yyaxis(axT2, 'left'); plot(t, cT, 'b-'); ylabel('Cyto', 'Color', 'b');
        set(axT2, 'YColor', 'b');
        yyaxis(axT2, 'right'); plot(t, mT, 'r-'); ylabel('Mito (Log)', 'Color', 'r');
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


function [stT, stP] = gen_trains(meanRate, nSpk, freq, dur)
% Tonic
stT = 0:(1/meanRate):dur;

% Phasic
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
