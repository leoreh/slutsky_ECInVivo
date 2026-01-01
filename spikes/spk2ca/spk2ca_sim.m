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
%                     'isiGain'       : (num) Gain for non-linear ISI weighting {0}
%                                       Weights spikes by 1 + Gain * ((10ms - ISI)/10ms)^2 for ISI < 10ms.
%                     'dt'            : (num) Simulation time step {0.001} (s)
%                     'dur'           : (num) Simulation duration {300} (s)
%                     'flgPlot'       : (log) Plot results {true}
%
%   OUTPUTS:
%       sim         - (struct) Simulation results structure.
%                     .tbl        : (table) Results table (Rate, TauC, Selectivity, Gain) plus traces.
%                     .res      : (struct) Summary matrices and optimal parameters.
%                     .params   : (struct) Input parameters.
%
%   See also: SPK2CA, PLOT_RASTER

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
addParameter(p, 'isiGain', 0, @isnumeric);
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
cSel   = cell(nTotal, 1);     % Vector (Sweep over Kd)
cGain  = cell(nTotal, 1);     % Vector (Sweep over Kd)
vOptSel  = zeros(nTotal, 1);  % Scalar (Optimal Kd)
vOptGain = zeros(nTotal, 1);  % Scalar (Optimal Kd)

% Trace Storage
cStT = cell(nTotal, 1);
cStP = cell(nTotal, 1);
cCT  = cell(nTotal, 1);
cCP  = cell(nTotal, 1);
cMT  = cell(nTotal, 1);
cMP  = cell(nTotal, 1);

cnt = 0;

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

fprintf('Starting Sweep [n=%d, %d Taus, %d Rates, %d Kds, Gain=%.2f]...\n', ...
    par.n, nTau, nRate, nKd, par.isiGain);

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

        % 2. Compute Cytosolic Calcium (spk2cyto)
        % This replaces spk2vec + filter and includes ISI gain logic
        cT = spk2cyto(stT, 't', t, 'dt', par.dt, 'isiGain', par.isiGain, 'tauC', currTau);
        cP = spk2cyto(stP, 't', t, 'dt', par.dt, 'isiGain', par.isiGain, 'tauC', currTau);

        % 4. Kd Sweep
        vSel  = zeros(1, nKd);
        vGain = zeros(1, nKd);

        for iK = 1:nKd
            valK = par.sweep_Kd(iK);
            Kn = valK ^ par.n;

            % Matrix Accumulation (Integration)
            mT = cyto2mito(cT, 'Kd', valK, 'n', par.n, 'dt', par.dt, 'flgVmax', true);
            mP = cyto2mito(cP, 'Kd', valK, 'n', par.n, 'dt', par.dt, 'flgVmax', true);

            ssT = mean(mT(idxSS));
            ssP = mean(mP(idxSS));

            % Metrics
            vGain(iK) = ssP / (ssT + eps);
            vSel(iK)  = (ssP - ssT) / (ssP + ssT + eps);
        end

        % Handle NaNs
        vSel(isnan(vSel)) = 0;

        % 5. Store Metrics
        strRate = sprintf('%.3g Hz', currRate);
        strTau  = sprintf('%.3g s', currTau);

        cRate{cnt} = strRate;
        cTau{cnt}  = strTau;
        cSel{cnt}  = vSel;
        cGain{cnt} = vGain;

        % 6. Pick Best Kd and Recalculate Final Traces
        peakVal = max(vSel);
        thresh = 0.95 * peakVal;
        idxOpt = find(vSel >= thresh, 1, 'first');
        if isempty(idxOpt), [~, idxOpt] = max(vSel); end

        bestKd = par.sweep_Kd(idxOpt);
        Kn = bestKd ^ par.n;

        % Recalculate Mito Flux for Best Kd (Using Helper)
        mT = cyto2mito(cT, 'Kd', bestKd, 'n', par.n, 'dt', par.dt, 'flgVmax', true);
        mP = cyto2mito(cP, 'Kd', bestKd, 'n', par.n, 'dt', par.dt, 'flgVmax', true);

        % Store Traces
        cStT{cnt} = stT;
        cStP{cnt} = stP;
        cCT{cnt}  = cT;
        cCP{cnt}  = cP;
        cMT{cnt}  = mT; %#ok<*NASGU>
        cMP{cnt}  = mP;

        % Store Scalar Optimal Metrics
        vOptSel(cnt)  = vSel(idxOpt);
        vOptGain(cnt) = vGain(idxOpt);

        % Store Best Params
        best.Kd   = bestKd;
        best.sel  = vSel(idxOpt);
        best.gain = vGain(idxOpt);
        best.rate = currRate;
        best.tau  = currTau;

        res.bestParams{iT, iR} = best;
        res.peakSel(iT, iR)    = peakVal;
        res.peakGain(iT, iR)   = best.gain;
    end
end

% Construct Simulation Table
% Add OptSelectivity and OptGain as scalar columns
tbl = table(cRate, cTau, cSel, cGain, vOptSel, vOptGain, cStT, cStP, cCT, cCP, cMT, cMP, ...
    'VariableNames', {'Rate', 'TauC', 'Selectivity', 'Gain', ...
    'OptSelectivity', 'OptGain', ...
    'stT', 'stP', 'cT', 'cP', 'mT', 'mP'});
tbl.Rate = categorical(tbl.Rate);
tbl.TauC = categorical(tbl.TauC);

% Reorder Categories properly
[~, idx] = natsort(unique(cRate));
tbl.Rate = reordercats(tbl.Rate, unique(cRate(idx))); %#ok
tbl.TauC = reordercats(tbl.TauC, unique(cTau));

sim.tbl    = tbl;
sim.res    = res;
sim.params = par;

fprintf('Done.\n');

%% ========================================================================
%  VISUALIZATION
%  ========================================================================

if par.flgPlot

    % --- Tranfer Function ---
    plot_transferFn(sim, t);

    % --- Parameter Sensitivity ---
    if nKd > 1
        tblGUI_xy(sim.params.sweep_Kd, sim.tbl, ...
            'yVar', 'Selectivity', ...
            'grpVar', 'Rate', ...
            'tileVar', 'none', ...
            'xLbl', 'Kd');
    end

    % --- Summary Analysis ---
    if nTau > 1
        [f2, ~] = plot_axSize('axShape', 'wide', 'axWidth', 800, 'hFig', []);
        set(f2, 'Name', 'TauC Analysis');

        tl2 = tiledlayout(f2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Selectivity vs Tau (Lines=Rate)
        axA = nexttile(tl2); hold(axA, 'on');
        for iR = 1:nRate
            plot(axA, par.tauC, res.peakSel(:, iR), 'o-', 'LineWidth', 1.5);
        end
        title(axA, 'Peak Selectivity vs TauC');
        xlabel(axA, 'TauC (s)'); ylabel(axA, 'Max Selectivity');
        grid(axA, 'on');

        % Panel B: Gain vs Tau (Lines=Rate)
        axB = nexttile(tl2); hold(axB, 'on');
        for iR = 1:nRate
            plot(axB, par.tauC, res.peakGain(:, iR), 's-', 'LineWidth', 1.5);
        end
        title(axB, 'Gain at Max Selectivity');
        xlabel(axB, 'TauC (s)'); ylabel(axB, 'Gain (P/T)');
        grid(axB, 'on');

        l2 = legend(axB, arrayfun(@(x) sprintf('R=%.3g', x), par.sweep_rate, 'UniformOutput', false));
        l2.Location = 'bestoutside';
    end

end

end     % EOF


%% ========================================================================
%  HELPER: GEN_SPKTIMES
%  ========================================================================
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






%% ========================================================================
%  HELPER: PLOT_TRANSFERFN
%  ========================================================================
function plot_transferFn(sim, t)
tbl = sim.tbl;

% Create Figure
hFig = figure('Name', 'Dynamics Explorer', 'Units', 'pixels', 'Position', [100 100 1200 800], 'Color', 'w');

% --- Header Controls ---
pnlH = uipanel('Parent', hFig, 'BorderType', 'none', 'BackgroundColor', 'w', ...
    'Units', 'normalized', 'Position', [0, 0.94, 1, 0.06]);

% Manually Position Controls (No uiflowcontainer)
uicontrol('Parent', pnlH, 'Style', 'text', 'String', 'Rate:', ...
    'BackgroundColor', 'w', 'FontWeight', 'bold', ...
    'Units', 'normalized', 'Position', [0.01, 0.2, 0.05, 0.6], 'HorizontalAlignment', 'right');

ddRate = uicontrol('Parent', pnlH, 'Style', 'popupmenu', 'String', categories(tbl.Rate), ...
    'Units', 'normalized', 'Position', [0.07, 0.2, 0.10, 0.6], ...
    'Callback', @updatePanels);

uicontrol('Parent', pnlH, 'Style', 'text', 'String', 'Tau:', ...
    'BackgroundColor', 'w', 'FontWeight', 'bold', ...
    'Units', 'normalized', 'Position', [0.20, 0.2, 0.05, 0.6], 'HorizontalAlignment', 'right');

ddTau = uicontrol('Parent', pnlH, 'Style', 'popupmenu', 'String', categories(tbl.TauC), ...
    'Units', 'normalized', 'Position', [0.26, 0.2, 0.10, 0.6], ...
    'Callback', @updatePanels);

uicontrol('Parent', pnlH, 'Style', 'text', 'String', 'Signal:', ...
    'BackgroundColor', 'w', 'FontWeight', 'bold', ...
    'Units', 'normalized', 'Position', [0.40, 0.2, 0.05, 0.6], 'HorizontalAlignment', 'right');

ddSig = uicontrol('Parent', pnlH, 'Style', 'popupmenu', 'String', {'Mito', 'Cyto', 'Both'}, ...
    'Units', 'normalized', 'Position', [0.46, 0.2, 0.10, 0.6], ...
    'Callback', @updatePanels);

% --- Body Panels ---
pnlBody = uipanel('Parent', hFig, 'BorderType', 'none', 'BackgroundColor', 'w', ...
    'Units', 'normalized', 'Position', [0, 0, 1, 0.94]);

% Create Axes Containers (Panels)
pnlTL = uipanel('Parent', pnlBody, 'Title', 'Tonic Raster', 'BackgroundColor', 'w', ...
    'Units', 'normalized', 'Position', [0.0, 0.5, 0.5, 0.5]);
pnlTR = uipanel('Parent', pnlBody, 'Title', 'Phasic Raster', 'BackgroundColor', 'w', ...
    'Units', 'normalized', 'Position', [0.5, 0.5, 0.5, 0.5]);
pnlBL = uipanel('Parent', pnlBody, 'Title', 'Tonic Signals', 'BackgroundColor', 'w', ...
    'Units', 'normalized', 'Position', [0.0, 0.0, 0.5, 0.5]);
pnlBR = uipanel('Parent', pnlBody, 'Title', 'Phasic Signals', 'BackgroundColor', 'w', ...
    'Units', 'normalized', 'Position', [0.5, 0.0, 0.5, 0.5]);

% Create Axes
axTL = axes('Parent', pnlTL);
axTR = axes('Parent', pnlTR);
axBL = axes('Parent', pnlBL);
axBR = axes('Parent', pnlBR);

% Initialize
updatePanels();

    function updatePanels(~, ~)
        % Get Selection
        idxR = get(ddRate, 'Value');
        catsR = categories(tbl.Rate);
        selRate = catsR{idxR};

        idxT = get(ddTau, 'Value');
        catsT = categories(tbl.TauC);
        selTau = catsT{idxT};

        idxS = get(ddSig, 'Value');
        optsS = {'Mito', 'Cyto', 'Both'};
        selSig = optsS{idxS};

        % Filter Table (Expected 1 row)
        subT = tbl(tbl.Rate == selRate & tbl.TauC == selTau, :);

        if isempty(subT)
            warning('No data for Rate=%s, Tau=%s', selRate, selTau);
            return;
        end

        % --- Plots ---

        % 1. Tonic Raster (TL)
        cla(axTL);
        plot_raster({subT.stT{1}}, 'hAx', axTL, 'xLim', [t(1) t(end)]);
        title(axTL, 'Tonic Raster');
        xlabel(axTL, 'Time (s)');

        % 2. Phasic Raster (TR)
        cla(axTR);
        plot_raster({subT.stP{1}}, 'hAx', axTR, 'xLim', [t(1) t(end)], 'clr', 'r');
        title(axTR, 'Phasic Raster');
        xlabel(axTR, 'Time (s)');

        % 3. Tonic Signals (BL)
        cla(axBL);
        hold(axBL, 'on');

        if strcmp(selSig, 'Mito') || strcmp(selSig, 'Both')
            plot(axBL, t, subT.mT{1}, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Mito');
        end
        if strcmp(selSig, 'Cyto') || strcmp(selSig, 'Both')
            plot(axBL, t, subT.cT{1}, 'b-', 'LineWidth', 0.5, 'DisplayName', 'Cyto');
        end

        hold(axBL, 'off');
        title(axBL, sprintf('Tonic (Sel=%.2f)', subT.OptSelectivity)); % Uses Scalar
        xlabel(axBL, 'Time (s)');
        legend(axBL, 'Location', 'best');
        grid(axBL, 'on');

        % 4. Phasic Signals (BR)
        cla(axBR);
        hold(axBR, 'on');

        if strcmp(selSig, 'Mito') || strcmp(selSig, 'Both')
            plot(axBR, t, subT.mP{1}, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Mito');
        end
        if strcmp(selSig, 'Cyto') || strcmp(selSig, 'Both')
            plot(axBR, t, subT.cP{1}, 'b-', 'LineWidth', 0.5, 'DisplayName', 'Cyto');
        end

        hold(axBR, 'off');
        title(axBR, sprintf('Phasic (Gain=%.2f)', subT.OptGain)); % Uses Scalar
        xlabel(axBR, 'Time (s)');
        legend(axBR, 'Location', 'best');
        grid(axBR, 'on');

        % Link Axes
        linkaxes([axTL, axTR, axBL, axBR], 'x');
    end

end


