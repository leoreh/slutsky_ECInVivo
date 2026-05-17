%% spontCa_threshSweep.m  Amp threshold sensitivity sweep.
%
% Independently sweeps the cyto and mito amp thresholds; at each step
% refilters events, recomputes cell-level rate/amp aggregates, and fits
% the T1 genotype contrast inline (skipping the full spontCa_stats path
% for speed). Two filter modes are run side by side:
%
%   all       - threshold only (flgPair = 'all').
%   triggered - threshold + flgPair = 'paired' (keep only events whose
%               pairLag falls inside winPair, the paired flag built in
%               spontCa2_metrics).
%
% USE
%   Run as a script after cache/spontCa_tbl.mat is loaded and metrics
%   are computed (tblCell, tblEvent, fs must be in workspace). The LOAD
%   section below is commented out by default - rely on the workspace
%   that mcu_spontCa already set up.


%% ========================================================================
%  LOAD
%  ========================================================================

spDir = fileparts(which('spontCa_detect'));
load(fullfile(spDir, 'cache', 'spontCa_tbl.mat'), 'tblCell', 'tblEvent', 'fs');
[tblCell, tblEvent] = spontCa2_metrics(tblCell, tblEvent, fs, 'aggFcn', 'mean');


%% ========================================================================
%  SWEEP
%  ========================================================================

threshGrid = [0, logspace(-2, -0.3, 30)];
nTh        = numel(threshGrid);
pairModes  = {'all', 'paired'};
modeLabels = {'all', 'triggered'};
nModes     = numel(pairModes);

% Per-mode storage. Each entry is (nTh x 2) for [rate, amp].
[betaCyto, seCyto, pCyto] = deal(cell(nModes, 1));
[betaMito, seMito, pMito] = deal(cell(nModes, 1));
for iM = 1:nModes
    betaCyto{iM} = nan(nTh, 2);  seCyto{iM} = nan(nTh, 2);  pCyto{iM} = nan(nTh, 2);
    betaMito{iM} = nan(nTh, 2);  seMito{iM} = nan(nTh, 2);  pMito{iM} = nan(nTh, 2);
end

tblCell0  = tblCell;
tblEvent0 = tblEvent;

for iM = 1:nModes
    flgPair = pairModes{iM};
    fprintf('=== filter mode: %s   flgPair=%s ===\n', ...
        modeLabels{iM}, flgPair);

    for iT = 1:nTh
        th = threshGrid(iT);
        fprintf('  step %d / %d   threshold %.4f\n', iT, nTh, th);

        % --- Cyto threshold sweep (drop low-amp cyto, keep mito) ---
        [tblC, tblE] = spontCa_filter(tblCell0, tblEvent0, ...
            'minAmp', [th; -inf], 'minEvents', [1; 1], ...
            'flgPair', flgPair, 'verbose', false);
        if ~isempty(tblC)
            [tblC, ~] = spontCa2_metrics(tblC, tblE, fs, ...
                'mode', 'cellOnly', 'aggFcn', 'mean');
            [betaCyto{iM}(iT,1), seCyto{iM}(iT,1), pCyto{iM}(iT,1)] = fitT1(tblC, 'rate', 'Cyto');
            [betaCyto{iM}(iT,2), seCyto{iM}(iT,2), pCyto{iM}(iT,2)] = fitT1(tblC, 'amp',  'Cyto');
        end

        % --- Mito threshold sweep (drop low-amp mito, keep cyto) ---
        [tblC, tblE] = spontCa_filter(tblCell0, tblEvent0, ...
            'minAmp', [-inf; th], 'minEvents', [1; 1], ...
            'flgPair', flgPair, 'verbose', false);
        if ~isempty(tblC)
            [tblC, ~] = spontCa2_metrics(tblC, tblE, fs, ...
                'mode', 'cellOnly', 'aggFcn', 'mean');
            [betaMito{iM}(iT,1), seMito{iM}(iT,1), pMito{iM}(iT,1)] = fitT1(tblC, 'rate', 'Mito');
            [betaMito{iM}(iT,2), seMito{iM}(iT,2), pMito{iM}(iT,2)] = fitT1(tblC, 'amp',  'Mito');
        end
    end
end


%% ========================================================================
%  PLOT
%  ========================================================================

figure('Color', 'w', 'Position', [100 100 1100 750]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

plotPanel(threshGrid, betaCyto, seCyto, pCyto, 1, 'Cyto amp threshold', 'Cyto rate', true);
plotPanel(threshGrid, betaCyto, seCyto, pCyto, 2, 'Cyto amp threshold', 'Cyto amp',  false);
plotPanel(threshGrid, betaMito, seMito, pMito, 1, 'Mito amp threshold', 'Mito rate', false);
plotPanel(threshGrid, betaMito, seMito, pMito, 2, 'Mito amp threshold', 'Mito amp',  false);

sgtitle('Threshold sweep   solid: all events    dashed: triggered only', ...
    'FontSize', 11);

% Link x-axes, set tight xlim.
allAxs = findall(gcf, 'Type', 'axes');
linkaxes(allAxs, 'x');
xlim(allAxs(1), [min(threshGrid), max(threshGrid)]);


%% ========================================================================
%  HELPERS
%  ========================================================================

function [beta, se, p] = fitT1(tblC, quantity, compartment)
% Inline T1 OLS contrast on log-transformed response (mirrors
% spontCa_stats.fitGenoOLS for dist='Log-Normal' but skips the dispatch
% overhead).
beta = NaN;  se = NaN;  p = NaN;
sub = tblC(tblC.compartment == compartment, :);
sub = sub(sub.(quantity) > 1e-8 & isfinite(sub.(quantity)), :);
if height(sub) < 4, return; end
tblFit = table(log(sub.(quantity)), sub.genotype, ...
    'VariableNames', {'y', 'genotype'});
try
    mdl = fitlm(tblFit, 'y ~ genotype');
    c   = mdl.Coefficients;
    idx = find(contains(c.Properties.RowNames, 'genotype') & ...
               ~contains(c.Properties.RowNames, ':'), 1);
    if isempty(idx), return; end
    beta = c.Estimate(idx);
    se   = c.SE(idx);
    p    = c.pValue(idx);
catch
end
end


function plotPanel(x, betaCell, seCell, pCell, col, xlbl, ttl, withLegend)
% Overlay 'all' (solid) and 'triggered' (dashed) curves on dual y-axes.
nexttile;


yyaxis left;
hAll  = errorbar(x, betaCell{1}(:, col), seCell{1}(:, col), '-',  'Color', 'k', 'LineWidth', 1);
hold on;
hTrig = errorbar(x, betaCell{2}(:, col), seCell{2}(:, col), '--', 'Color', 'k', 'LineWidth', 1);
yline(0, ':', 'Color', [0.5 0.5 0.5]);
ylabel('\beta (KO - Ctrl)');

yyaxis right;
plot(x, pCell{1}(:, col), '-',  'Color', 'r', 'LineWidth', 1);
hold on;
plot(x, pCell{2}(:, col), '--', 'Color', 'r', 'LineWidth', 1);
yline(0.05, ':', 'Color', 'r');
set(gca, 'YScale', 'log');
ylim([1e-4 1]);
ylabel('p');

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';

xlabel(xlbl);
title(ttl);
grid on;

if withLegend
    yyaxis left;
    legend([hAll, hTrig], {'all', 'triggered'}, 'Location', 'best', 'Box', 'off');
end
end
