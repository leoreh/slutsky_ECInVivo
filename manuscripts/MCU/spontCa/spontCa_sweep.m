function out = spontCa_sweep(tblCell, tblEvent, fs, varargin)
% SPONTCA_SWEEP  Threshold sweep with within- and cross-compartment readouts.
%
% out = SPONTCA_SWEEP(tblCell, tblEvent, fs, ...) sweeps an amp threshold
% in three directions:
%   cyto     - cyto threshold only (mito events kept untouched).
%   mito     - mito threshold only (cyto events kept untouched).
%   coupled  - same threshold applied to both compartments.
% At each step it refilters via spontCa_filter, refreshes cell aggregates
% via spontCa2_metrics(mode='cellOnly'), and fits:
%
%   WITHIN-compartment   log(var) ~ genotype           per compartment.
%                        Reports beta = KO - Ctrl and SE.
%
%   CROSS-compartment    log(Mito) ~ log(Cyto) * genotype
%                        on the unstacked cell-level table. Reports
%                        per-genotype slope + SE, interaction beta + SE,
%                        and per-genotype R^2 (the primary fit-quality
%                        readout for the new sweep dimension).
%
% Note: scanning a sweep for the threshold that maximizes R^2 is a
% post-hoc selection. For a published number, fix the threshold a
% priori and report the sweep as exploratory.
%
% NAME-VALUE
%   'var'        - column in tblCell to analyze. Default 'amp'.
%                  Supported by current metrics: amp, dur, flux, fluxRate,
%                  ampRate, rate, load, pairFlux, pairAmp, fluxOther.
%   'threshGrid' - amp thresholds. Default [0, logspace(-2, -0.3, 30)].
%   'flgPair'    - 'all' | 'paired' | 'unpaired'. Default 'all'.
%                  Passed through to spontCa_filter.
%   'aggFcn'     - 'mean' | 'median'. Default 'mean'.
%                  Passed through to spontCa2_metrics.
%   'flgPlot'    - draw the 2x3 figure. Default true.
%
% OUTPUT  (struct)
%   .threshGrid                          nTh x 1
%   .cyto.within.beta, .se               nTh x 2 (col 1 Cyto, col 2 Mito)
%   .cyto.cross.slopeCtrl, .slopeCtrlSE  nTh x 1
%   .cyto.cross.slopeKO,   .slopeKOSE    nTh x 1
%   .cyto.cross.dSlope,    .dSlopeSE     KO - Ctrl interaction
%   .cyto.cross.R2Ctrl, .R2KO            per-genotype OLS R^2
%   .cyto.counts.nCells, .nCyto, .nMito  surviving cell pairs / events
%   .mito.<same shape>                   mito-threshold sweep mirror
%   .coupled.<same shape>                same-threshold-both sweep
%   .baseline.nCells, .nCyto, .nMito     unfiltered totals (for % scaling)
%   .opts                                echo of name-value inputs
%
% See also: SPONTCA_FILTER, SPONTCA2_METRICS.


%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tblCell',  @istable);
addRequired(p, 'tblEvent', @istable);
addRequired(p, 'fs',       @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'var',        'amp',                       @(x) ischar(x) || isstring(x));
addParameter(p, 'threshGrid', [0, logspace(-2, -0.3, 30)], @isnumeric);
addParameter(p, 'flgPair',    'all',                       @(x) any(strcmpi(x, {'all','paired','unpaired'})));
addParameter(p, 'aggFcn',     'mean',                      @(x) any(strcmpi(x, {'mean','median'})));
addParameter(p, 'flgPlot',    true,                        @islogical);
parse(p, tblCell, tblEvent, fs, varargin{:});

varName    = char(p.Results.var);
threshGrid = p.Results.threshGrid(:);
flgPair    = p.Results.flgPair;
aggFcn     = p.Results.aggFcn;
flgPlot    = p.Results.flgPlot;
nTh        = numel(threshGrid);

if ~ismember(varName, tblCell.Properties.VariableNames)
    error('spontCa_sweep:badVar', 'var ''%s'' not found in tblCell.', varName);
end


%% ========================================================================
%  ALLOCATE
%  ========================================================================

dims = {'cyto', 'mito', 'coupled'};
result = struct();
for k = 1:numel(dims)
    d = dims{k};
    result.(d).within.beta       = nan(nTh, 2);
    result.(d).within.se         = nan(nTh, 2);
    result.(d).cross.slopeCtrl   = nan(nTh, 1);
    result.(d).cross.slopeCtrlSE = nan(nTh, 1);
    result.(d).cross.slopeKO     = nan(nTh, 1);
    result.(d).cross.slopeKOSE   = nan(nTh, 1);
    result.(d).cross.dSlope      = nan(nTh, 1);
    result.(d).cross.dSlopeSE    = nan(nTh, 1);
    result.(d).cross.R2Ctrl      = nan(nTh, 1);
    result.(d).cross.R2KO        = nan(nTh, 1);
    result.(d).counts.nCells     = nan(nTh, 1);
    result.(d).counts.nCyto      = nan(nTh, 1);
    result.(d).counts.nMito      = nan(nTh, 1);
end


%% ========================================================================
%  SWEEP
%  ========================================================================

tblCell0  = tblCell;
tblEvent0 = tblEvent;

nCellsBase = height(tblCell0) / 2;
nCytoBase  = sum(tblEvent0.compartment == 'Cyto');
nMitoBase  = sum(tblEvent0.compartment == 'Mito');

fprintf('=== spontCa_sweep   var=%s   flgPair=%s   aggFcn=%s ===\n', ...
    varName, flgPair, aggFcn);

for iT = 1:nTh
    th = threshGrid(iT);
    fprintf('  step %d / %d   threshold %.4f\n', iT, nTh, th);

    for iDim = 1:numel(dims)
        switch iDim
            case 1, minAmp = [th;   -inf];   % cyto only
            case 2, minAmp = [-inf;  th ];   % mito only
            case 3, minAmp = [th;    th ];   % coupled
        end
        d = dims{iDim};

        [tblC, tblE, info] = spontCa_filter(tblCell0, tblEvent0, ...
            'minAmp',    minAmp,    'minEvents', [1; 1], ...
            'flgPair',   flgPair,   'verbose',   false);

        result.(d).counts.nCells(iT) = info.nCellsAfter;
        result.(d).counts.nCyto(iT)  = info.nEventsAfter.Cyto;
        result.(d).counts.nMito(iT)  = info.nEventsAfter.Mito;

        if isempty(tblC), continue; end

        [tblC, ~] = spontCa2_metrics(tblC, tblE, fs, ...
            'mode', 'cellOnly', 'aggFcn', aggFcn);

        [bC, seCw] = fitWithin(tblC, varName, 'Cyto');
        [bM, seMw] = fitWithin(tblC, varName, 'Mito');
        [slc, seSlC, slk, seSlK, dsl, dslSe, r2c, r2k] = fitCross(tblC, varName);

        result.(d).within.beta(iT, :)    = [bC, bM];
        result.(d).within.se(iT, :)      = [seCw, seMw];
        result.(d).cross.slopeCtrl(iT)   = slc;
        result.(d).cross.slopeCtrlSE(iT) = seSlC;
        result.(d).cross.slopeKO(iT)     = slk;
        result.(d).cross.slopeKOSE(iT)   = seSlK;
        result.(d).cross.dSlope(iT)      = dsl;
        result.(d).cross.dSlopeSE(iT)    = dslSe;
        result.(d).cross.R2Ctrl(iT)      = r2c;
        result.(d).cross.R2KO(iT)        = r2k;
    end
end


%% ========================================================================
%  PACK OUTPUT
%  ========================================================================

out.threshGrid = threshGrid;
out.cyto       = result.cyto;
out.mito       = result.mito;
out.coupled    = result.coupled;
out.baseline   = struct('nCells', nCellsBase, ...
                        'nCyto',  nCytoBase, ...
                        'nMito',  nMitoBase);
out.opts       = p.Results;

if flgPlot
    plotSweep(out, varName);
end

end     % SPONTCA_SWEEP


%% ========================================================================
%  HELPERS
%  ========================================================================

function [beta, se] = fitWithin(tblC, varName, comp)
% OLS log(var) ~ genotype on a single compartment subset.
% Returns KO - Ctrl beta and SE.
beta = NaN; se = NaN;
sub = tblC(tblC.compartment == comp, :);
sub = sub(sub.(varName) > 1e-8 & isfinite(sub.(varName)), :);
if height(sub) < 4, return; end
tblFit = table(log(sub.(varName)), sub.genotype, ...
    'VariableNames', {'y', 'genotype'});
try
    mdl = fitlm(tblFit, 'y ~ genotype');
    rn  = mdl.Coefficients.Properties.RowNames;
    idx = find(contains(rn, 'genotype') & ~contains(rn, ':'), 1);
    if ~isempty(idx)
        beta = mdl.Coefficients.Estimate(idx);
        se   = mdl.Coefficients.SE(idx);
    end
catch
end
end


function [slC, seC, slK, seK, dSl, dSlSe, r2C, r2K] = fitCross(tblC, varName)
% Cell-level cross-compartment fit. log(Mito) ~ log(Cyto) on the table
% unstacked by compartment. Per-genotype slope, SE, R^2 come from
% per-genotype single fits; interaction beta and SE come from the
% combined model.
slC = NaN; seC = NaN; slK = NaN; seK = NaN;
dSl = NaN; dSlSe = NaN; r2C = NaN; r2K = NaN;

cols = {'genotype', 'sbjID', 'compartment', varName};
if ~all(ismember(cols, tblC.Properties.VariableNames)), return; end

sub = tblC(:, cols);
sub = sub(sub.(varName) > 1e-8 & isfinite(sub.(varName)), :);
if isempty(sub), return; end

try
    tblWide = unstack(sub, varName, 'compartment');
catch
    return;
end
if ~all(ismember({'Cyto', 'Mito'}, tblWide.Properties.VariableNames)), return; end

tblWide = tblWide(isfinite(tblWide.Cyto) & isfinite(tblWide.Mito) & ...
                  tblWide.Cyto > 0 & tblWide.Mito > 0, :);
if height(tblWide) < 8, return; end

tblWide.logCyto = log(tblWide.Cyto);
tblWide.logMito = log(tblWide.Mito);

% Per-genotype slope, SE, R^2.
genoLvls = unique(string(tblWide.genotype));
for i = 1:numel(genoLvls)
    s = tblWide(string(tblWide.genotype) == genoLvls(i), :);
    if height(s) < 4, continue; end
    try
        m  = fitlm(s, 'logMito ~ logCyto');
        sl = m.Coefficients.Estimate(2);
        sE = m.Coefficients.SE(2);
        r2 = m.Rsquared.Ordinary;
    catch
        continue;
    end
    if strcmpi(genoLvls(i), 'Control')
        slC = sl; seC = sE; r2C = r2;
    else
        slK = sl; seK = sE; r2K = r2;
    end
end

% Combined fit: interaction beta + SE.
try
    mdl = fitlm(tblWide, 'logMito ~ logCyto * genotype');
    rn  = mdl.Coefficients.Properties.RowNames;
    iI  = find(contains(rn, 'logCyto') & contains(rn, 'genotype'), 1);
    if ~isempty(iI)
        dSl   = mdl.Coefficients.Estimate(iI);
        dSlSe = mdl.Coefficients.SE(iI);
    end
catch
end
end


function plotSweep(out, varName)
% 3 rows x 4 cols.
%   Rows: cyto threshold / mito threshold / coupled threshold.
%   Cols: within beta, cross slope, cross R^2, surviving fraction.
clrComp = [0.20 0.50 0.85;   ...   % Cyto compartment
           0.55 0.20 0.65];        % Mito compartment
clrGeno = [0.00 0.00 0.00;   ...   % Control
           0.85 0.33 0.10];        % MCU-KO
clrCell = [0.30 0.30 0.30];        % cell-pair retention line

x = out.threshGrid;

hFig = figure('Color', 'w', 'Position', [60 60 1700 950]);
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

dims  = {'cyto', 'mito', 'coupled'};
xLbls = {'Cyto amp threshold', 'Mito amp threshold', 'Coupled amp threshold'};

base = out.baseline;

for iRow = 1:3
    d = out.(dims{iRow});

    % --- Col 1: within-compartment beta on var --------------------------
    nexttile;
    h1 = plotBand(x, d.within.beta(:,1), d.within.se(:,1), clrComp(1,:));
    h2 = plotBand(x, d.within.beta(:,2), d.within.se(:,2), clrComp(2,:));
    yline(0, ':', 'Color', [0.5 0.5 0.5]);
    set(gca, 'XScale', 'log'); grid on; box off;
    xlabel(xLbls{iRow});
    ylabel(sprintf('\\beta_{KO-Ctrl}  log(%s)', varName));
    title('Within-compartment contrast');
    legObj = [h1, h2]; legObj = legObj(~cellfun(@isempty, {h1, h2}));
    if iRow == 1 && numel(legObj) == 2
        legend(legObj, {'Cyto', 'Mito'}, 'Location', 'best', 'Box', 'off');
    end

    % --- Col 2: cross-compartment slope per genotype --------------------
    nexttile;
    h1 = plotBand(x, d.cross.slopeCtrl, d.cross.slopeCtrlSE, clrGeno(1,:));
    h2 = plotBand(x, d.cross.slopeKO,   d.cross.slopeKOSE,   clrGeno(2,:));
    yline(0, ':', 'Color', [0.5 0.5 0.5]);
    yline(1, ':', 'Color', [0.5 0.5 0.5]);
    set(gca, 'XScale', 'log'); grid on; box off;
    xlabel(xLbls{iRow});
    ylabel('slope  log(Mito) ~ log(Cyto)');
    title('Cross-compartment slope');
    legObj = [h1, h2]; legObj = legObj(~cellfun(@isempty, {h1, h2}));
    if iRow == 1 && numel(legObj) == 2
        legend(legObj, {'Control', 'MCU-KO'}, 'Location', 'best', 'Box', 'off');
    end

    % --- Col 3: cross-compartment R^2 per genotype ----------------------
    nexttile;
    h1 = plot(x, d.cross.R2Ctrl, '-', 'Color', clrGeno(1,:), 'LineWidth', 1.5);
    hold on;
    h2 = plot(x, d.cross.R2KO,   '-', 'Color', clrGeno(2,:), 'LineWidth', 1.5);
    [vC, iC] = max(d.cross.R2Ctrl);
    [vK, iK] = max(d.cross.R2KO);
    if isfinite(vC)
        plot(x(iC), vC, 'o', 'MarkerFaceColor', clrGeno(1,:), ...
             'MarkerEdgeColor', 'none', 'MarkerSize', 7, ...
             'HandleVisibility', 'off');
    end
    if isfinite(vK)
        plot(x(iK), vK, 'o', 'MarkerFaceColor', clrGeno(2,:), ...
             'MarkerEdgeColor', 'none', 'MarkerSize', 7, ...
             'HandleVisibility', 'off');
    end
    set(gca, 'XScale', 'log'); grid on; box off;
    xlabel(xLbls{iRow});
    ylabel('R^2');
    title('Cross-compartment fit quality');
    ylim([0 1]);
    if iRow == 1
        legend([h1, h2], {'Control', 'MCU-KO'}, 'Location', 'best', 'Box', 'off');
    end

    % --- Col 4: surviving fraction (cells + events) ---------------------
    nexttile;
    pctCells = 100 * d.counts.nCells / max(1, base.nCells);
    pctCyto  = 100 * d.counts.nCyto  / max(1, base.nCyto);
    pctMito  = 100 * d.counts.nMito  / max(1, base.nMito);
    h1 = plot(x, pctCells, '-', 'Color', clrCell,    'LineWidth', 1.8);  hold on;
    h2 = plot(x, pctCyto,  '-', 'Color', clrComp(1,:), 'LineWidth', 1.3);
    h3 = plot(x, pctMito,  '-', 'Color', clrComp(2,:), 'LineWidth', 1.3);
    yline(50, ':', 'Color', [0.5 0.5 0.5]);
    set(gca, 'XScale', 'log'); grid on; box off;
    xlabel(xLbls{iRow});
    ylabel('% remaining');
    title('Surviving fraction');
    ylim([0 110]);
    if iRow == 1
        legend([h1, h2, h3], ...
            {sprintf('cell pairs (N_0=%d)', base.nCells), ...
             sprintf('cyto events (N_0=%d)', base.nCyto), ...
             sprintf('mito events (N_0=%d)', base.nMito)}, ...
             'Location', 'best', 'Box', 'off');
    end
end

allAxs = findall(hFig, 'Type', 'axes');
linkaxes(allAxs, 'x');

sgtitle(sprintf('spontCa\\_sweep   var = %s   flgPair = %s   aggFcn = %s', ...
    varName, out.opts.flgPair, out.opts.aggFcn), 'FontSize', 11);
end


function h = plotBand(x, y, se, clr)
% Line with shaded +/- 1 SE band, returns line handle for legends.
h = [];
x = x(:); y = y(:); se = se(:);
ok = isfinite(y);
if ~any(ok), return; end
xV = x(ok); yV = y(ok); seV = se(ok);
if any(isfinite(seV) & seV > 0)
    seV(~isfinite(seV)) = 0;
    fill([xV; flipud(xV)], [yV + seV; flipud(yV - seV)], clr, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;
end
h = plot(xV, yV, '-', 'Color', clr, 'LineWidth', 1.5);
hold on;
end
