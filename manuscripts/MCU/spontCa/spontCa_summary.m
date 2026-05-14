function [hFig, stats] = spontCa_summary(tbl, varargin)
% SPONTCA_SUMMARY One-figure + one-stats-table summary of the spontCa pipeline.
%
%   [HFIG, STATS] = SPONTCA_SUMMARY(TBL) auto-detects whether TBL is
%   tblEvent (has 'flux', 'pairFlux') or tblCell (has 'fluxRate', 'rate'),
%   builds a summary figure on log-scale axes, and runs the canonical
%   log-normal LME/OLS battery.
%
%   EVENT FIGURE (3x4): rows 1-2 = Cyto then Mito histograms of
%   amp/flux/pairFlux/pairAmp, row 3 = pairFlux vs flux and pairAmp vs amp
%   scatters, per compartment. pairAmp is the duration-free analog of
%   pairFlux (peak partner trace in the same window vs trace integral).
%
%   CELL FIGURE (3x4): same quantities as the event figure but on
%   per-cell means. Rows 1-2 = Cyto then Mito histograms of
%   amp/flux/pairFlux/pairAmp, row 3 = per-cell scatters of pairFlux vs
%   flux and pairAmp vs amp, per compartment. rate, fluxRate, and dur
%   stay in the stats table only.
%
%   tf, the cross-compartment paired scatter, and the cell-level
%   cyto-vs-mito scatter are omitted: tf algebraically restates
%   pairFlux/flux; the paired scatter is largely the same events as the
%   mito pairFlux~flux scatter; the cell-level cyto-vs-mito does not
%   condition on event-level temporal coupling.
%
%   OPTIONAL KEY-VALUE
%       'dist'       (char)    Response distribution               {'Log-Normal'}
%       'flgPlot'    (logical) Build the figure                    {true}
%       'flgStats'   (logical) Run the LME / OLS battery           {true}
%       'flgPrint'   (logical) Print pasteable Markdown summary    {true}
%       'flgLogAxis' (logical) Log-scale histograms and scatters   {true}
%       'cClrs'      (2x3)     Genotype colors. Default mcu_cfg.clr.grp.
%       'fSave'      (char)    PNG path                            {''}
%
%   See also: PLOT_HIST, PLOT_SCAT, LME_ANALYSE, MCU_CFG

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'dist',       'Log-Normal', @ischar);
addParameter(p, 'flgPlot',    true,  @islogical);
addParameter(p, 'flgStats',   true,  @islogical);
addParameter(p, 'flgPrint',   true,  @islogical);
addParameter(p, 'flgLogAxis', true,  @islogical);
addParameter(p, 'cClrs',      [],    @(x) isnumeric(x) && size(x,2)==3);
addParameter(p, 'fSave',      '',    @(x) ischar(x) || isstring(x));
parse(p, tbl, varargin{:});
args = p.Results;

cfg = mcu_cfg;
if isempty(args.cClrs), args.cClrs = cfg.clr.grp; end

varN = string(tbl.Properties.VariableNames);
if ismember("fluxRate", varN)
    mode = 'cell';                                  % tblCell has fluxRate; tblEvent does not
elseif all(ismember(["flux", "pairFlux", "pairIdx"], varN))
    mode = 'event';
else
    error('spontCa_summary:badInput', ['Input must be tblEvent ' ...
        '(flux, pairFlux, pairIdx) or tblCell (fluxRate).']);
end


%% ========================================================================
%  FIGURE
%  ========================================================================

hFig = [];
if args.flgPlot
    if strcmp(mode, 'event')
        hFig = plotEvent(tbl, args.cClrs, args.flgLogAxis);
    else
        hFig = plotCell(tbl, args.cClrs, args.flgLogAxis);
    end
    if ~isempty(args.fSave)
        exportgraphics(hFig, char(args.fSave), 'Resolution', 200);
    end
end


%% ========================================================================
%  STATS
%  ========================================================================

stats = struct();
stats.mode = mode;
stats.dist = args.dist;
stats.meta = describe(tbl, mode);

if args.flgStats
    if strcmp(mode, 'event')
        [stats.genotype, stats.scaling] = statsEvent(tbl, args.dist);
    else
        [stats.genotype, stats.scaling] = statsCell(tbl, args.dist);
    end
end

if args.flgPrint && args.flgStats
    printMarkdown(stats);
end

end     % SPONTCA_SUMMARY


%% ========================================================================
%  PLOTTING - EVENT
%  ========================================================================

function hFig = plotEvent(tbl, cClrs, flgLogAxis)

hFig = figure('Name', 'spontCa event-level summary', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.04 0.06 0.92 0.85]);
tl = tiledlayout(hFig, 3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, 'Event level summary', 'FontWeight', 'bold');

comps  = {'Cyto', 'Mito'};
quants = {'amp', 'flux', 'pairFlux', 'pairAmp'};
qLbl   = {'amp (dF/F)', 'flux (dF/F\cdots)', 'pairFlux (dF/F\cdots)', ...
          'pairAmp (dF/F)'};

scaleHist = 'linear'; if flgLogAxis, scaleHist = 'log'; end

% Rows 1-2: histograms.
for iC = 1:2
    tblC = tbl(tbl.compartment == comps{iC}, :);
    for iQ = 1:4
        ax = nexttile(tl, (iC-1)*4 + iQ);
        plot_hist(tblC, quants{iQ}, 'g', 'genotype', ...
            'hAx', ax, 'c', cClrs, 'flgKDE', true, 'flgStat', true, ...
            'scale', scaleHist);
        xlabel(ax, qLbl{iQ});
        ylabel(ax, 'pdf');
        title(ax, sprintf('%s | %s (n=%d)', quants{iQ}, comps{iC}, height(tblC)));
        if flgLogAxis, set(ax, 'XScale', 'log'); end
        if iC == 1 && iQ == 1
            legend(ax, 'Location', 'best', 'Box', 'off');
        end
    end
end

% Row 3: pairFlux ~ flux and pairAmp ~ amp scatters, per compartment.
epsThr = 1e-8;
scatOrder = { 'Cyto', 'flux', 'pairFlux'; ...
              'Mito', 'flux', 'pairFlux'; ...
              'Cyto', 'amp',  'pairAmp'; ...
              'Mito', 'amp',  'pairAmp' };
for iS = 1:4
    ax = nexttile(tl, 8 + iS);
    cmp = scatOrder{iS,1}; xVar = scatOrder{iS,2}; yVar = scatOrder{iS,3};
    tblC = tbl(tbl.compartment == cmp, :);
    if flgLogAxis
        set(ax, 'XScale', 'log', 'YScale', 'log');
        tblC = tblC(tblC.(xVar) > epsThr & tblC.(yVar) > epsThr, :);
    end
    plot_scat(tblC, xVar, yVar, 'g', 'genotype', ...
        'hAx', ax, 'c', cClrs, 'fitType', 'Linear', 'flgStats', true, ...
        'sz', 10, 'alpha', 0.35);
    xlabel(ax, xVar); ylabel(ax, yVar);
    title(ax, sprintf('%s vs %s | %s', yVar, xVar, cmp));
end

end


%% ========================================================================
%  PLOTTING - CELL
%  ========================================================================

function hFig = plotCell(tbl, cClrs, flgLogAxis)

hFig = figure('Name', 'spontCa cell-level summary', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.04 0.06 0.92 0.85]);
tl = tiledlayout(hFig, 3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, 'Cell level summary', 'FontWeight', 'bold');

comps  = {'Cyto', 'Mito'};
quants = {'amp', 'flux', 'pairFlux', 'pairAmp'};
qLbl   = {'amp (dF/F)', 'flux (dF/F\cdots)', ...
          'pairFlux (dF/F\cdots)', 'pairAmp (dF/F)'};

scaleHist = 'linear'; if flgLogAxis, scaleHist = 'log'; end

% Rows 1-2: histograms (same set as event level, on per-cell means).
for iC = 1:2
    tblC = tbl(tbl.compartment == comps{iC}, :);
    for iQ = 1:4
        if ~ismember(quants{iQ}, tbl.Properties.VariableNames), continue; end
        ax = nexttile(tl, (iC-1)*4 + iQ);
        plot_hist(tblC, quants{iQ}, 'g', 'genotype', ...
            'hAx', ax, 'c', cClrs, 'flgKDE', true, 'flgStat', true, ...
            'scale', scaleHist);
        xlabel(ax, qLbl{iQ});
        ylabel(ax, 'pdf');
        title(ax, sprintf('%s | %s (n=%d)', quants{iQ}, comps{iC}, height(tblC)));
        if flgLogAxis, set(ax, 'XScale', 'log'); end
        if iC == 1 && iQ == 1
            legend(ax, 'Location', 'best', 'Box', 'off');
        end
    end
end

% Row 3: per-cell coupling. Each point is one cell, axes are mean
% per-event values. Two scatter pairs: pairFlux vs flux, pairAmp vs amp.
hasCouple = all(ismember({'flux', 'pairFlux'}, tbl.Properties.VariableNames));
hasAmpCouple = all(ismember({'amp', 'pairAmp'}, tbl.Properties.VariableNames));
if hasCouple || hasAmpCouple
    epsThr = 1e-8;
    scatOrder = { 'Cyto', 'flux', 'pairFlux'; ...
                  'Mito', 'flux', 'pairFlux'; ...
                  'Cyto', 'amp',  'pairAmp'; ...
                  'Mito', 'amp',  'pairAmp' };
    for iS = 1:4
        cmp = scatOrder{iS,1}; xVar = scatOrder{iS,2}; yVar = scatOrder{iS,3};
        if ~all(ismember({xVar, yVar}, tbl.Properties.VariableNames)), continue; end
        ax = nexttile(tl, 8 + iS);
        tblC = tbl(tbl.compartment == cmp, :);
        if flgLogAxis
            set(ax, 'XScale', 'log', 'YScale', 'log');
            tblC = tblC(tblC.(xVar) > epsThr & tblC.(yVar) > epsThr, :);
        end
        plot_scat(tblC, xVar, yVar, 'g', 'genotype', ...
            'hAx', ax, 'c', cClrs, 'fitType', 'Linear', 'flgStats', true, ...
            'sz', 28, 'alpha', 0.65);
        xlabel(ax, sprintf('mean per-event %s', xVar));
        ylabel(ax, sprintf('mean per-event %s', yVar));
        title(ax, sprintf('%s vs %s (per cell) | %s', yVar, xVar, cmp));
    end
end

end


%% ========================================================================
%  STATS - EVENT
%  ========================================================================

function [tblGen, tblSc] = statsEvent(tbl, dist)

quants = {'amp', 'flux', 'pairFlux', 'pairAmp', 'dur'};
quants = quants(ismember(quants, tbl.Properties.VariableNames));
comps  = {'Cyto', 'Mito'};
useLme = true;

% Genotype contrast. Cyto dur is degenerate (= dt for every event), skip.
rows = {};
for iQ = 1:numel(quants)
    for iC = 1:numel(comps)
        if strcmp(quants{iQ}, 'dur') && strcmp(comps{iC}, 'Cyto'), continue; end
        subTbl = tbl(tbl.compartment == comps{iC}, :);
        frml = sprintf('%s ~ genotype + (1|sbjID)', quants{iQ});
        r = runFit(subTbl, frml, 'genotype', useLme, dist);
        rows(end+1, :) = {quants{iQ}, comps{iC}, r.beta, r.SE, r.t, ...
            r.df, r.p, r.transform};         %#ok<*AGROW>
    end
end
tblGen = cell2table(rows, 'VariableNames', ...
    {'quantity', 'compartment', 'beta', 'SE', 't', 'df', 'p', 'transform'});

% Within-compartment scaling: pairFlux ~ flux and pairAmp ~ amp.
scalingPairs = { 'pairFlux', 'flux'; 'pairAmp', 'amp' };
rows = {};
for iM = 1:size(scalingPairs, 1)
    yV = scalingPairs{iM, 1}; xV = scalingPairs{iM, 2};
    if ~all(ismember({yV, xV}, tbl.Properties.VariableNames)), continue; end
    for iC = 1:numel(comps)
        subTbl = tbl(tbl.compartment == comps{iC}, :);
        frml = sprintf('%s ~ %s * genotype + (1|sbjID)', yV, xV);
        rMain  = runFit(subTbl, frml, sprintf('main:%s', xV), useLme, dist);
        rInter = runFit(subTbl, frml, 'inter',                useLme, dist);
        rows(end+1, :) = {sprintf('%s ~ %s', yV, xV), comps{iC}, ...
            rMain.beta, rMain.p, rInter.beta, rInter.p};
    end
end
tblSc = cell2table(rows, 'VariableNames', ...
    {'model', 'compartment', 'beta_main', 'p_main', 'beta_inter', 'p_inter'});

end


%% ========================================================================
%  STATS - CELL
%  ========================================================================

function [tblGen, tblSc] = statsCell(tbl, dist)

quants = {'rate', 'amp', 'flux', 'fluxRate', 'pairFlux', 'pairAmp', 'dur'};
quants = quants(ismember(quants, tbl.Properties.VariableNames));
comps  = {'Cyto', 'Mito'};
useLme = false;

% Genotype contrast. Cyto dur is degenerate (mean over events that all
% have dur = dt), skip.
rows = {};
for iQ = 1:numel(quants)
    for iC = 1:numel(comps)
        if strcmp(quants{iQ}, 'dur') && strcmp(comps{iC}, 'Cyto'), continue; end
        subTbl = tbl(tbl.compartment == comps{iC}, :);
        frml = sprintf('%s ~ genotype', quants{iQ});
        r = runFit(subTbl, frml, 'genotype', useLme, dist);
        rows(end+1, :) = {quants{iQ}, comps{iC}, r.beta, r.SE, r.t, ...
            r.df, r.p, r.transform};
    end
end
tblGen = cell2table(rows, 'VariableNames', ...
    {'quantity', 'compartment', 'beta', 'SE', 't', 'df', 'p', 'transform'});

% Per-cell scaling: pairFlux ~ flux and pairAmp ~ amp.
scalingPairs = { 'pairFlux', 'flux'; 'pairAmp', 'amp' };
rows = {};
for iM = 1:size(scalingPairs, 1)
    yV = scalingPairs{iM, 1}; xV = scalingPairs{iM, 2};
    if ~all(ismember({yV, xV}, tbl.Properties.VariableNames)), continue; end
    for iC = 1:numel(comps)
        subTbl = tbl(tbl.compartment == comps{iC}, :);
        frml = sprintf('%s ~ %s * genotype', yV, xV);
        rMain  = runFit(subTbl, frml, sprintf('main:%s', xV), useLme, dist);
        rInter = runFit(subTbl, frml, 'inter',                useLme, dist);
        rows(end+1, :) = {sprintf('%s ~ %s', yV, xV), comps{iC}, ...
            rMain.beta, rMain.p, rInter.beta, rInter.p};
    end
end
if isempty(rows)
    tblSc = cell2table(cell(0, 6), 'VariableNames', ...
        {'model', 'compartment', 'beta_main', 'p_main', 'beta_inter', 'p_inter'});
else
    tblSc = cell2table(rows, 'VariableNames', ...
        {'model', 'compartment', 'beta_main', 'p_main', 'beta_inter', 'p_inter'});
end

end


%% ========================================================================
%  HELPERS - FIT + EXTRACT
%  ========================================================================

function r = runFit(tbl, frml, target, useLme, dist)
r = emptyResult();
isLog = strcmpi(dist, 'Log-Normal');

yVar = strtrim(extractBefore(frml, '~'));
yVar = char(yVar);

if isLog
    keep = tbl.(yVar) > 1e-8 & isfinite(tbl.(yVar));
    tbl  = tbl(keep, :);
end
if height(tbl) < 4, r.note = 'too few rows'; return; end

try
    if useLme
        [mdl, ~, ~] = lme_analyse(tbl, frml, ...
            'dist', dist, ...
            'flgPlot', false, 'verbose', false, 'flgStnd', false);
        r.transform = ternary(isLog, 'log', 'none');
    else
        tblFit = tbl;
        if isLog
            tblFit.(yVar) = log(tbl.(yVar));
            for vN = string(tbl.Properties.VariableNames)
                v = char(vN);
                if strcmp(v, yVar), continue; end
                col = tbl.(v);
                if ~isnumeric(col) || size(col, 2) ~= 1, continue; end
                if all(col > 0 & isfinite(col)) && any(col ~= 1)
                    tblFit.(v) = log(col);
                end
            end
            r.transform = 'log';
        else
            r.transform = 'none';
        end
        mdl = fitlm(tblFit, frml);
    end
    r = extractCoef(mdl, target, r);
catch ME
    r.note = ME.message;
end
end

function r = extractCoef(mdl, target, r)
coefs = mdl.Coefficients;
if istable(coefs)
    names = string(coefs.Properties.RowNames);
else
    names = string(coefs.Name);
end

if startsWith(target, 'main:')
    xName = extractAfter(target, 'main:');
    idx = find(names == xName, 1);
elseif strcmp(target, 'inter')
    isInter = contains(names, ':') & contains(names, 'genotype');
    idx = find(isInter, 1);
elseif strcmp(target, 'genotype')
    isGeno = contains(names, 'genotype') & ~contains(names, ':');
    idx = find(isGeno, 1);
else
    idx = find(names == string(target), 1);
end
if isempty(idx), return; end

if istable(coefs)
    r.beta = coefs.Estimate(idx);
    r.SE   = coefs.SE(idx);
    r.t    = coefs.tStat(idx);
    r.p    = coefs.pValue(idx);
    if ismember('DF', coefs.Properties.VariableNames)
        r.df = coefs.DF(idx);
    elseif isa(mdl, 'LinearModel')
        r.df = mdl.DFE;
    end
else
    r.beta = coefs.Estimate(idx);
    r.SE   = coefs.SE(idx);
    r.t    = coefs.tStat(idx);
    r.p    = coefs.pValue(idx);
    r.df   = coefs.DF(idx);
end
r.name = char(names(idx));
end

function r = emptyResult()
r = struct('beta', NaN, 'SE', NaN, 't', NaN, 'df', NaN, 'p', NaN, ...
    'transform', '', 'name', '', 'note', '');
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end


%% ========================================================================
%  HELPERS - DATA RESHAPING
%  ========================================================================

function meta = describe(tbl, mode)
meta = struct();
meta.mode = mode;
meta.nRows = height(tbl);
if ismember('compartment', tbl.Properties.VariableNames)
    for c = {'Cyto', 'Mito'}
        meta.byCompartment.(c{1}) = sum(tbl.compartment == c{1});
    end
end
if ismember('genotype', tbl.Properties.VariableNames)
    gLvls = categories(tbl.genotype);
    for iG = 1:numel(gLvls)
        g = gLvls{iG};
        meta.byGenotype.(matlab.lang.makeValidName(g)) = sum(tbl.genotype == g);
    end
end
if ismember('sbjID', tbl.Properties.VariableNames)
    meta.nCells = numel(unique(tbl.sbjID));
end
end


%% ========================================================================
%  HELPERS - PRINT MARKDOWN
%  ========================================================================

function printMarkdown(stats)
fprintf('\n### %s level summary\n\n', upper(stats.mode));

m = stats.meta;
if isfield(m, 'byCompartment') && isfield(m, 'byGenotype')
    cnts = m.byCompartment; gen = m.byGenotype;
    fprintf('rows: %d (Cyto %d, Mito %d).  ', m.nRows, cnts.Cyto, cnts.Mito);
    f = fieldnames(gen);
    parts = arrayfun(@(i) sprintf('%s %d', strrep(f{i},'_','-'), gen.(f{i})), ...
        1:numel(f), 'UniformOutput', false);
    fprintf('by genotype: %s.', strjoin(parts, ', '));
    if isfield(m, 'nCells'), fprintf('  cells: %d.', m.nCells); end
    fprintf('\n\n');
end

if isfield(stats, 'genotype')
    fprintf('Genotype contrast (KO vs Control, per compartment)\n\n');
    printGenoTable(stats.genotype);
    fprintf('\n');
end
if isfield(stats, 'scaling')
    fprintf('Within-compartment scaling (interaction with genotype)\n\n');
    printScalingTable(stats.scaling);
    fprintf('\n');
end
end

function printGenoTable(t)
fprintf('| quantity | compartment | beta(KO-Ctrl) | p | sig | transform |\n');
fprintf('|---|---|---:|---:|:---:|:---:|\n');
for iR = 1:height(t)
    fprintf('| %s | %s | %s | %s | %s | %s |\n', ...
        t.quantity{iR}, t.compartment{iR}, ...
        fmtBeta(t.beta(iR)), fmtP(t.p(iR)), sigStars(t.p(iR)), ...
        cellChar(t.transform(iR)));
end
end

function printScalingTable(t)
fprintf('| model | compartment | beta(main) | p(main) | beta(x:KO) | p(int) | sig(int) |\n');
fprintf('|---|---|---:|---:|---:|---:|:---:|\n');
for iR = 1:height(t)
    fprintf('| %s | %s | %s | %s | %s | %s | %s |\n', ...
        t.model{iR}, t.compartment{iR}, ...
        fmtBeta(t.beta_main(iR)), fmtP(t.p_main(iR)), ...
        fmtBeta(t.beta_inter(iR)), fmtP(t.p_inter(iR)), ...
        sigStars(t.p_inter(iR)));
end
end

function s = fmtBeta(b)
if isnan(b), s = '—'; return; end
if abs(b) < 1e-3 || abs(b) >= 1e3
    s = sprintf('%+.2e', b);
else
    s = sprintf('%+.3f', b);
end
end

function s = fmtP(p)
if isnan(p), s = '—';
elseif p < 1e-4, s = '<1e-4';
elseif p < 1e-3, s = sprintf('%.4f', p);
elseif p < 0.01, s = sprintf('%.3f', p);
else,            s = sprintf('%.2f',  p);
end
end

function s = sigStars(p)
if isnan(p),    s = '';
elseif p<0.001, s = '***';
elseif p<0.01,  s = '**';
elseif p<0.05,  s = '*';
elseif p<0.1,   s = '.';
else,           s = 'ns';
end
end

function s = cellChar(c)
if iscell(c), c = c{1}; end
if isstring(c), c = char(c); end
if isempty(c), s = 'none'; else, s = c; end
end
