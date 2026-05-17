function stats = spontCa_stats(tblCell, tblEvent, varargin)
% SPONTCA_STATS  Predefined statistical tables for the spontCa pipeline.
%
%   STATS = SPONTCA_STATS(TBLCELL, TBLEVENT, ...) runs four predefined
%   tests and packages the results as four tables:
%
%       T1 (cellGeno)  Cell-level genotype contrasts (KO vs Control)
%                      for each scalar metric per compartment. OLS, log
%                      response.
%
%       T2 (eventGeno) Event-level genotype contrasts for each event
%                      metric per compartment. LME with random cell
%                      intercept (1|sbjID), log response.
%
%       T3 (cellRel)   Cell-level cross-compartment relationships
%                      (compensation framework). Specific pairs:
%                          mito.amp   ~ cyto.amp   * genotype
%                          mito.flux  ~ cyto.flux  * genotype
%                          mito.load  ~ cyto.load  * genotype
%                      Reports main slope (Ctrl), genotype interaction,
%                      and joint F-test on all genotype coefficients
%                      (the compensation-model null).
%
%       T4 (eventRel)  Event-level within-compartment paired scaling.
%                      Specific pairs per compartment:
%                          pairFlux  ~ flux * genotype
%                          pairAmp   ~ amp  * genotype
%                          fluxOther ~ flux * genotype
%                      LME with random cell intercept. Reports main
%                      slope, interaction, and the joint F-test.
%
%                      Row semantics under the single-pass pairing:
%                          pairFlux/pairAmp - partner event's event-
%                              level flux/amp on the partner's own scale
%                              (cyto row: paired mito; mito row: max-amp
%                              cyto among claimers).
%                          fluxOther - partner-trace window integral,
%                              detection-free (cyto row: mito trace in
%                              [s - winPair(1), s + winPair(2)]; mito
%                              row: cyto trace in [m - winPair(2), m +
%                              winPair(1)], causally flipped).
%                      Cyto and Mito subset slopes are NOT numerically
%                      reciprocal - different samples, different aggreg-
%                      ation, different flux units across compartments.
%
%   OPTIONAL (Name-Value)
%       'flgPrint' (logical) Print Markdown summary. Default true.
%       'dist'     (char)    Response distribution for LME. Default
%                            'Log-Normal'. OLS rows always log-transform
%                            the response.
%
%   See also: SPONTCA2_METRICS, MCU_SPONTCA.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tblCell',  @istable);
addRequired(p, 'tblEvent', @istable);
addParameter(p, 'flgPrint', true,         @islogical);
addParameter(p, 'dist',     'Log-Normal', @ischar);
parse(p, tblCell, tblEvent, varargin{:});
args = p.Results;


%% ========================================================================
%  RUN ALL FOUR TABLES
%  ========================================================================

stats           = struct();
stats.cellGeno  = genoCell(tblCell,  args.dist);
stats.eventGeno = genoEvent(tblEvent, args.dist);
stats.cellRel   = relCell(tblCell);
stats.eventRel  = relEvent(tblEvent, args.dist);

if args.flgPrint
    printAll(stats);
end

end     % SPONTCA_STATS


%% ========================================================================
%  T1: cell-level genotype contrasts (OLS, log response)
%  ========================================================================

function tbl = genoCell(tblCell, dist)
quants = {'rate', 'amp', 'dur', 'flux', 'fluxRate', 'ampRate', ...
          'pairFlux', 'pairAmp', 'fluxOther', 'load'};
quants = quants(ismember(quants, tblCell.Properties.VariableNames));
comps  = {'Cyto', 'Mito'};
rows = cell(0, 6);
for iQ = 1:numel(quants)
    for iC = 1:numel(comps)
        if strcmp(quants{iQ}, 'dur') && strcmp(comps{iC}, 'Cyto'), continue; end
        sub  = tblCell(tblCell.compartment == comps{iC}, :);
        frml = sprintf('%s ~ genotype', quants{iQ});
        r    = fitGenoOLS(sub, frml, quants{iQ}, dist);
        rows(end+1, :) = {quants{iQ}, comps{iC}, r.beta, r.SE, r.p, r.transform}; %#ok<AGROW>
    end
end
tbl = cell2table(rows, 'VariableNames', ...
    {'quantity', 'compartment', 'beta_KOvsCtrl', 'SE', 'p', 'transform'});
end


%% ========================================================================
%  T2: event-level genotype contrasts (LME, log response)
%  ========================================================================

function tbl = genoEvent(tblEvent, dist)
quants = {'amp', 'dur', 'flux', 'pairFlux', 'pairAmp', 'fluxOther'};
quants = quants(ismember(quants, tblEvent.Properties.VariableNames));
comps  = {'Cyto', 'Mito'};
rows = cell(0, 6);
for iQ = 1:numel(quants)
    for iC = 1:numel(comps)
        if strcmp(quants{iQ}, 'dur') && strcmp(comps{iC}, 'Cyto'), continue; end
        sub  = tblEvent(tblEvent.compartment == comps{iC}, :);
        frml = sprintf('%s ~ genotype + (1|sbjID)', quants{iQ});
        r    = fitGenoLME(sub, frml, quants{iQ}, dist);
        rows(end+1, :) = {quants{iQ}, comps{iC}, r.beta, r.SE, r.p, r.transform}; %#ok<AGROW>
    end
end
tbl = cell2table(rows, 'VariableNames', ...
    {'quantity', 'compartment', 'beta_KOvsCtrl', 'SE', 'p', 'transform'});
end


%% ========================================================================
%  T3: cell-level cross-compartment relationships (OLS, log-log)
%  ========================================================================

function tbl = relCell(tblCell)
quants = {'amp', 'flux', 'load'};
quants = quants(ismember(quants, tblCell.Properties.VariableNames));

isC = tblCell.compartment == 'Cyto';
isM = tblCell.compartment == 'Mito';
wide = tblCell(isC, {'sbjID', 'genotype'});
for q = quants
    qC = ['cyto' upper(q{1}(1)) q{1}(2:end)];
    qM = ['mito' upper(q{1}(1)) q{1}(2:end)];
    wide.(qC) = tblCell.(q{1})(isC);
    wide.(qM) = tblCell.(q{1})(isM);
end
wide.genotype = setcats(wide.genotype, {'Control', 'MCU-KO'});

rows = cell(0, 11);
for iM = 1:numel(quants)
    q = quants{iM};
    xCol = ['cyto' upper(q(1)) q(2:end)];
    yCol = ['mito' upper(q(1)) q(2:end)];
    wFit = wide;
    wFit.(xCol) = log(wide.(xCol));
    wFit.(yCol) = log(wide.(yCol));
    ok = isfinite(wFit.(xCol)) & isfinite(wFit.(yCol));
    wFit = wFit(ok, :);
    if height(wFit) < 4, continue; end
    mdl = fitlm(wFit, sprintf('%s ~ %s * genotype', yCol, xCol));
    rr  = relOLSCoefs(mdl, xCol);
    rows(end+1, :) = {sprintf('mito.%s ~ cyto.%s', q, q), ...
        rr.bMain, rr.pMain, rr.bInter, rr.pInter, rr.bKO, ...
        rr.F, rr.df1, rr.df2, rr.p, mdl.Rsquared.Ordinary}; %#ok<AGROW>
end
tbl = cell2table(rows, 'VariableNames', ...
    {'model', 'beta_Ctrl', 'p_main', 'beta_inter', 'p_inter', 'beta_KO', ...
     'F_geno', 'df1', 'df2', 'p_geno', 'R2'});
end


%% ========================================================================
%  T4: event-level within-compartment paired scaling (LME, log-log)
%  ========================================================================

function tbl = relEvent(tblEvent, dist)
scaling = { 'pairFlux',  'flux'; ...
            'pairAmp',   'amp';  ...
            'fluxOther', 'flux' };
comps   = {'Cyto', 'Mito'};
rows = cell(0, 11);
for iM = 1:size(scaling, 1)
    yV = scaling{iM, 1};
    xV = scaling{iM, 2};
    if ~all(ismember({yV, xV}, tblEvent.Properties.VariableNames)), continue; end
    for iC = 1:numel(comps)
        sub  = tblEvent(tblEvent.compartment == comps{iC}, :);
        frml = sprintf('%s ~ %s * genotype + (1|sbjID)', yV, xV);
        rr   = fitScaleLME(sub, frml, yV, xV, dist);
        rows(end+1, :) = {sprintf('%s ~ %s', yV, xV), comps{iC}, ...
            rr.bMain, rr.pMain, rr.bInter, rr.pInter, rr.bKO, ...
            rr.F, rr.df1, rr.df2, rr.p}; %#ok<AGROW>
    end
end
tbl = cell2table(rows, 'VariableNames', ...
    {'model', 'compartment', 'beta_main', 'p_main', 'beta_inter', 'p_inter', ...
     'beta_KO', 'F_geno', 'df1', 'df2', 'p_geno'});
end


%% ========================================================================
%  HELPERS - FIT + EXTRACT
%  ========================================================================

function r = fitGenoOLS(tbl, frml, yVar, dist)
r = emptyGenoResult();
isLog = strcmpi(dist, 'Log-Normal');
if isLog
    keep = tbl.(yVar) > 1e-8 & isfinite(tbl.(yVar));
    tbl  = tbl(keep, :);
end
if height(tbl) < 4, return; end
tblFit = tbl;
if isLog
    tblFit.(yVar) = log(tbl.(yVar));
    r.transform = 'log';
else
    r.transform = 'none';
end
try
    mdl = fitlm(tblFit, frml);
    c = mdl.Coefficients;
    nm = string(c.Properties.RowNames);
    idx = find(contains(nm, 'genotype') & ~contains(nm, ':'), 1);
    if isempty(idx), return; end
    r.beta = c.Estimate(idx);
    r.SE   = c.SE(idx);
    r.p    = c.pValue(idx);
catch
end
end


function r = fitGenoLME(tbl, frml, yVar, dist)
r = emptyGenoResult();
isLog = strcmpi(dist, 'Log-Normal');
if isLog
    keep = tbl.(yVar) > 1e-8 & isfinite(tbl.(yVar));
    tbl  = tbl(keep, :);
end
if height(tbl) < 4, return; end
try
    [mdl, ~, ~] = lme_analyse(tbl, frml, ...
        'dist', dist, 'flgPlot', false, 'verbose', false, 'flgStnd', false);
    r.transform = ternary(isLog, 'log', 'none');
    c = mdl.Coefficients;
    if istable(c)
        nm = string(c.Properties.RowNames);
    else
        nm = string(c.Name);
    end
    idx = find(contains(nm, 'genotype') & ~contains(nm, ':'), 1);
    if isempty(idx), return; end
    r.beta = c.Estimate(idx);
    r.SE   = c.SE(idx);
    r.p    = c.pValue(idx);
catch
end
end


function rr = relOLSCoefs(mdl, xCol)
% Pull main slope, interaction slope, and joint F on all genotype coefs.
rr = emptyRelResult();
c = mdl.Coefficients;
nm = string(c.Properties.RowNames);
iM = find(nm == xCol, 1);
iI = find(contains(nm, xCol) & contains(nm, "genotype"), 1);
if ~isempty(iM)
    rr.bMain = c.Estimate(iM);
    rr.pMain = c.pValue(iM);
end
if ~isempty(iI)
    rr.bInter = c.Estimate(iI);
    rr.pInter = c.pValue(iI);
    rr.bKO    = rr.bMain + rr.bInter;
end
genoIdx = find(contains(nm, 'genotype'));
if isempty(genoIdx), return; end
H = zeros(numel(genoIdx), numel(nm));
for k = 1:numel(genoIdx), H(k, genoIdx(k)) = 1; end
try
    [p, F, df1] = coefTest(mdl, H);
    rr.F   = F;
    rr.df1 = df1;
    rr.df2 = mdl.DFE;
    rr.p   = p;
catch
end
end


function rr = fitScaleLME(tbl, frml, yVar, xVar, dist)
rr = emptyRelResult();
isLog = strcmpi(dist, 'Log-Normal');
if isLog
    keep = tbl.(yVar) > 1e-8 & isfinite(tbl.(yVar)) & ...
           tbl.(xVar) > 1e-8 & isfinite(tbl.(xVar));
    tbl  = tbl(keep, :);
end
if height(tbl) < 8, return; end
try
    [mdl, ~, ~] = lme_analyse(tbl, frml, ...
        'dist', dist, 'flgPlot', false, 'verbose', false, 'flgStnd', false);
    c = mdl.Coefficients;
    if istable(c)
        nm = string(c.Properties.RowNames);
    else
        nm = string(c.Name);
    end
    iM = find(nm == xVar, 1);
    iI = find(contains(nm, xVar) & contains(nm, ":"), 1);
    if ~isempty(iM)
        rr.bMain = c.Estimate(iM);
        rr.pMain = c.pValue(iM);
    end
    if ~isempty(iI)
        rr.bInter = c.Estimate(iI);
        rr.pInter = c.pValue(iI);
        rr.bKO    = rr.bMain + rr.bInter;
    end
    genoIdx = find(contains(nm, 'genotype'));
    if isempty(genoIdx), return; end
    H = zeros(numel(genoIdx), numel(nm));
    for k = 1:numel(genoIdx), H(k, genoIdx(k)) = 1; end
    [p, F, df1] = coefTest(mdl, H);
    rr.F   = F;
    rr.df1 = df1;
    rr.df2 = mdl.DFE;
    rr.p   = p;
catch
end
end


function r = emptyGenoResult()
r = struct('beta', NaN, 'SE', NaN, 'p', NaN, 'transform', '');
end


function r = emptyRelResult()
r = struct('bMain', NaN, 'pMain', NaN, 'bInter', NaN, 'pInter', NaN, ...
    'bKO', NaN, 'F', NaN, 'df1', NaN, 'df2', NaN, 'p', NaN);
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end


%% ========================================================================
%  HELPERS - PRINT
%  ========================================================================

function printAll(stats)
fprintf('\n### T1 - Cell-level genotype contrasts (KO vs Control)\n\n');
printGenoTable(stats.cellGeno);

fprintf('\n### T2 - Event-level genotype contrasts (KO vs Control)\n\n');
printGenoTable(stats.eventGeno);

fprintf('\n### T3 - Cell-level cross-compartment relationships (compensation)\n\n');
printRelTable(stats.cellRel);

fprintf('\n### T4 - Event-level within-compartment paired scaling\n\n');
printRelTable(stats.eventRel);
end


function printGenoTable(t)
fprintf('| quantity | compartment | beta(KO-Ctrl) | SE | p | sig | transform |\n');
fprintf('|---|---|---:|---:|---:|:---:|:---:|\n');
for iR = 1:height(t)
    fprintf('| %s | %s | %s | %s | %s | %s | %s |\n', ...
        t.quantity{iR}, t.compartment{iR}, ...
        fmtBeta(t.beta_KOvsCtrl(iR)), fmtBeta(t.SE(iR)), ...
        fmtP(t.p(iR)), sigStars(t.p(iR)), ...
        cellChar(t.transform(iR)));
end
end


function printRelTable(t)
hasComp = ismember('compartment', t.Properties.VariableNames);
if hasComp
    fprintf('| model | compartment | beta_main | p_main | beta_inter | p_inter | beta_KO | F(geno) | df | p(geno) | sig(geno) |\n');
    fprintf('|---|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|\n');
else
    fprintf('| model | beta_main | p_main | beta_inter | p_inter | beta_KO | F(geno) | df | p(geno) | sig(geno) | R2 |\n');
    fprintf('|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---:|\n');
end
for iR = 1:height(t)
    df_str = sprintf('%g,%g', t.df1(iR), t.df2(iR));
    if hasComp
        fprintf('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', ...
            t.model{iR}, t.compartment{iR}, ...
            fmtBeta(t.beta_main(iR)), fmtP(t.p_main(iR)), ...
            fmtBeta(t.beta_inter(iR)), fmtP(t.p_inter(iR)), ...
            fmtBeta(t.beta_KO(iR)), fmtBeta(t.F_geno(iR)), ...
            df_str, fmtP(t.p_geno(iR)), sigStars(t.p_geno(iR)));
    else
        fprintf('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %.3f |\n', ...
            t.model{iR}, ...
            fmtBeta(t.beta_Ctrl(iR)), fmtP(t.p_main(iR)), ...
            fmtBeta(t.beta_inter(iR)), fmtP(t.p_inter(iR)), ...
            fmtBeta(t.beta_KO(iR)), fmtBeta(t.F_geno(iR)), ...
            df_str, fmtP(t.p_geno(iR)), sigStars(t.p_geno(iR)), ...
            t.R2(iR));
    end
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
if isnan(p),     s = '—';
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
