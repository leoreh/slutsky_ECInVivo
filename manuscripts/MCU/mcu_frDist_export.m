function kdMat = mcu_frDist_export(tbl, outPath)
% MCU_FRDIST_EXPORT  Publication figure of RS FR distributions across days.
%
%   KDMAT = MCU_FRDIST_EXPORT(TBL, OUTPATH)
%
%   SUMMARY:
%       Two-panel vector figure (one panel per genotype) overlaying kernel-
%       density estimates of per-unit firing rates across days on a log
%       frequency axis. Exports as PDF and EPS to OUTPATH.
%
%   INPUTS:
%       tbl      - (table) Must contain fields: fr, genotype, day.
%       outPath  - (char)  Output directory.
%
%   OUTPUTS:
%       kdMat    - (matrix) Column 1 is the FR evaluation grid in Hz.
%                  Subsequent columns are the peak-normalized KDEs,
%                  ordered genotype-outer / day-inner (Control days
%                  first, then MCU-KO days), e.g. Control_BSL,
%                  Control_BAC1, ..., MCU-KO_BSL, MCU-KO_BAC1, ....
%                  Missing day x genotype combinations are NaN. Paste
%                  directly into a GraphPad Prism XY dataset.
%
%   AXIS DIMENSIONS:
%       Width  : 5 cm per panel
%       Height : 4 cm
%
%   TYPOGRAPHY:
%       Font name   : Arial
%       Tick labels : 10 pt
%       Axis labels : 12 pt
%
%   DEPENDENCIES:
%       plot_hist.
%
%   HISTORY:
%       Created: 21 Apr 2026

%% ========================================================================
%  FIGURE GEOMETRY  [units: centimeters]
%  ========================================================================

axW   = 5.0;
axH   = 4.0;
lMarg = 1.5;
rMarg = 0.3;
bMarg = 1.2;
tMarg = 0.8;
hGap  = 1.2;


%% ========================================================================
%  TYPOGRAPHY
%  ========================================================================

fontName   = 'Arial';
fontSzTick = 10;
fontSzLbl  = 12;


%% ========================================================================
%  DATA PREPARATION
%  ========================================================================

grps = categories(removecats(tbl.genotype));
days = categories(removecats(tbl.day));
nGrp = numel(grps);
nDay = numel(days);

% Per-day palette: BSL near-black (reference); BAC1-3 purple -> blue
colDay = [ ...
    0.10 0.10 0.10;       % BSL
    0.55 0.20 0.65;       % BAC1 - purple
    0.30 0.35 0.75;       % BAC2 - indigo
    0.15 0.55 0.85;       % BAC3 - blue
    0.60 0.60 0.60];      % WASH - gray

if nDay > size(colDay, 1)
    colDay = [colDay; lines(nDay - size(colDay, 1))];
else
    colDay = colDay(1:nDay, :);
end

% Common bin edges so panels are directly comparable
validFR = tbl.fr(~isnan(tbl.fr) & tbl.fr > 0);
bins    = logspace(log10(min(validFR)), log10(max(validFR)), 40);


%% ========================================================================
%  CREATE FIGURE
%  ========================================================================

figW = lMarg + nGrp * axW + (nGrp - 1) * hGap + rMarg;
figH = bMarg + axH + tMarg;

hFig = figure('Visible', 'on', 'Color', 'w');
set(hFig, 'Units', 'centimeters');
hFig.Position(3:4) = [figW, figH];

hAx  = gobjects(nGrp, 1);
hKDE = gobjects(nDay, 1);
kdX  = [];
kdY  = [];

for iGrp = 1:nGrp
    xL = lMarg + (iGrp - 1) * (axW + hGap);
    hAx(iGrp) = axes('Units', 'centimeters', ...
        'Position', [xL, bMarg, axW, axH]); %#ok<LAXES>
    set(hAx(iGrp), 'FontName', fontName, 'FontSize', fontSzTick);
    hold(hAx(iGrp), 'on');

    subT = tbl(tbl.genotype == grps{iGrp}, :);
    [hH, hK, ~] = plot_hist(subT, 'fr', ...
        'g', 'day', 'c', colDay, 'bins', bins, ...
        'scale', 'log', 'flgKDE', true, 'flgStat', false, ...
        'alpha', 0, 'hAx', hAx(iGrp));

    % Drop histogram patches (kept only for bin scaffolding) so axis
    % autoscale is driven by the KDE curves alone.
    delete(hH(isgraphics(hH)));

    % Peak-normalize each KDE curve to 1 (so overlaid distributions are
    % compared on shape/position rather than absolute density) and capture
    % XData/YData into kdMat for Prism export. Map this subset's day
    % index into the global day order so missing day x genotype
    % combinations land as NaN columns.
    subDays = categories(removecats(subT.day));
    [~, dayIdx] = ismember(subDays, days);
    for iL = 1:numel(hK)
        if isgraphics(hK(iL))
            yd = get(hK(iL), 'YData');
            yd = yd(:) / max(yd);
            set(hK(iL), 'YData', yd);
            if isempty(kdX)
                xd  = get(hK(iL), 'XData');
                kdX = xd(:);
                kdY = nan(numel(kdX), nGrp * nDay);
            end
            col = (iGrp - 1) * nDay + dayIdx(iL);
            kdY(:, col) = yd;
        end
    end

    if iGrp == nGrp, hKDE = hK; end

    set(hAx(iGrp), 'XScale', 'log');
    xlabel(hAx(iGrp), 'Firing rate (Hz)', ...
        'FontName', fontName, 'FontSize', fontSzLbl);
    ylabel(hAx(iGrp), 'Normalized density', ...
        'FontName', fontName, 'FontSize', fontSzLbl);
    title(hAx(iGrp), grps{iGrp}, ...
        'FontName', fontName, 'FontSize', fontSzLbl, 'FontWeight', 'normal');

    xlim(hAx(iGrp), [0.01, 18])
    xticks(hAx(iGrp), [0.01, 0.1, 1, 10])

    ylim(hAx(iGrp), [0, 1.05]);

    % Legend (after loop, so hAx and hKDE are both fully populated)
    validK = isgraphics(hK);
    hLgd = legend(hAx(iGrp), hK(validK), days(validK), ...
        'Location', 'northwest', 'Box', 'off', ...
        'FontName', fontName, 'FontSize', fontSzTick);
    hLgd.ItemTokenSize = [10, 18];

end

% Assemble Prism-ready matrix and print column labels
if isempty(kdX)
    kdMat = [];
else
    kdMat = [kdX, kdY];
    fprintf('[MCU_FRDIST_EXPORT] kdMat columns:\n   1: FR (Hz)\n');
    col = 1;
    for iGrp = 1:nGrp
        for iDay = 1:nDay
            col = col + 1;
            fprintf('   %2d: %s_%s\n', col, grps{iGrp}, days{iDay});
        end
    end
end


%% ========================================================================
%  EXPORT
%  ========================================================================

% tag     = 'mcu_frDist';
% filePDF = fullfile(outPath, [tag, '.pdf']);
% fileEPS = fullfile(outPath, [tag, '.eps']);

% exportgraphics(hFig, filePDF, 'ContentType', 'vector', 'BackgroundColor', 'white');
% exportgraphics(hFig, fileEPS, 'ContentType', 'vector', 'BackgroundColor', 'white');

% fprintf('[MCU_FRDIST_EXPORT]: Saved to:\n   %s\n   %s\n', filePDF, fileEPS);

end     % EOF
