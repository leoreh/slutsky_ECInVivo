function ripp_gui2_export(figH)
% RIPP_GUI2_EXPORT  Save a publication-quality vector figure of the current SWR event.
%
%   RIPP_GUI2_EXPORT(figH)
%
%   SUMMARY:
%       Reads the render data cached in the ripp_gui2 figure's UserData,
%       reconstructs the three-panel view in a clean figure with no UI
%       components, places each axis at an exact centimeter position, and
%       exports the result as both a PDF and an EPS to basepath.
%
%       The exported figure can be opened and edited in Adobe Illustrator
%       (File → Open the PDF or EPS) because all graphics are vector.
%
%   AXIS DIMENSIONS:
%       Width           : 6 cm  (all three panels)
%       LFP height      : 3 cm
%       Raster height   : 3 cm
%       KDE strip height: 1.5 cm
%
%   TYPOGRAPHY:
%       Font name: Arial
%       Tick labels: 10 pt
%       Axis labels: 12 pt
%
%   INPUTS:
%       figH  - (Handle) Figure handle returned by ripp_gui2.
%
%   OUTPUTS (written to basepath):
%       <basename>_ripp<idx>.pdf
%       <basename>_ripp<idx>.eps
%
%   DEPENDENCIES:
%       plot_raster.
%
%   HISTORY:
%       Created: 24 Mar 2026

%% ========================================================================
%  READ GUI STATE
%  ========================================================================

data = figH.UserData;
rd   = data.renderData;

if isempty(fieldnames(rd))
    warning('[RIPP_GUI2_EXPORT]: No render data found. Navigate to an event first.');
    return;
end


%% ========================================================================
%  FIGURE GEOMETRY  [all units: centimeters]
%  ========================================================================

% Axis dimensions
axW    = 6.0;    % width (all panels)
axHraw = 3.0;    % LFP panel height
axHspk = 3.0;    % Raster panel height
axHmua = 1.5;    % Spike-density strip height

% Margins and inter-axis gap
lMarg = 1.5;    % left   (room for y-labels)
rMarg = 0.3;    % right
bMarg = 1.2;    % bottom (room for x-label + CoM triangle)
tMarg = 0.5;    % top
gap   = 0.25;   % vertical gap between adjacent panels

% Derived figure size
figW = lMarg + axW + rMarg;
figH = bMarg + axHmua + gap + axHspk + gap + axHraw + tMarg;

% Axis bottom-left corners [x, y] — origin at bottom-left of figure
xLeft  = lMarg;
yMua   = bMarg;
ySpk   = bMarg + axHmua + gap;
yRaw   = bMarg + axHmua + gap + axHspk + gap;


%% ========================================================================
%  TYPOGRAPHY
%  ========================================================================

fontName   = 'Arial';
fontSzTick = 10;      % tick-label font size [pt]
fontSzLbl  = 12;      % axis-label font size [pt]


%% ========================================================================
%  CREATE EXPORT FIGURE
%  ========================================================================

hExpFig = figure('Visible', 'on', 'Color', 'w');
set(hExpFig, 'Units', 'centimeters');
hExpFig.Position(3:4) = [figW, figH];

% Helper: set axis to exact cm position and apply typography
    function hAx = makeAx(posVec)
        hAx = axes('Units', 'centimeters', 'Position', posVec); %#ok<LAXES>
        set(hAx, 'FontName', fontName, 'FontSize', fontSzTick);
        hold(hAx, 'on');
    end

hAxRaw = makeAx([xLeft, yRaw, axW, axHraw]);
hAxSpk = makeAx([xLeft, ySpk, axW, axHspk]);
hAxMua = makeAx([xLeft, yMua, axW, axHmua]);


%% ========================================================================
%  PANEL 1: RAW LFP
%  ========================================================================

plot(hAxRaw, rd.tVec, rd.lfpSeg, 'k');
xline(hAxRaw, [rd.stT, rd.pkT, rd.enT], '--b');
xlim(hAxRaw, [rd.winStart, rd.winEnd]);
ylim(hAxRaw, rd.yLimRaw);
grid(hAxRaw, 'on');
set(hAxRaw, 'XTickLabel', []);   % shared x-axis: suppress tick labels
ylabel(hAxRaw, 'LFP (\muV)', 'FontName', fontName, 'FontSize', fontSzLbl);


%% ========================================================================
%  PANEL 2: RASTER (firing-order sorted, as displayed in the GUI)
%  ========================================================================

% All spikes — black
plot_raster(rd.spkSorted, 'hAx', hAxSpk, ...
    'xLim', [rd.winStart, rd.winEnd], 'flgLbls', false, 'clr', [0 0 0]);

% Burst spikes overlaid — red
if any(~cellfun(@isempty, rd.brstSorted))
    plot_raster(rd.brstSorted, 'hAx', hAxSpk, ...
        'xLim', [rd.winStart, rd.winEnd], 'flgLbls', false, 'clr', [1 0 0]);
end

xline(hAxSpk, [rd.stT, rd.pkT, rd.enT], '--b');
xlim(hAxSpk, [rd.winStart, rd.winEnd]);
set(hAxSpk, 'XTickLabel', []);   % shared x-axis: suppress tick labels
ylabel(hAxSpk, 'Units (RS)', 'FontName', fontName, 'FontSize', fontSzLbl);


%% ========================================================================
%  PANEL 3: SPIKE-DENSITY STRIP
%  ========================================================================

% Ripple extent shading
patch(hAxMua, [rd.stT, rd.enT, rd.enT, rd.stT], [0, 0, 1, 1], ...
    [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Normalized KDE curve
if ~isempty(rd.muaKDE)
    plot(hAxMua, rd.kdePts, rd.muaKDE, 'k', 'LineWidth', 1.5);
end

% SWR peak
xline(hAxMua, rd.pkT, '--b');

% CoM marker — orange upward triangle
if ~isnan(rd.com)
    plot(hAxMua, rd.com, 0, '^', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', [0.85, 0.33, 0.10], ...
        'MarkerEdgeColor', 'none', ...
        'Clipping', 'off');
end

xlim(hAxMua, [rd.winStart, rd.winEnd]);
ylim(hAxMua, [0, 1.15]);
grid(hAxMua, 'on');
xlabel(hAxMua, 'Time (s)',       'FontName', fontName, 'FontSize', fontSzLbl);
ylabel(hAxMua, 'Spike Density',  'FontName', fontName, 'FontSize', fontSzLbl);


%% ========================================================================
%  FINALISE & SAVE
%  ========================================================================

% Link x-axes so zooming in Illustrator is consistent
linkaxes([hAxRaw, hAxSpk, hAxMua], 'x');

% File names
tag      = sprintf('%s_ripp%04d', data.basename, rd.idx);
filePDF  = fullfile(data.basepath, [tag, '.pdf']);
fileEPS  = fullfile(data.basepath, [tag, '.eps']);

% Export as vector (requires R2020a+)
exportgraphics(hExpFig, filePDF, 'ContentType', 'vector', 'BackgroundColor', 'white');
exportgraphics(hExpFig, fileEPS, 'ContentType', 'vector', 'BackgroundColor', 'white');

fprintf('[RIPP_GUI2_EXPORT]: Saved to:\n   %s\n   %s\n', filePDF, fileEPS);

end     % EOF
