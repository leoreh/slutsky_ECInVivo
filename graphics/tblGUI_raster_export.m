function tblGUI_raster_export(figH)
% TBLGUI_RASTER_EXPORT  Save a publication-quality vector figure of the current raster view.
%
%   TBLGUI_RASTER_EXPORT(figH)
%
%   SUMMARY:
%       Reads the render data cached in the tblGUI_raster figure's
%       UserData, reconstructs the raster in a clean figure with no UI
%       components, places the axis at an exact centimeter position, and
%       exports as both PDF and EPS.
%
%       The exported figure can be opened and edited in Adobe Illustrator
%       (File > Open the PDF or EPS) because all graphics are vector.
%
%   AXIS DIMENSIONS:
%       Width  : 6 cm
%       Height : 3.5 cm
%
%   TYPOGRAPHY:
%       Font name   : Arial
%       Tick labels  : 10 pt
%       Axis labels  : 12 pt
%
%   INPUTS:
%       figH  - (Handle) Figure handle returned by tblGUI_raster.
%
%   OUTPUTS (written to current directory):
%       raster_export.pdf
%       raster_export.eps
%
%   DEPENDENCIES:
%       plot_raster.
%
%   HISTORY:
%       Created: 29 Mar 2026

%% ========================================================================
%  READ GUI STATE
%  ========================================================================

data = figH.UserData;
rd   = data.renderData;

if isempty(fieldnames(rd))
    warning('[TBLGUI_RASTER_EXPORT]: No render data found. Plot must be drawn first.');
    return;
end


%% ========================================================================
%  FIGURE GEOMETRY  [all units: centimeters]
%  ========================================================================

% Axis dimensions
axW = 6.0;     % width
axH = 3.5;     % height

% Margins
lMarg = 1.5;   % left   (room for y-labels)
rMarg = 0.3;   % right
bMarg = 1.2;   % bottom (room for x-label)
tMarg = 0.3;   % top

% Figure size
figW = lMarg + axW + rMarg;
figHt = bMarg + axH + tMarg;


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
hExpFig.Position(3:4) = [figW, figHt];

hAx = axes('Units', 'centimeters', ...
    'Position', [lMarg, bMarg, axW, axH]); %#ok<LAXES>
set(hAx, 'FontName', fontName, 'FontSize', fontSzTick);
hold(hAx, 'on');

% All spikes (black)
if ~isempty(rd.spikes)
    plot_raster(rd.spikes, 'hAx', hAx, ...
        'plotType', 'vertline', 'clr', [0 0 0], ...
        'xLim', [rd.xLo, rd.xHi], 'flgLbls', false);
end

% Burst spikes overlay (red)
if rd.flgBrst && any(~cellfun('isempty', rd.brstSpks))
    plot_raster(rd.brstSpks, 'hAx', hAx, ...
        'plotType', 'vertline', 'clr', [1 0 0], ...
        'xLim', [rd.xLo, rd.xHi], 'flgLbls', false);
end

hold(hAx, 'off');

% Axis limits and formatting
xlim(hAx, [rd.xLo, rd.xHi]);
set(hAx, 'YDir', 'normal');
xlabel(hAx, 'Time (s)', 'FontName', fontName, 'FontSize', fontSzLbl);
ylabel(hAx, 'Unit No.', 'FontName', fontName, 'FontSize', fontSzLbl);


%% ========================================================================
%  SAVE
%  ========================================================================

filePDF = fullfile(pwd, 'raster_export.pdf');
fileEPS = fullfile(pwd, 'raster_export.eps');

exportgraphics(hExpFig, filePDF, 'ContentType', 'vector', 'BackgroundColor', 'white');
exportgraphics(hExpFig, fileEPS, 'ContentType', 'vector', 'BackgroundColor', 'white');

fprintf('[TBLGUI_RASTER_EXPORT]: Saved to:\n   %s\n   %s\n', filePDF, fileEPS);

end     % EOF
