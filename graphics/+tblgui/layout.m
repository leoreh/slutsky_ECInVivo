function [gl, gPlot, gCtrl, gActions] = layout(parent, varargin)
% TBLGUI.LAYOUT  Standard two-column GUI scaffold (controls | plot).
%
%   [gl, gPlot, gCtrl, gActions] = tblgui.layout(parent) builds a 1x2
%   uigridlayout inside PARENT and returns:
%       gl        the top-level grid
%       gPlot     the plot-side cell (host a uiaxes, nested grid, or tiledlayout)
%       gCtrl     a scrollable vertical grid for stacked controls
%                 (append rows with tblgui.labeledControl)
%       gActions  a fixed grid pinned at the bottom of the control column for
%                 action buttons (e.g. Select / Save / "Push Units"), so they
%                 never collide with the scrollable controls above
%
%   PARENT must be a uifigure or a uifigure container (uigridlayout cannot be
%   parented into a legacy figure).
%
%   Name-value:
%       'CtrlWidth'  control-column width in pixels {220}
%       'CtrlSide'   'left' (default) or 'right'

p = inputParser;
addParameter(p, 'CtrlWidth', 220, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'CtrlSide', 'left', @(x) any(strcmpi(x, {'left', 'right'})));
parse(p, varargin{:});
w = p.Results.CtrlWidth;

if strcmpi(p.Results.CtrlSide, 'left')
    colW = {w, '1x'};  cCtrl = 1;  cPlot = 2;
else
    colW = {'1x', w};  cCtrl = 2;  cPlot = 1;
end

gl = uigridlayout(parent, [1, 2], 'ColumnWidth', colW, 'RowHeight', {'1x'}, ...
    'Padding', 4, 'ColumnSpacing', 6);

% Control column: scrollable controls on top, fixed action area at the bottom.
gCol = uigridlayout(gl, [2, 1], 'RowHeight', {'1x', 'fit'}, ...
    'Padding', 0, 'RowSpacing', 4);
gCol.Layout.Row = 1;
gCol.Layout.Column = cCtrl;

gCtrl = uigridlayout(gCol, [1, 1], 'RowHeight', {'fit'}, 'ColumnWidth', {'1x'}, ...
    'Padding', 4, 'RowSpacing', 4, 'Scrollable', 'on');
gCtrl.Layout.Row = 1;

gActions = uigridlayout(gCol, [1, 1], 'RowHeight', {'fit'}, 'ColumnWidth', {'1x'}, ...
    'Padding', 4, 'RowSpacing', 4);
gActions.Layout.Row = 2;

gPlot = uigridlayout(gl, [1, 1], 'Padding', 0);
gPlot.Layout.Row = 1;
gPlot.Layout.Column = cPlot;

end     % EOF
