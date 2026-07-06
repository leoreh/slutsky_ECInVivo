function h = gui_labeledControl(grid, kind, labelText, varargin)
% GUI_LABELEDCONTROL  Append a control (optionally label + control) to a
% one-column vertical uigridlayout, growing the grid by a row per element.
%
%   h = gui_labeledControl(grid, kind, labelText, ...) adds an optional
%   bold label row (when labelText is non-empty) followed by the control, and
%   returns the control handle. It hides the legacy-to-uifigure component
%   mapping so callers never name uidropdown / uibutton / uieditfield etc.
%
%   kind (case-insensitive):
%       'label'    a standalone uilabel (uses labelText as its text)
%       'dropdown' uidropdown
%       'checkbox' uicheckbox
%       'button'   uibutton (push)
%       'toggle'   uibutton('state', ...)
%       'edit'     uieditfield('text')
%       'editnum'  uieditfield('numeric')
%       'panel'    uipanel (BorderType none) - e.g. holds a gui_filterPanel
%       'spacer'   an empty '1x' row that pushes preceding rows to the top
%
%   Name-value:
%       'RowHeight'  height of the control's row {'fit'}. Use '1x' for an
%                    expanding region such as a filter-checkbox panel.
%   Remaining name-value pairs are forwarded to the created component.

% Pull the RowHeight option out of varargin (default 'fit').
rh = 'fit';
keys = varargin(1:2:end);
ix = find(strcmpi(keys, 'RowHeight'), 1);
if ~isempty(ix)
    rh = varargin{2 * ix};
    varargin([2 * ix - 1, 2 * ix]) = [];
end

if strcmpi(kind, 'spacer')
    appendRow(grid, '1x');
    h = gobjects(0);
    return;
end

% Optional label row above the control.
if ~isempty(labelText) && ~strcmpi(kind, 'label')
    rL = appendRow(grid, 'fit');
    t = uilabel(grid, 'Text', labelText, 'FontWeight', 'bold');
    t.Layout.Row = rL;
    t.Layout.Column = 1;
end

r = appendRow(grid, rh);
switch lower(kind)
    case 'label',    h = uilabel(grid, 'Text', labelText, varargin{:});
    case 'dropdown', h = uidropdown(grid, varargin{:});
    case 'checkbox', h = uicheckbox(grid, varargin{:});
    case 'button',   h = uibutton(grid, varargin{:});
    case 'toggle',   h = uibutton(grid, 'state', varargin{:});
    case 'edit',     h = uieditfield(grid, 'text', varargin{:});
    case 'editnum',  h = uieditfield(grid, 'numeric', varargin{:});
    case 'panel',    h = uipanel(grid, 'BorderType', 'none', varargin{:});
    otherwise
        error('gui_labeledControl:kind', 'Unknown kind "%s".', kind);
end
h.Layout.Row = r;
h.Layout.Column = 1;

end

function r = appendRow(grid, h)
% Grow the grid's RowHeight by one and return the new row index. A running
% count is kept on the grid's UserData (the figure holds GUI state elsewhere).
n = 0;
if isstruct(grid.UserData) && isfield(grid.UserData, 'nRow')
    n = grid.UserData.nRow;
end
n = n + 1;
if n == 1
    grid.RowHeight = {h};
else
    grid.RowHeight = [grid.RowHeight, {h}];
end
ud = grid.UserData;
if ~isstruct(ud), ud = struct(); end
ud.nRow = n;
grid.UserData = ud;
r = n;
end     % EOF
