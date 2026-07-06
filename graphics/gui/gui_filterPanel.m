function chk = gui_filterPanel(parent, cats, cb, varargin)
% GUI_FILTERPANEL  Build a column of category checkboxes in a container.
%
%   chk = gui_filterPanel(parent, cats, cb) clears PARENT (a uipanel or
%   other uifigure container) and fills it with one uicheckbox per category
%   in CATS (cellstr), each wired to the ValueChanged callback CB. Returns the
%   checkbox handles (1 x nCats, gobjects when empty). The checkboxes live in
%   a scrollable grid so long category lists do not overflow.
%
%   Name-value:
%       'InitVal'  logical/numeric vector of initial checked states {all true}

p = inputParser;
addParameter(p, 'InitVal', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
parse(p, varargin{:});
initVal = p.Results.InitVal;

delete(allchild(parent));
cats = cellstr(cats);
n = numel(cats);
chk = gobjects(1, n);
if n == 0
    return;
end

if isempty(initVal), initVal = true(1, n); end
initVal = logical(initVal);

g = uigridlayout(parent, [n + 1, 1], ...
    'RowHeight', [repmat({22}, 1, n), {'1x'}], 'ColumnWidth', {'1x'}, ...
    'Padding', [2, 2, 2, 2], 'RowSpacing', 2, 'Scrollable', 'on');

for k = 1:n
    chk(k) = uicheckbox(g, 'Text', cats{k}, 'Value', initVal(k), ...
        'ValueChangedFcn', cb);
    chk(k).Layout.Row = k;
    chk(k).Layout.Column = 1;
end

end     % EOF
