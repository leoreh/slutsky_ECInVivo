function h = gui_viewPanel(parent, api)
% GUI_VIEWPANEL  Window navigator for the view (no-target) mode of guiPath.
%
%   h = gui_viewPanel(parent, api) builds a Prev/Next block that steps the
%   window across the session when nothing is curated (CURATE = None). It is the
%   sibling of gui_eventPanel / gui_statePanel and shares their drop-in contract
%   (H.grid + H.refresh), but carries no accept/reject/assign/save: with no
%   target there is nothing to edit, only a window to move.
%
%   The host owns all state. API is a struct of callbacks:
%       api.prev()   - step the window back  (by ~one window width)
%       api.next()   - step the window forward
%   Any missing field is ignored (the control becomes a no-op).
%
%   The returned struct H carries the handles and:
%       H.refresh(t0, Tend)   - update the position readout [s].
%
%   HISTORY:
%       17 Jul 2026 - created for the CURATE = None view mode (window stepping).

narginchk(2, 2)
if ~isstruct(api), error('gui_viewPanel:api', 'api must be a struct of callbacks'); end

g = uigridlayout(parent, [2, 2], 'RowHeight', {'fit', 'fit'}, ...
    'ColumnWidth', {'1x', '1x'}, 'Padding', 4, 'RowSpacing', 4, 'ColumnSpacing', 4);

% row 1: prev / next  (arrows match the host keyboard shortcuts)
bPrev = uibutton(g, 'Text', 'Prev (<-)', 'ButtonPushedFcn', @(~, ~) safecall(api, 'prev'));
bPrev.Layout.Row = 1; bPrev.Layout.Column = 1;
bNext = uibutton(g, 'Text', 'Next (->)', 'ButtonPushedFcn', @(~, ~) safecall(api, 'next'));
bNext.Layout.Row = 1; bNext.Layout.Column = 2;

% row 2: position readout (spans both columns)
hStat = uilabel(g, 'Text', '', 'WordWrap', 'on', 'FontColor', [0.3 0.3 0.3]);
hStat.Layout.Row = 2; hStat.Layout.Column = [1, 2];

h = struct('grid', g, 'status', hStat);
h.refresh = @(t0, Tend) refreshPanel(h, t0, Tend);

end     % MAIN

% ------------------------------------------------------------------------
function refreshPanel(h, t0, Tend)
% show the window centre against the session length
if ~isvalid(h.status), return; end
h.status.Text = sprintf('t0 = %.2f / %.0f s', t0, Tend);
end

% ------------------------------------------------------------------------
function safecall(api, name, varargin)
% invoke api.(name) if present and non-empty; otherwise no-op
if isfield(api, name) && ~isempty(api.(name))
    api.(name)(varargin{:});
end
end

% EOF
