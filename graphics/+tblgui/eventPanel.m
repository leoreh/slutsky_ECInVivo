function h = eventPanel(parent, api)
% TBLGUI.EVENTPANEL  Compact, swappable event-curation stepper.
%
%   h = tblgui.eventPanel(parent, api) builds a self-contained control block
%   for stepping through and triaging a list of events: an event-index field
%   (+ "/ N"), Prev/Next, Accept/Reject, Save, and a status line. It is meant
%   to drop into the pinned action area of tblgui.layout (or any uifigure
%   container) so different event types (EDs, ripples, states) can swap the
%   same widget without touching the host viewer.
%
%   The host owns all state. API is a struct of callbacks invoked on action:
%       api.prev()        - step to the previous event
%       api.next()        - step to the next event
%       api.accept()      - accept the current event
%       api.reject()      - reject the current event
%       api.save()        - persist the curation
%       api.setIdx(v)     - jump to 1-based event index v
%   Any missing field is simply ignored (the control becomes a no-op), so a
%   host can wire only the callbacks it needs.
%
%   The returned struct H carries the control handles and:
%       H.refresh(idx, total, accepted, statusText)
%   which the host calls after every state change to sync the display (index,
%   total, accept/reject color cue, and a free-form status string).
%
%   The same api callbacks should be driven by the host's keyboard handler so
%   buttons and shortcuts behave identically.
%
%   HISTORY:
%       23 Jun 2026 - created for the ED curation GUI redesign

narginchk(2, 2)
if ~isstruct(api), error('tblgui:eventPanel:api', 'api must be a struct of callbacks'); end

clrA = [0.10 0.55 0.10];
clrR = [0.65 0.15 0.15];

g = uigridlayout(parent, [5, 2], ...
    'RowHeight', {'fit', 'fit', 'fit', 'fit', 'fit'}, ...
    'ColumnWidth', {'1x', '1x'}, ...
    'Padding', 4, 'RowSpacing', 4, 'ColumnSpacing', 4);

% row 1: event index + total
hIdx = uieditfield(g, 'numeric', 'Limits', [1, Inf], 'RoundFractionalValues', 'on', ...
    'Value', 1, 'ValueChangedFcn', @(s, ~) safecall(api, 'setIdx', round(s.Value)));
hIdx.Layout.Row = 1; hIdx.Layout.Column = 1;
hTot = uilabel(g, 'Text', '/ 0', 'VerticalAlignment', 'center');
hTot.Layout.Row = 1; hTot.Layout.Column = 2;

% row 2: prev / next  (arrows match the host keyboard shortcuts)
bPrev = uibutton(g, 'Text', 'Prev (←)', 'ButtonPushedFcn', @(~, ~) safecall(api, 'prev'));
bPrev.Layout.Row = 2; bPrev.Layout.Column = 1;
bNext = uibutton(g, 'Text', 'Next (→)', 'ButtonPushedFcn', @(~, ~) safecall(api, 'next'));
bNext.Layout.Row = 2; bNext.Layout.Column = 2;

% row 3: accept / reject
bAcc = uibutton(g, 'Text', 'Accept (↑)', 'BackgroundColor', clrA, 'FontColor', 'w', ...
    'ButtonPushedFcn', @(~, ~) safecall(api, 'accept'));
bAcc.Layout.Row = 3; bAcc.Layout.Column = 1;
bRej = uibutton(g, 'Text', 'Reject (↓)', 'BackgroundColor', clrR, 'FontColor', 'w', ...
    'ButtonPushedFcn', @(~, ~) safecall(api, 'reject'));
bRej.Layout.Row = 3; bRej.Layout.Column = 2;

% row 4: save (spans both columns)
bSave = uibutton(g, 'Text', 'Save (Ctrl+S)', 'ButtonPushedFcn', @(~, ~) safecall(api, 'save'));
bSave.Layout.Row = 4; bSave.Layout.Column = [1, 2];

% row 5: status (spans both columns)
hStat = uilabel(g, 'Text', '', 'WordWrap', 'on', 'FontColor', [0.3 0.3 0.3]);
hStat.Layout.Row = 5; hStat.Layout.Column = [1, 2];

h = struct('grid', g, 'idx', hIdx, 'total', hTot, ...
    'accept', bAcc, 'reject', bRej, 'save', bSave, 'status', hStat, ...
    'clrAccept', clrA, 'clrReject', clrR);
h.refresh = @(idx, total, accepted, statusText) refreshPanel(h, idx, total, accepted, statusText);

end     % MAIN

% ------------------------------------------------------------------------
function refreshPanel(h, idx, total, accepted, statusText) %#ok<INUSD>
% sync the index field + total. The accept/reject state is shown in the plot
% panels themselves, so no status prefix is added here; statusText is shown as
% given (the host may pass '' for none).
if ~isvalid(h.idx), return; end
total = max(0, round(total));
h.idx.Limits = [1, max(1, total)];
if total > 0
    h.idx.Value = min(max(1, round(idx)), total);
end
h.total.Text = sprintf('/ %d', total);
if nargin < 5 || isempty(statusText), statusText = ''; end
h.status.Text = statusText;
h.status.FontColor = [0.3 0.3 0.3];
end

% ------------------------------------------------------------------------
function safecall(api, name, varargin)
% invoke api.(name) if present and non-empty; otherwise no-op
if isfield(api, name) && ~isempty(api.(name))
    api.(name)(varargin{:});
end
end

% EOF
