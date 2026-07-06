function gui_notify(parent, msg, kind)
% GUI_NOTIFY  Alert dialog for uifigure GUIs (replaces msgbox / errordlg).
%
%   gui_notify(parent, msg) shows an information alert on the uifigure that
%   owns PARENT. kind (optional): 'info' (default) | 'warning' | 'error' |
%   'success'.

if nargin < 3 || isempty(kind), kind = 'info'; end
fig = ancestor(parent, 'figure');

switch lower(kind)
    case {'error', 'err'},      icon = 'error';   ttl = 'Error';
    case {'warning', 'warn'},   icon = 'warning'; ttl = 'Warning';
    case 'success',             icon = 'success'; ttl = 'Success';
    otherwise,                  icon = 'info';    ttl = 'Info';
end

% uialert requires a visible uifigure. Fall back to the command window when
% the figure is hidden (e.g. headless smoke tests) so callers never crash.
if ~isempty(fig) && isprop(fig, 'Visible') && strcmp(fig.Visible, 'on')
    uialert(fig, msg, ttl, 'Icon', icon);
else
    fprintf('[%s] %s\n', ttl, msg);
end

end     % EOF
