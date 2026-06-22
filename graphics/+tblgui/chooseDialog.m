function sel = chooseDialog(parent, prompt, options)
% TBLGUI.CHOOSEDIALOG  Pick one option (DPI-robust listdlg replacement).
%
%   sel = tblgui.chooseDialog(parent, prompt, options) returns the chosen
%   option (char) or '' if cancelled. For up to four options it uses
%   uiconfirm (button choice); for more it shows a small modal uifigure with
%   a dropdown. Both render correctly regardless of display scaling, unlike
%   the legacy listdlg.

options = cellstr(options);
options = options(:)';
fig = ancestor(parent, 'figure');
sel = '';
if isempty(options)
    return;
end

if numel(options) <= 4
    choice = uiconfirm(fig, prompt, 'Assign', ...
        'Options', [options, {'Cancel'}], ...
        'DefaultOption', 1, 'CancelOption', numel(options) + 1);
    if ~strcmp(choice, 'Cancel')
        sel = choice;
    end
    return;
end

% Many options: modal dropdown dialog.
d = uifigure('Name', 'Assign', 'Position', [100, 100, 300, 150]);
d.WindowStyle = 'modal';
cleaner = onCleanup(@() delete(d(isvalid(d))));
movegui(d, 'center');

gl = uigridlayout(d, [3, 2], 'RowHeight', {'fit', 'fit', 'fit'}, ...
    'ColumnWidth', {'1x', '1x'});

lb = uilabel(gl, 'Text', prompt);
lb.Layout.Row = 1; lb.Layout.Column = [1, 2];

dd = uidropdown(gl, 'Items', options);
dd.Layout.Row = 2; dd.Layout.Column = [1, 2];

accepted = false;
bOk = uibutton(gl, 'Text', 'OK', 'ButtonPushedFcn', @(~, ~) onOk());
bOk.Layout.Row = 3; bOk.Layout.Column = 1;

bCancel = uibutton(gl, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) uiresume(d));
bCancel.Layout.Row = 3; bCancel.Layout.Column = 2;

d.CloseRequestFcn = @(~, ~) uiresume(d);
uiwait(d);

if accepted && isvalid(d)
    sel = dd.Value;
end

    function onOk()
        accepted = true;
        uiresume(d);
    end

end     % EOF
