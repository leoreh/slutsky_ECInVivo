function txt = gui_inputDialog(parent, prompt, dflt)
% GUI_INPUTDIALOG  Ask for one line of text (DPI-robust inputdlg replacement).
%
%   txt = gui_inputDialog(parent, prompt, dflt) shows a small dialog centred on
%   PARENT's window, with an edit field pre-filled with DFLT, and returns the
%   entered text (trimmed) or '' if cancelled. Like gui_chooseDialog, it renders
%   correctly regardless of display scaling, unlike the legacy inputdlg.
%
%   The window is torn down (delete + drawnow) BEFORE this returns, and is
%   deliberately not WindowStyle 'modal'. An app-modal uifigure blocks every
%   other MATLAB window, so if the caller opens a second dialog (an alert, a
%   confirm) while this one is still being destroyed, the two block each other
%   and the session hangs with no way out. uiwait already stops the caller.

if nargin < 3, dflt = ''; end
txt = '';

d = uifigure('Name', 'Input', 'Position', [100, 100, 380, 200]);
cleaner = onCleanup(@() delete(d(isvalid(d))));   % error path only
fig = ancestor(parent, 'figure');
if ~isempty(fig) && isvalid(fig)
    pp = fig.Position;
    d.Position(1:2) = pp(1:2) + (pp(3:4) - d.Position(3:4)) / 2;
else
    movegui(d, 'center');
end

gl = uigridlayout(d, [3, 2], 'RowHeight', {'1x', 'fit', 'fit'}, ...
    'ColumnWidth', {'1x', '1x'});

lb = uilabel(gl, 'Text', prompt, 'WordWrap', 'on');
lb.Layout.Row = 1; lb.Layout.Column = [1, 2];

ed = uieditfield(gl, 'text', 'Value', char(dflt));
ed.Layout.Row = 2; ed.Layout.Column = [1, 2];

accepted = false;
bOk = uibutton(gl, 'Text', 'OK', 'ButtonPushedFcn', @(~, ~) onOk());
bOk.Layout.Row = 3; bOk.Layout.Column = 1;

bCancel = uibutton(gl, 'Text', 'Cancel', ...
    'ButtonPushedFcn', @(~, ~) uiresume(d));
bCancel.Layout.Row = 3; bCancel.Layout.Column = 2;

d.CloseRequestFcn = @(~, ~) uiresume(d);
uiwait(d);

if accepted && isvalid(d)
    txt = strtrim(ed.Value);
end

% take the window down here, not on the way out of the workspace: the caller
% may open its next dialog immediately, and a half-destroyed window would still
% be on screen and taking focus
delete(d(isvalid(d)));
drawnow;

    function onOk()
        accepted = true;
        uiresume(d);
    end

end     % EOF
