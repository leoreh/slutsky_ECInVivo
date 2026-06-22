function [active, allCats] = selectedCats(chk)
% TBLGUI.SELECTEDCATS  Checked category labels from a row of checkboxes.
%
%   [active, allCats] = tblgui.selectedCats(chk) takes an array of checkbox
%   handles and returns ACTIVE (cellstr of the checked labels) and ALLCATS
%   (cellstr of every label, regardless of state). ALLCATS gives callers the
%   full, stable category ordering used for color assignment and broadcasts.
%
%   Works for both uifigure uicheckbox (Text/Value) and legacy uicontrol
%   checkboxes (String/Value).

active = {};
allCats = {};
if isempty(chk), return; end

chk = chk(isgraphics(chk));
if isempty(chk), return; end

labels = arrayfun(@(h) string(localLabel(h)), chk);
vals   = arrayfun(@(h) logical(h.Value), chk);
labels = labels(:)';        % force row
vals   = vals(:)';

allCats = cellstr(labels);
active  = cellstr(labels(vals));

end

function s = localLabel(h)
% uicheckbox exposes Text; legacy uicontrol checkbox exposes String.
if isprop(h, 'Text')
    s = h.Text;
else
    s = h.String;
end
end     % EOF
