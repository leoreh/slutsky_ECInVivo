function [clrMat, idxOf] = groupColors(fullCats, varargin)
% TBLGUI.GROUPCOLORS  Stable per-category colors.
%
%   [clrMat, idxOf] = tblgui.groupColors(fullCats) returns an N x 3 color
%   matrix for the full (unfiltered) category list FULLCATS, plus a function
%   idxOf(name) that maps a category label to its row in CLRMAT. Looking up
%   colors by name keeps a category's color fixed as other categories are
%   filtered in and out.
%
%   Name-value:
%       'BaseColors'   (M x 3) palette to cycle through (e.g. cfg.clr.unit).
%                      When empty, uses lines(). {[]}
%       'Distinguish'  (logical) use distinguishable_colors instead of
%                      lines() when no BaseColors given. {false}

p = inputParser;
addParameter(p, 'BaseColors', [], @(x) isempty(x) || (isnumeric(x) && size(x, 2) == 3));
addParameter(p, 'Distinguish', false, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
base = p.Results.BaseColors;

fullCats = cellstr(fullCats);
nFull = max(numel(fullCats), 1);

if ~isempty(base)
    clrMat = base(mod(0:nFull-1, size(base, 1)) + 1, :);
elseif p.Results.Distinguish
    clrMat = distinguishable_colors(nFull);
else
    clrMat = lines(nFull);
end

idxOf = @(name) localIdx(fullCats, name);

end

function r = localIdx(fullCats, name)
r = find(strcmp(fullCats, char(name)), 1);
if isempty(r), r = 1; end
end     % EOF
