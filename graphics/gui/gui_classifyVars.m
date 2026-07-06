function [numVars, catVars] = gui_classifyVars(tbl, varargin)
% GUI_CLASSIFYVARS  Split table variable names into numeric and categorical.
%
%   [numVars, catVars] = gui_classifyVars(tbl) returns cellstr lists of
%   the numeric variable names and the categorical-like (categorical /
%   string / logical) variable names of TBL.
%
%   Name-value pairs:
%       'Exclude'    (cellstr) numeric names to drop (e.g. id columns) {{}}
%       'SortNames'  (logical) alphabetically sort both lists {false}
%       'Shape'      'any' (default) returns every numeric column; 'vector'
%                    returns only scalar-per-row numerics (one column wide).
%                    Bar / scatter aggregate one value per row, so they pass
%                    'vector' to exclude matrix columns (e.g. traces, ACGs).
%
%   This is the single source for the numeric-vs-categorical split that was
%   previously duplicated across the guiTbl_* family.

p = inputParser;
addParameter(p, 'Exclude', {}, @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(p, 'SortNames', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Shape', 'any', @(x) any(strcmpi(x, {'any', 'vector'})));
parse(p, varargin{:});
exclude = cellstr(p.Results.Exclude);

allVars = tbl.Properties.VariableNames;
if strcmpi(p.Results.Shape, 'vector')
    isNum = varfun(@(x) isnumeric(x) && size(x, 2) == 1, tbl, 'OutputFormat', 'uniform');
else
    isNum = varfun(@isnumeric, tbl, 'OutputFormat', 'uniform');
end
numVars = allVars(isNum);
catVars = allVars(varfun(@(x) iscategorical(x) || isstring(x) || islogical(x), ...
    tbl, 'OutputFormat', 'uniform'));

if ~isempty(exclude)
    numVars = setdiff(numVars, exclude, 'stable');
end

if p.Results.SortNames
    numVars = sort(numVars);
    catVars = sort(catVars);
end

end     % EOF
