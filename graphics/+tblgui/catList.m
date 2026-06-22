function cats = catList(col)
% TBLGUI.CATLIST  Categories present in a grouping column.
%
%   cats = tblgui.catList(col) coerces a logical / string / categorical
%   column to categorical and returns, as a cellstr, the categories that
%   actually appear in the data (preserving category order).
%
%   Centralizes the coerce-then-list idiom previously repeated in every
%   tblGUI_* checkbox and plot loop.

if islogical(col) || ~iscategorical(col)
    col = categorical(col);
end
cats = categories(col);
cats = cats(ismember(cats, unique(col)));

end     % EOF
