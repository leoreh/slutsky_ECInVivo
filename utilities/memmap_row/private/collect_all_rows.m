function in_cell = collect_all_rows(in_cell)
% find in in_cell memmap_row objects, collect all data from them

% convert in place (memory expansive)
for iCell = 1:numel(in_cell)
    if isa(in_cell{iCell},'memmap_row')
        in_cell{iCell} = in_cell{iCell}(:)'; %transpose to keep it as row after collecting
    end
end
end