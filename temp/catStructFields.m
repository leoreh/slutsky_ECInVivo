cats = [];
s = frBins;

fldNames = fieldnames(s);
fldData = struct2cell(s);
fldSize = cellfun(@size, fldData, 'uni', false);
fldDim = cellfun(@isvector, fldData, 'uni', true);
fldType = cellfun(@class, fldData, 'uni', false);
x = cellfun(@isnumeric, fldData, 'uni', false);

for ifld = 1 : nflds
    
    switch fldType{ifld}
        case 'struct'
            continue
        otherwise
            x = cat(2, frBins.(fldNames{ifld}));
    end
    
end



cellfun(@vertcat, struct2cell(fr),struct2cell(T),'uni',0),fieldnames(S),1

>> Z = cell2struct(cellfun(@vertcat,struct2cell(S),struct2cell(T),'uni',0),fieldnames(S),1);