function cats = catfields(s, catdef, copen)

% concatenates fields of a struct array. if a field consists of different
% data types it will be concatenated as a cell. if array dimensions are not
% compatiple for concatenation along the requested dimension they will be
% padded
%
% INPUT:
%   s               struct array to concatenate fields from
%   catdef          default behavior of concatenation. can be 'cell',
%                   addim', or a scalar. 
%   copen           logical. if true, for fields consisting of cells, will
%                   try to concatenate the contents of the cells rather
%                   then the cells themselves
% 
% OUTPUT:
%   cats            struct with concatenated fields
%
% DEPENDENCIES:
%   cell2padmat     for padding arrays when needed
%
% TO DO LIST:
%   add option to cat by adding dimension (done)
%   add option to cat by cell2nanmat (done - 14 apr 24)
%   add option to specify catdef for specific fields
%   add option to input separate structs and handle unique fields
%
% 28 feb 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    copen = false;
end

% prepare
cats = struct();
flds = fieldnames(s(1));

% iterate over each field
for ifld = 1 : length(flds)

    % extract data for the current field
    fieldData = {s.(flds{ifld})};
    fieldTypes = cellfun(@class, fieldData, 'uni', false);
    
    % determine dimension for concatenation
    if strcmp(catdef, 'addim')
        dim = max(cellfun(@ndims, fieldData)) + 1;
    else
        dim = catdef;
    end

    % if default concatenation is 'cell', or the field consists of different
    % data types, concatenate as cell
    if strcmp(catdef, 'cell') || ~all(strcmp(fieldTypes, fieldTypes{1}))            
        cats.(flds{ifld}) = fieldData;
        continue;
    end
    
    % if fields consist of cells, try to cat across dimension or revert to
    % concatenating as cell. if cellOpen will try to cat the contents of
    % the cells rather then the cells themselves
    if all(strcmp(fieldTypes, 'cell'))
        if copen
            try
                cellContents = {fieldData{:}};
                
                % determine dimension for concatenation
                if strcmp(catdef, 'addim')
                    cellDim = max(cellfun(@ndims, cellContents)) + 1;
                else
                    cellDim = catdef;
                end

                cats.(flds{ifld}) = cell2padmat(cellContents, cellDim);
                continue
            end
        end

        try
            cats.(flds{ifld}) = (cat(dim, fieldData{:}));
        catch
            cats.(flds{ifld}) = fieldData;
        end
        continue
    end

    % if field consists of chars / strings, cat as cell 
    if all(cellfun(@(x) ischar(x) || isstring(x), fieldData))
        try
            cats.(flds{ifld}) = cell2padmat(fieldData, dim)
        catch
            cats.(flds{ifld}) = fieldData;
        end
        continue
    end

    % if the field is made of structs, assure all structs contain the
    % same fields and concatenate them recursively
    if all(cellfun(@isstruct, fieldData))

        % check if all structs have the same fields
        fieldSets = cellfun(@fieldnames, fieldData, 'uni', false);
        refSet = fieldSets{1};
        uniSet = all(cellfun(@(fields) isequal(fields, refSet), fieldSets));
        uniSet = false; % manual override, currently not working (e.g., v(:).psd
        
        if uniSet

            % recursively handle structs
            recurCat = arrayfun(@(x) catfields(x, 'catdef', catdef), [fieldData{:}], 'uni', false);
            cats.(flds{ifld}) = cell2nanmat(recurCat, catdef);
        else
            
            % store as cells if struct fields are not uniform
            cats.(flds{ifld}) = fieldData;
        end

        continue
    end
        
    % if the field is made of numeric / logical arrays concatenate
    % using cell2nanmat
    if all(cellfun(@(x) isnumeric(x) || islogical(x), fieldData))

        cats.(flds{ifld}) = cell2padmat(fieldData, dim);
    
    % if non of the above criteria were met, concatenate as cells
    else
        cats.(flds{ifld}) = fieldData;

    end

end

end

% EOF