function cats = catfields(s, varargin)

% cats fields of a struct
%
% INPUT:
%   catdef          default behavior of concatenation. can be 'long',
%                   'symmetric', 'addim', 'cell', or a scalar. if the default is
%                   not possible will revert to cell 
%   force           logical. if true, will try to concatenate fields even
%                   if requested default is not available. will even try to
%                   transpose array to allow for concatenation
%
% OUTPUT:
%   cats            struct with concatenated fields
%
% DEPENDENCIES:
%
% TO DO LIST:
%   add option to cat by adding dimension (done)
%   add option to transponse
%   add option to specify dim for specific fields
%   add option to input separate structs and handle unique fields
%   add option to cat by cell2nanmat
%
% 28 feb 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'catdef', 1);
addOptional(p, 'force', true, @islogical);

parse(p, varargin{:})
catdef      = p.Results.catdef;
force       = p.Results.force;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cats

% create 2d cell array of struct x field
flds = fieldnames(s(1));    % make sure this will work for dissimilar structs or add perp step of concat
for istruct = 1 : length(s)
    vals{istruct, :} = struct2cell(s(istruct));       
    for ifld = 1 : length(flds)
        valarray{istruct, ifld} = vals{istruct}{ifld};
        
        % check data type of field
        tf{istruct, ifld} = class(valarray{istruct, ifld});
    end
end

% cat
for ifld = 1 : length(flds)
    
    % if all struct iterate
    if all(strcmp({tf{:, ifld}}, 'struct'))
        cats.(flds{ifld}) = catfields([valarray{:, ifld}], 'catdef', catdef);
        continue
    end
    
    % if any struct cat as cell
    if any(strcmp({tf{:, ifld}}, 'struct'))
        cats.(flds{ifld}) = valarray(:, ifld);
        continue
    end
    
    % if catdef cell cat as cell
    if strcmp(catdef, 'cell')
        cats.(flds{ifld}) = valarray(:, ifld);
        continue
    end

    valsz = cellfun(@size, valarray(:, ifld), 'uni', false);
    ndim = cellfun(@numel, valsz, 'uni', true);
    
    % cat as cell if arrays of same field have different number of
    % dimensions
    if ~all(ndim == ndim(1))
        cats.(flds{ifld}) = valarray(:, ifld);
        continue
    end
    
    % find dimensions that can be used for cat
    valsz = vertcat(valsz{:});
    eqdim = false(1, ndim(1));
    for idim = 1 : ndim
        eqdim(idim) = all(valsz(:, idim) == valsz(1, idim));
    end
    
    % cat as cell if more than one dimension is unequal
    if sum(~eqdim) > 1
        cats.(flds{ifld}) = valarray(:, ifld);
        continue
    end
    
    % cat to the only dimension available
    if sum(~eqdim) == 1
        cats.(flds{ifld}) = cat(find(~eqdim), valarray{:, ifld});
        continue
    end
    
    % cat by adding another dimension if possible
    if strcmp(catdef, 'addim') && sum(eqdim) == unique(ndim)
        cats.(flds{ifld}) = cat(length(eqdim) + 1, valarray{:, ifld});
        cats.(flds{ifld}) = squeeze(cats.(flds{ifld}));
        continue
    end
    
    % sort the dimensions available for concatenation
    [~, sdim] = sort(valsz(1, :));
    sdim = sdim(eqdim);

    % cat according to user selection / default if multiple options exist
    if strcmp(catdef, 'long')
        cats.(flds{ifld}) = cat(sdim(end), valarray{:, ifld});
    elseif strcmp(catdef, 'symmetric')
        cats.(flds{ifld}) = cat(sdim(1), valarray{:, ifld});
    elseif isnumeric(catdef)
        if eqdim(catdef)
            cats.(flds{ifld}) = cat(catdef, valarray{:, ifld});
        else
            cats.(flds{ifld}) = cat(find(~eqdim, 1), valarray{:, ifld});
        end
    end
end

end

% EOF