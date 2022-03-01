function cats = catfields(s, varargin)

% cats fields of a struct
%
% INPUT:
%   catdef          default behavior of concatenation. can be 'long',
%                   'symmetric', or a scalar. if the default is not
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
%   add option to cat by adding dimension
%   add option to transponse
%   add option to specify dim for specific fields
%   add option to input separate structs and handle unique fields
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
fields = fieldnames(s(1));    % make sure this will work for dissimilar structs or add perp step of concat
for istruct = 1 : length(s)
    vals{istruct, :} = struct2cell(s(istruct));       
    for ifield = 1 : length(fields)
        valarray{istruct, ifield} = vals{istruct}{ifield};
        
        % check data type of field
        tf{istruct, ifield} = class(valarray{istruct, ifield});
    end
end

% cat
for ifield = 1 : length(fields)
    
    % if all struct iterate
    if all(strcmp({tf{:, ifield}}, 'struct'))
        cats.(fields{ifield}) = catfields([valarray{:, ifield}], 'catdef', catdef);
        continue
    end
    
    % if any struct cat as cell
    if any(strcmp({tf{:, ifield}}, 'struct'))
        cats.(fields{ifield}) = valarray(:, ifield);
        continue
    end
       
    valsz = cellfun(@size, valarray(:, ifield), 'uni', false);
    ndim = cellfun(@numel, valsz, 'uni', true);
    
    % cat as cell if arrays of same field have different number of
    % dimensions
    if ~all(ndim == ndim(1))
        cats.(fields{ifield}) = valarray(:, ifield);
        continue
    end
    
    % find dimensions that can be used for cat
    valsz = vertcat(valsz{:});
    eqdim = false(1, ndim(1));
    for idim = 1 : ndim
        eqdim(idim) = all(valsz(:, idim) == valsz(1, idim));
    end
    
    % cat as cell if no dimensions are equal
    if ~any(eqdim)
        cats.(fields{ifield}) = valarray(:, ifield);
        continue
    end
    
    % cat to the only dimension available
    if sum(eqdim) == 1
        cats.(fields{ifield}) = cat(find(~eqdim), valarray{:, ifield});
        continue
    end
    
    % cat according to user selection / default if multiple options exist
    [~, sdim] = sort(valsz(1, eqdim));
    switch catdef
        case 'long'
            cats.(fields{ifield}) = cat(sdim(end), valarray{:, ifield});
        case 'symmetric'
            cats.(fields{ifield}) = cat(sdim(1), valarray{:, ifield});
        otherwise
            if eqdim(catdef)
                cats.(fields{ifield}) = cat(catdef, valarray{:, ifield});
            else
                cats.(fields{ifield}) = cat(find(~eqdim, 1), valarray{:, ifield});
            end
    end
end
cats

end

% EOF