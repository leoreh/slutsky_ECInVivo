function [mat, cpad] = cell2nanmat(c, dim, paddef)

% converts a cell array to a matrix, paddig as needed.
%
% INPUT:
%   c               1D cell of numeric matrices
%   dim             dimension to concatenate along
%   paddef          defines the value used for padding (nan for numeric,
%                   false for logical, space for char / string, empty cell
%                   for cell array)
%
% OUTPUT
%   mat             nan-padded and concatenated matrix
%   cpad            1D cell of matrices after padding
%
% 09 may 20 LH  
% 13 jun 21 LH      adapted for arrays with n dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get data type from each cell
cellTypes = cellfun(@class, c, 'uni', false);
if ~all(strcmp(cellTypes, cellTypes{1}))
    error('Mixed data types or unsupported data type');
end

% set default padding if not specified
if nargin < 3
    if all(cellfun(@isnumeric, c))
        paddef = NaN;
    elseif all(cellfun(@islogical, c))
        paddef = false;
    elseif all(cellfun(@ischar, c))
        paddef = ' '; 
    elseif all(cellfun(@(x) isa(x, 'string'), c))
        paddef = "";  
    elseif all(cellfun(@(x) iscell(x), c))
        paddef = cell(1, 1);
    end
end

% get maximum size across each dimension
maxDim = max(cellfun(@ndims, c, 'uni', true));
maxSize = zeros(1, maxDim);
for icell = 1 : length(c)
    sz = size(c{icell});
    sz = [sz, ones(1, maxDim - numel(sz))];  
    maxSize = max(maxSize, sz);
end

% initialize cell array for padded arrays
cpad = cell(size(c));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pad necessary dimensions for each cell
for icell = 1 : length(c)
    
    x = c{icell}; 
    sz = [size(x), ones(1, maxDim - ndims(x))];

    % determine padding size for each dimension
    padSize = maxSize - sz;
    padSize(dim) = 0; % no padding in the concatenation dimension
    
    % adjust padSize to ensure all elements are nonnegative
    padSize = max(padSize, 0);

    if isnumeric(x) || islogical(x) 
        cpad{icell} = padarray(x, padSize, paddef, 'post');
    
    elseif iscell(x)
        cpad{icell} = repmat(paddef, maxSize);
        idx = arrayfun(@(x) 1 : x, sz, 'uni', false);
        cpad{icell}(idx{:}) = x;

    elseif ischar(x) || isa(x, 'string')
        
        % manually pad character and string arrays, use original if no
        % padding is needed
        if dim <= numel(sz)
            
            if ischar(x)
                cpad{icell} = [x repmat(paddef, size(x, 1), max(padSize(2 : end), 0))];
            elseif isa(x, 'string')
                cpad{icell} = [x; repmat(paddef, max(padSize(2 : end), 0), 1)];
            end

        else
            cpad{icell} = x;  
        end
    end
end

% concatenate padded arrays along the specified dimension
mat = cat(dim, cpad{:});

end
 
% EOF
