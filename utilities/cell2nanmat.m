function [mat, cpad] = cell2nanmat(c, dim)

% converts a cell array of numeric arrays to a matrix.
% numeric arrays are nan-padded along dimensions necessary for concatenation.
% 'c' must be a 1D cell array.
%
% INPUT:
%   c               1D cell of numeric matrices
%   dim             dimension to concatenate along
%
% OUTPUT
%   mat             nan-padded and concatenated matrix
%   cpad            1D cell of matrices after nan padding
%
% 09 may 20 LH  
% 13 jun 21 LH      adapted for arrays with n dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(dim)
    dim = 1;
end

% ensure c is a column array
c = c(:);

% find the maximum size in each dimension across all arrays
maxSize = max(cell2mat(cellfun(@(x) padarray(size(x),...
    [0 max(0, dim - numel(size(x)))], 'post'), c, 'uni', false)), [], 1);

% initialize cell array for padded arrays
cpad = cell(size(c));

% iterate through each cell to pad necessary dimensions
for icell = 1 : length(c)
    
    x = c{icell}; % Current array
    sz = size(x); % Current size

    % ensure 'sz' covers up to 'dim' dimensions
    if numel(sz) < dim
        sz = [sz ones(1, dim - numel(sz))];
    end

    % determine padding size for each dimension
    padSize = maxSize - sz;
    padSize(dim) = 0; % no padding in the concatenation dimension
    
   % adjust padSize to ensure all elements are nonnegative
    padSize = max(padSize, 0);

    % pad array with NaNs in necessary dimensions
    cpad{icell} = padarray(x, padSize, NaN, 'post');
end

% concatenate padded arrays along the specified dimension
mat = cat(dim, cpad{:});

% EOF


