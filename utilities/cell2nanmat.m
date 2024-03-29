function mat = cell2nanmat(c, dim)

% converts a cell array of n vectors with maximum length m to a matrix n x
% m. vectors shorter than m are nan padded. can also concatenate an
% array of 2d matrices along the dim dimension by nan padding the other
% dimension. c must be a 1d cell array. 
%
% INPUT:
%   c           cell of numeric vectors
%   dim         dimension to concatenate along
%
% OUTPUT
%   mat         
%
% 09 may 20 LH  updates:
% 13 jun 21         adapted for arrays with n dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(dim)
    dim = 1;
end

c = c(:);

% validate cells have the same number of dimensions.
szc = cellfun(@size, c, 'uni', false);
ndim = unique(cellfun(@length, szc, 'uni', true));
if length(ndim) > 1 
    error('all arrays must have the same number of dimensions')
end
cellvec = all(cellfun(@isvector, c(~isempty(c)), 'uni', true));

if cellvec
    c = cellfun(@(x) x(:), c, 'uni', false);
end

% max length along dimension
maxlength = max(cellfun('size', c, 1));

if cellvec
    mat = cellfun(@(x) [x(:); nan(maxlength - length(x), 1)], c,...
        'UniformOutput', false);
else
    for icell = 1 : length(c)
        if dim == 1
            mat{icell} = [c{icell}; nan(maxlength - size(c{icell}, 1), size(c{icell}, 2))];
        elseif dim == 2
            mat{icell} = [c{icell}, nan(maxlength - size(c{icell}, 2), size(c{icell}, 1))];
        end
    end
end
if dim == 1
    mat = cell2mat(mat);
else
    mat = cell2mat(mat');
end

% EOF


