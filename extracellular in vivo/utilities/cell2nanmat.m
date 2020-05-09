function mat = cell2nanmat(c)

% converts a cell of n vectors with maximum length m to a matrix n x m.
% vectors shorter than m are nan padded.
%
% INPUT:
%   c           cell of numeric vectors
%
% OUTPUT
%   mat         
%
% 09 may 20 LH      

maxlength = 1;
for i = 1 : length(c)
    maxlength = max([maxlength, length(c{i})]);
end

mat = cellfun(@(x)[x(:); NaN(maxlength-length(x), 1)], c,...
    'UniformOutput', false);
mat = cell2mat(mat);
