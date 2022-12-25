function [tileMat, tileIdx, tiles] = vec2tileMat(vec, ntiles, include)

% organizes a vector as a matrix (nan-padded) according to
% percentiles (columns). also returns indices as a cell array
%
% INPUT
%   vec         numeric
%   ntiles      scalar. number of percentiles to divide vec {2}
%   include     logical same size as vec indicating which elements to
%               include when calculating the percentiles and organizing
%               tileMat. the indices in tileIdx will be according to the
%               original vec (considering all elements). 
%
% OUTPUT
%   tileMat     mat of n / ntiles (rows) x ntiles (columns) where n is the
%               sum of include. the mat values are those of vec sorted per
%               column
%   tileIdx     cell (1 x ntiles) of indices to vec
%   tiles       vec of tile borders
%
% CALLS
%   cell2nanmat
%
% TO DO LIST
%
% 22 oct 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(ntiles)
    ntiles = 2;
end

include = include(:);

% get percentiles value
tiles = prctile(vec(include), [1 / ntiles : 1 / ntiles : 1 - 1 / ntiles] * 100);
tiles = [0, tiles, Inf];

clear tileMat tileIdx
for itile = 1 : ntiles
    tileIdx{itile} = find((vec >= tiles(itile) & vec < tiles(itile + 1)) &...
        include);
    [tileMat{itile}, sidx] = sort(vec(tileIdx{itile}));
    tileIdx{itile} = tileIdx{itile}(sidx);
end

tileMat = cell2nanmat(tileMat, 2);


end

% EOF