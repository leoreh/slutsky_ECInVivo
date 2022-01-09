function [chunks] = n2nchunks(varargin)

% INPUT:
%   n           numeric. number of elements to split
%   nchunks     numeric. number of chunks to create
%   timepoints  numeric. instead of even chunksize for all chunks, replace
%               the closest index with speicific values. 
%
% OUTPUT
%   chunks      mat n x 2 
%
% 08 jan 22 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'n', [], @isnumeric);
addOptional(p, 'nchunks', 4, @isnumeric);
addOptional(p, 'timepoints', [], @isnumeric);

parse(p, varargin{:})
n = p.Results.n;
nchunks = p.Results.nchunks;
timepoints = p.Results.timepoints;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the chunks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chunks = n2chunks('n', n, 'chunksize', ceil(n / nchunks));
chunks = chunks';

for ipnt = 1 : length(timepoints)
    [~, pntIdx] = min(abs(timepoints(ipnt) - chunks(:)));
    if mod(pntIdx, 2) == 0
        chunks(pntIdx) = floor(timepoints(ipnt));
        if pntIdx < numel(chunks)
            chunks(pntIdx + 1) = ceil(timepoints(ipnt)) + 1;
        end
    else
        chunks(pntIdx) = ceil(timepoints(ipnt)) + 1;
        if pntIdx > 1
            chunks(pntIdx - 1) = floor(timepoints(ipnt));
        end
    end
    
end
chunks = chunks';

end






