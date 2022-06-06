function [chunks] = n2chunks(varargin)

% gets a number of elements and divides it to chunks of chunksize. handles
% start/end points. allows clipping of certain elements and overlap between
% chunks. update dec 21: fixed problem that occured when two or more clips
% were in the same chunk. also allows specifying exact points for chunk
% transition
%
% INPUT:
%   n           numeric. number of elements to split
%   chunksize   numeric. number of elements in a chunk {1e6}. 
%   nchunks     numeric. number of chunks to create. will override
%               chunksize.
%   overlap     2-element vector defining the overlap between chunks. 
%               first element is reduced from chunk start and second
%               element is added to chunk end. if single value is specified
%               than the overlap will be symmetrical start to end. for
%               example, overlap = [100 150] for chunksize = 1000
%               chunks may be [1 1150; 900 2150]. 
%   clip        mat n x 2 indicating samples to diregard from chunks.
%               for example: clip = [0 50; 700 Inf] will remove the first
%               50 samples and all samples between 700 and n
%   pnts        numeric. specific points that should serve as the start /
%               end of a chunk. will find the chunk closest to pnts
%               and replace its boundries with the specific time point
%
% OUTPUT
%   chunks      mat n x 2 
%
% EXAMPLE
%   chunks = n2chunks('n', 80000, 'nchunks', 8, 'pnts', 74980, 'overlap', [1000 0])
%
% TO DO LIST:
%
% 22 apr 20 LH  updates:
% 11 aug 20 LH  overlap
% 14 dec 21 LH  fixed clipping
% 12 jan 22 LH  timepoints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'n', [], @isnumeric);
addOptional(p, 'chunksize', [], @isnumeric);
addOptional(p, 'nchunks', [], @isnumeric);
addOptional(p, 'overlap', [0 0], @isnumeric);
addOptional(p, 'clip', [], @isnumeric);
addOptional(p, 'pnts', [], @isnumeric);

parse(p, varargin{:})
n               = p.Results.n;
chunksize       = p.Results.chunksize;
nchunks         = p.Results.nchunks;
overlap         = p.Results.overlap;
clip            = p.Results.clip;
pnts            = p.Results.pnts;

if numel(overlap) == 1
    overlap = [overlap overlap];
elseif numel(overlap) > 2
    error('overlap must be a 2-element vector')
end
if isempty(overlap)
    overlap = [0 0];
end

if ~isempty(nchunks)
    chunksize = ceil(n / nchunks);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partition into chunks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% partition into chunks
if isempty(chunksize)       % one large chunk
    chunks = [1 n];
else                        % load file in chunks
    nchunks = ceil(n / chunksize);
    chunks = [1 : chunksize : chunksize * nchunks;...
        chunksize : chunksize : chunksize * nchunks]';
    chunks(nchunks, 2) = n;
end

% insert overlap
chunks(:, 1) = chunks(:, 1) - overlap(1);
chunks(:, 2) = chunks(:, 2) + overlap(2);
chunks(nchunks, 2) = n;
chunks(1, 1) = 1;

% assimilate clip into chunks
for iclip = 1 : size(clip, 1)
    
    % special care if clip includes first samples
    if clip(iclip, 1) == 1
        rmIdx = chunks(:, 1) <= clip(iclip, 2);
        chunks(rmIdx, :) = [];
        chunks(1, 1) = clip(iclip, 2);
    end
    
    % change chunk end to clip start
    clip_start = find(clip(iclip, 1) < chunks(:, 2), 1, 'first');
    chunk_end = chunks(clip_start, 2);
    chunks(clip_start, 2) = clip(iclip, 1) - 1;
    
    if iclip == size(clip, 1)
        % change chunk start to clip end
        chunks(clip_start + 1, 1) = clip(iclip, 2) - 1;
    else
        % add another chunk from clip end to clip + 1 start
        chunks = [chunks(1 : clip_start, :);...
            clip(iclip, 2) + 1, chunk_end;...           
            chunks(clip_start + 1 : end, :)];
    end
    
    % remove chunks that are after clip. this occurs when clip is greater
    % than chunksize.
    rmidx = find(chunks(:, 1) > chunks(:, 2));
    if ~isempty(rmidx)
        replaceidx = find(chunks(rmidx, 1) < chunks(:, 2), 1);
        chunks(replaceidx, 1) = chunks(rmidx, 1);
        chunks(rmidx : replaceidx - 1, :) = [];
    end
end

% remove chunks that are greater than nsamps. this can occur if clip
% includes Inf
chunks(find(chunks > n) : end, :) = [];

% replace boundries with specific time points
chunks = chunks';
for ipnt = 1 : length(pnts)
    [~, pntIdx] = min(abs(pnts(ipnt) - chunks(:)));
    if mod(pntIdx, 2) == 0
        chunks(pntIdx) = floor(pnts(ipnt));
        if pntIdx < numel(chunks)
            chunks(pntIdx + 1) = ceil(pnts(ipnt)) + 1;
        end
    else
        chunks(pntIdx) = ceil(pnts(ipnt)) + 1;
        if pntIdx > 1
            chunks(pntIdx - 1) = floor(pnts(ipnt));
        end
    end
   
end
chunks = chunks';

end
