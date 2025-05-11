function [chunks] = n2chunks(varargin)

% n2chunks Divides a number of elements into chunks.
%
% This function takes a total number of elements and divides it into
% chunks of a specified size. It handles start and end points, allows for
% clipping of certain elements (i.e., excluding them from any chunk), 
% and supports overlap between chunks. It can also adjust chunk boundaries 
% to align with specified points.
% 
% INPUT (name-value pairs):
%   'n'           - Numeric. Total number of elements to split.
%   'chunksize'   - Numeric. Number of elements in each chunk. 
%                   Defaults to n if not specified (i.e., one chunk).
%   'nchunks'     - Numeric. Desired number of chunks. If provided, this
%                   will override 'chunksize'.
%   'overlap'     - Numeric scalar or 2-element vector. Defines the overlap
%                   between chunks. 
%                   If scalar: symmetrical overlap (e.g., overlap = 100 means
%                   chunk [1-1000] becomes [1-1100] and next chunk [1001-2000]
%                   becomes [901-2100], before boundary adjustments).
%                   If 2-element vector [pre post]: 'pre' elements are added
%                   to the start of the chunk (extending into the previous one)
%                   and 'post' elements are added to the end (extending into
%                   the next one). Example: overlap = [100 150] for
%                   chunksize = 1000 could result in chunks like 
%                   [1 1150], [901 2150], etc. (actual start/end may vary due
%                   to boundary conditions like total 'n').
%                   Default: [0 0] (no overlap).
%   'clip'        - Numeric matrix (m x 2). Specifies intervals to exclude
%                   from the chunks. Each row [start end] defines an interval
%                   to clip. 'Inf' can be used for end points.
%                   Example: clip = [0 50; 700 Inf] will ensure that elements
%                   1-50 and 700-n are not included in any chunk output.
%   'pnts'        - Numeric vector. Specific element indices that should ideally
%                   mark a transition between chunks. The function will adjust
%                   the nearest chunk boundary to align with these points.
%
% OUTPUT:
%   chunks      - Numeric matrix (k x 2). Each row defines a chunk as 
%                 [startIndex endIndex]. k is the resulting number of chunks.
% 
% DEPENDENCIES:
%   SubtractIntervals (from FMAToolbox, if 'clip' is used)
%
% EXAMPLE:
%   chunks = n2chunks('n', 80000, 'nchunks', 8, 'pnts', 74980, 'overlap', [1000 0])
%   % This would divide 80000 elements into 8 chunks, try to make a chunk
%   % boundary near 74980, and make each chunk overlap with the previous
%   % one by 1000 elements (except the first chunk).
%
% See also: class2bytes, loadBinary

% REVISIONS:
% 22 Apr 20 LH  Original version.
% 11 Aug 20 LH  Added 'overlap' functionality.
% 14 Dec 21 LH  Fixed clipping issues, especially with multiple clips in one chunk.
% 12 Jan 22 LH  Added 'pnts' functionality to align chunks with specific timepoints.
% 08 Sep 22 LH  Noted dependency on SubtractIntervals for 'clip'.
% [Current Date] [Your Initials] Updated comments for clarity and consistency.

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

if isscalar(overlap)
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
if ~isempty(clip)
    chunks = SubtractIntervals(chunks, clip);
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
