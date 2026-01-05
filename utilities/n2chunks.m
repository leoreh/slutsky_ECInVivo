function chunks = n2chunks(varargin)
% N2CHUNKS Divides a number of elements into chunks.
%
%   chunks = N2CHUNKS(...) divides a total number of elements into chunks
%   of a specified size. It handles clipping, overlapping, and boundary
%   adjustment to specific points.
%
%   INPUTS:
%       (Key-Value Pairs)
%       'n'           - (num) Total number of elements.
%       'chunksize'   - (num) Number of elements per chunk {n}.
%       'lastChunk'   - (str) Behavior for the last chunk:
%                       'keep'    : (Default) Retain as is (may be smaller).
%                       'extend'  : Merge with previous chunk.
%                       'exclude' : Remove if smaller than chunksize.
%       'nchunks'     - (num) Desired number of chunks (overrides chunksize).
%       'overlap'     - (num) Overlap between chunks (scalar or [pre post]).
%       'clip'        - (num) Intervals to exclude [start end].
%
%   OUTPUTS:
%       chunks        - (num) Matrix [k x 2] of [start, end] indices.
%
%   DEPENDENCIES:
%       SubtractIntervals (FMAToolbox, if clipping used)
%
%   See also: TIMES2RATE

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'n', [], @isnumeric);
addOptional(p, 'chunksize', [], @isnumeric);
addOptional(p, 'nchunks', [], @isnumeric);
addOptional(p, 'overlap', [0 0], @isnumeric);
addOptional(p, 'clip', [], @isnumeric);
addOptional(p, 'lastChunk', 'keep', @ischar);

parse(p, varargin{:})
n         = p.Results.n;
chunksize = p.Results.chunksize;
nchunks   = p.Results.nchunks;
overlap   = p.Results.overlap;
clip      = p.Results.clip;
lastChunk = p.Results.lastChunk;

% Validate overlap
if isscalar(overlap)
    overlap = [overlap overlap];
elseif numel(overlap) > 2
    error('overlap must be a 2-element vector')
end
if isempty(overlap)
    overlap = [0 0];
end

% Determine chunksize
if ~isempty(nchunks)
    chunksize = ceil(n / nchunks);
end


%% ========================================================================
%  PARTITION
%  ========================================================================

if isempty(chunksize)       % one large chunk
    chunks = [1 n];
else                        
    % Standard chunk generation
    nchunks = ceil(n / chunksize);
    chunks = [1 : chunksize : chunksize * nchunks;...
        chunksize : chunksize : chunksize * nchunks]';
    chunks(nchunks, 2) = n;
end


%% ========================================================================
%  OVERLAP
%  ========================================================================

if ~isempty(chunks)
    chunks(:, 1) = chunks(:, 1) - overlap(1);
    chunks(:, 2) = chunks(:, 2) + overlap(2);
end

if ~isempty(chunks)
    chunks(1, 1) = 1;
end


%% ========================================================================
%  LAST CHUNK LOGIC
%  ========================================================================

lastDur = chunks(end, 2) - chunks(end, 1) + 1;
if lastDur < chunksize
    if strcmp(lastChunk, 'exclude')
        chunks(end, :) = [];
    elseif strcmp(lastChunk, 'extend')
        if size(chunks, 1) > 1
            chunks(end-1, 2) = chunks(end, 2);
            chunks(end, :) = [];
        end
    end
    % 'keep' is default, do nothing
end

if ~strcmp(lastChunk, 'exclude')
    chunks(end, 2) = n;
end


%% ========================================================================
%  CLIPPING
%  ========================================================================

if ~isempty(clip)
    chunks = SubtractIntervals(chunks, clip);
end

% Remove chunks greater than n (can occur if clip includes Inf)
chunks(find(chunks > n) : end, :) = [];

end     % EOF
