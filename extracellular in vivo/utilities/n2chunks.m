function [chunks] = n2chunks(varargin)

% gets a number of elements and divides it to chunks of chunksize. handles
% start/end points. allows clipping of certain elements and overlap between
% chunks.
%
% INPUT:
%   n           numeric. number of elements to split
%   chunksize   numeric. number of elements in a chunk {1e6}. 
%   overlap     2-element vector defining the overlap between chunks. 
%               first element is reduced from chunk start and second
%               element is added to chunk end. if single value is specified
%               than the overlap will be symmetrical start to end. for
%               example, overlap = [100 150] than for chunksize = 1000
%               chunks may be [1 1150; 900 2150]. 
%   clip        mat n x 2 indicating samples to diregard from chunks.
%               for example: clip = [0 50; 700 Inf] will remove the first
%               50 samples and all samples between 700 and n
%
% OUTPUT
%   chunks      mat n x 2 
%
% CALLS:
%
% TO DO LIST:
%   # add an option to restrict minimum chunk size
%
% 22 apr 20 LH  updates:
% 11 aug 20 LH  overlap


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'n', [], @isnumeric);
addOptional(p, 'chunksize', [], @isnumeric);
addOptional(p, 'overlap', [0 0], @isnumeric);
addOptional(p, 'clip', [], @isnumeric);

parse(p, varargin{:})
n = p.Results.n;
chunksize = p.Results.chunksize;
overlap = p.Results.overlap;
clip = p.Results.clip;

if numel(overlap) == 1
    overlap = [overlap overlap];
elseif numel(overlap) > 2
    error('overlap must be a 2-element vector')
end

% validate
% if max(clip(:)) > n
%     error('')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partition into chunks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% partition into chunks
if isempty(chunksize)       % load entire file
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
    if isempty(chunksize)
        if size(clip, 1) == 1 && clip(1) == 0
            chunks(1, 1) = clip(1, 2);
        elseif size(clip, 1) == 1 && clip(2) == Inf
            chunks(1, 2) = clip(1, 1);
        else
            for j = 1 : size(clip, 1) - 1
                chunks = [chunks; clip(j, 2) clip(j + 1, 1)];
            end
        end
        if clip(1, 1) > 0
            chunks = [0 clip(1, 1); chunks];
        end
    else
        idx = zeros(1, 2);
        for j = 1 : size(clip, 1)
            idx(1) = find(chunks(:, 2) > clip(j, 1), 1, 'first');
            chunks(idx(1), 2) = clip(j, 1);
            if clip(j, 2) ~= Inf
                idx(2) = find(chunks(:, 1) > clip(j, 2), 1, 'first');
                chunks(idx(2), 1) = clip(j, 2);
                rmblk = idx(1) + 1 : idx(2) - 1;
            else
                rmblk = idx(1) + 1 : size(chunks, 1);
            end
            if ~isempty(rmblk)
                chunks(rmblk, :) = [];
            end
        end
    end
    chunks = chunks(any(chunks, 2), :);
end