function tdt2dat(basepath, store, blocks, chunksize, mapch, rmvch, clip)

% converts tank (TDT) to dat (neurosuite). Concatenates blocks.
% performs basic preprocessing.
%
% INPUT:
%   basepath    path to recording folder {pwd}.
%   store       stream. typically {'Raw1'} or 'Raw2'
%   blocks      vector. blocks to convert {all}. e.g. [1 2 4 5];
%   chunksize   load data in chunks {60} s. if empty will load entire block.
%   mapch       new order of channels {[]}.
%   rmvch       channels to remove (according to original order) {[]}
%   clip        array of mats indicating times to diregard from recording.
%               each cell corresponds to a block. for example:
%               clip{3} = [0 50; 700 Inf] will remove the first 50 s of
%               Block-3 and the time between 700 s and the end of Block-3
%
% CALLS:
%   TDTbin2mat
%
% TO DO LIST:
%   handle chunks better (e.g. linspace)
%
% 06 dec 18 LH.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

if nargin < 1 || isempty(basepath)
    basepath = pwd;
end
if nargin < 2 || isempty(store)
    store = 'Raw1';
end
if nargin < 3 || isempty(blocks)
    blocks = [];
end
if nargin < 4
    chunksize = 60;
end
if nargin < 5 || isempty(mapch)
    mapch = [];
end
if nargin < 6 || isempty(rmvch)
    rmvch = [];
end
if nargin < 6 || isempty(clip)
    clip{max(blocks)} = [];
else
    if length(clip) < max(blocks)
        clip{max(blocks)} = [];
    elseif length(clip) > max(blocks)
        error('clip array larger then number of blocks')
    end
end

% constants
scalef = 1e6;   % scale factor for int16 conversion [uV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rearrange mapch according to rmvch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(mapch)
    mapch = mapch(~ismember(mapch, rmvch));
    nch = length(mapch);
    if sum(mapch > nch)
        mat = [unique(mapch); (1 : nch)]';
        for j = 1 : nch
            mapch(j) = mat(mat(:, 1) == mapch(j), 2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get tank blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(basepath)
blockfiles = dir('Block*');
blocknames = {blockfiles.name};
fprintf(1, '\nFound %d blocks in %s\n\n', length(blocknames), basepath);

if isempty(blocknames)
    error('no blocks in dir %s.', basepath)
end
if ~isempty(blocks)
    blocknames = blocknames(blocks);
end
nblocks = length(blocknames);

% open dat file
[~, basename] = fileparts(basepath);
newname = [basename '.dat'];
fout = fopen(newname, 'w');
if(fout == -1)
    error('cannot open file');
end

% go over blocks
for i = 1 : nblocks
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract info
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    blockpath = fullfile(basepath, blocknames{i});
    fprintf(1, 'Working on %s\n', blocknames{i});
    
    heads = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store, 'HEADERS', 1);
    nsec(i) = heads.stores.(store).ts(end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % partition into chunks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isempty(chunksize)       % load entire block
        nchunks = 1;
        chunks = [0 0];
    else                        % load block in chunks
        nchunks = ceil(nsec(i) / chunksize);
        chunks = [0 : chunksize : chunksize * (nchunks - 1); chunksize : chunksize : chunksize * nchunks]';
        chunks(nchunks, 2) = 0;
        chunks(1, 1) = 0;
    end
    
    % clip unwanted times
    clipblk = clip{blocks(i)};
    if ~isempty(clipblk)
        if isempty(chunksize)
            if size(clipblk, 1) == 1 && clipblk(1) == 0
                chunks(1, 1) = clipblk(1, 2);
            elseif size(clipblk, 1) == 1 && clipblk(2) == Inf
                chunks(1, 2) = clipblk(1, 1);
            else
                for j = 1 : size(clipblk, 1) - 1
                    chunks = [chunks; clipblk(j, 2) clipblk(j + 1, 1)];
                end
            end
            if clipblk(1, 1) > 0
                chunks = [0 clipblk(1, 1); chunks];
            end
        else
            idx = zeros(1, 2);
            for j = 1 : size(clipblk, 1)
                idx(1) = find(chunks(:, 2) > clipblk(j, 1), 1, 'first');
                chunks(idx(1), 2) = clipblk(j, 1);
                if clipblk(j, 2) ~= Inf
                    idx(2) = find(chunks(:, 1) > clipblk(j, 2), 1, 'first');
                    chunks(idx(2), 1) = clipblk(j, 2);
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % go over chunks, remap, remove and write
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1 : nchunks
        data = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store, 'T1', chunks(j, 1), 'T2', chunks(j, 2));
        data = data.streams.(store).data;
        
        if ~isempty(rmvch)                      % remove channels
            data(rmvch, :) = [];
        end
        
        if ~isempty(mapch)                      % remap channels
            data = data(mapch, :);
        end
        
        fwrite(fout, data(:), 'int16');         % write data
        
        clear data
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rc = fclose(fout);
if rc == -1
    error('cannot save new file');
end

info = dir(newname);
fprintf(1, '\nCreated %s. \nFile size = %.2f MB\n', newname, info.bytes / 1e6);
fprintf(1, '\nElapsed time = %.2f seconds\n', toc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear info
info.filename = basename;
info.store = store;
info.blocks = blocks;
info.blockduration = nsec;
info.clip = clip;
info.rmvch = rmvch;
info.mapch = mapch;
info.fs = heads.stores.(store).fs;

save([basename, '.tdtInfo.mat'], 'info');

end

% EOF