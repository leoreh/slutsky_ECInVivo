function [datInfo, data] = tdt2dat(varargin)

% converts tank (TDT) to dat (neurosuite). Concatenates blocks.
% performs basic preprocessing. if data (second output) is requested than
% will output a matrix of data instead of a .dat file
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
%               clip{3} = [0 50; 700 Inf] will remove the first 50 s from
%               block-3 and the time between 700 s and the end of Block-3
%   saveVar     logical. save variable {true} or not (false).
%
% OUTPUT
%   info        struct with fields describing tdt params
%   data        mat (channels x samples). only if requested
%
% CALLS:
%   TDTbin2mat
%   natsort
%
% TO DO LIST:
%   handle chunks better (e.g. linspace)
%   add time limit to split files
%   separate chunks to standalone function (see n2chunks, not implemented)
%   include clip within blockduration
%
% 06 dec 18 LH      updates:
% 18 sep 19 LH      handle arguments
% 27 mar 20 LH      allow output mat instead of dat file
% 04 sep 20 LH      added nsamps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'store', 'Raw1', @ischar);
addOptional(p, 'blocks', [], @isnumeric);
addOptional(p, 'chunksize', 60, @isnumeric);
addOptional(p, 'mapch', [], @isnumeric);
addOptional(p, 'rmvch', [], @isnumeric);
addOptional(p, 'clip', {}, @iscell);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
store = p.Results.store;
blocks = p.Results.blocks;
chunksize = p.Results.chunksize;
mapch = p.Results.mapch;
rmvch = p.Results.rmvch;
clip = p.Results.clip;
saveVar = p.Results.saveVar;

if isempty(clip)
    clip{max(blocks)} = [];
else
    if length(clip) < max(blocks)
        clip{max(blocks)} = [];
    elseif length(clip) > max(blocks)
        error('clip array larger then number of blocks')
    end
end

[~, basename] = fileparts(basepath);

% constants
scalef = 1e6;   % scale factor for int16 conversion [uV]

% initialize
data = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rearrange mapch according to rmvch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(mapch)
    mapch = mapch(~ismember(mapch, rmvch));
    nch = length(mapch);
    if sum(mapch > nch)
        mat = [unique(mapch); (1 : nch)]';
        for ich = 1 : nch
            mapch(ich) = mat(mat(:, 1) == mapch(ich), 2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get tank blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(basepath)
blockfiles = dir('block*');
blocknames = {blockfiles.name};
blocknames = natsort(blocknames);
fprintf(1, '\nFound %d blocks in %s\n\n', length(blocknames), basepath);

if isempty(blocknames)
    error('no blocks in dir %s.', basepath)
end
if ~isempty(blocks)
    blocknames = blocknames(blocks);
end
nblocks = length(blocknames);

% get start/stop date time of blocks
for iblock = 1 : nblocks
    blockpath = fullfile(basepath, ['block-', num2str(blocks(iblock))]);
    blockHead = TDTbin2mat(blockpath, 'HEADERS', 1);
    if ~isnan(blockHead.startTime)
        blockStart{iblock} = datestr(datenum([1970, 1, 1, 0, 0, blockHead.startTime]) + hours(2), 'yymmdd_HHMMss');
    end
    if ~isnan(blockHead.stopTime)
        blockEnd{iblock} = datestr(datenum([1970, 1, 1, 0, 0, blockHead.stopTime]) + hours(2), 'yymmdd_HHMMss');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open dat file
if nargout == 1
    newname = [basename '.dat'];
    fout = fopen(newname, 'w');
    if(fout == -1)
        error('cannot open file');
    end
end

% go over blocks
nsec = zeros(1, nblocks);
nsamps = zeros(1, nblocks);
for i = 1 : nblocks
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract info
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    blockpath = fullfile(basepath, blocknames{i});
    fprintf(1, 'Working on %s\n', blocknames{i});
    
    blockHead = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store, 'HEADERS', 1);
    nsec(i) = blockHead.stores.(store).ts(end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % partition into chunks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    if isempty(chunksize)       % load entire block
        nchunks = 1;
        chunks = [0 0];
    else                        % load block in chunks
        nchunks = ceil(nsec(i) / chunksize);
        chunks = [0 : chunksize : chunksize * (nchunks - 1);...
            chunksize : chunksize : chunksize * nchunks]';
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
        nchunks = size(chunks, 1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % go over chunks, remap, remove and write
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1 : nchunks
        dat = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store,...
            'T1', chunks(j , 1) , 'T2' ,chunks(j , 2));
        dat = dat.streams.(store).data;
        
        if ~isempty(dat)
            if ~isempty(rmvch)                     % remove channels
                dat(rmvch, :) = [];
            end
            
            if ~isempty(mapch)                     % remap channels
                dat = dat(mapch, :);
            end
            
            if nargout == 1  
                fwrite(fout, dat(:), 'int16');     % write data
            elseif nargout == 2
                data = [data, dat];                % output data
            end
        end
        nsamps(i) = [nsamps(i) + size(dat, 2)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 1
    rc = fclose(fout);
    if rc == -1
        error('cannot save new file');
    end
datInfo = dir(newname);
fprintf(1, '\nCreated %s. \nFile size = %.2f MB\n', newname, datInfo.bytes / 1e6);
end
fprintf(1, '\nElapsed time = %.2f seconds\n', toc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear info
datInfo.filename = basename;
datInfo.store = store;
datInfo.blocks = blocknames;
datInfo.blockduration = nsec;
datInfo.clip = clip;
datInfo.rmvch = rmvch;
datInfo.mapch = mapch;
datInfo.fs = blockHead.stores.(store).fs;
datInfo.nsamps = nsamps; 
datInfo.nsec = nsec; 
datInfo.blockStart = blockStart;
datInfo.blockEnd = blockEnd;

if saveVar
    save([basename, '.' store, '.datInfo.mat'], 'datInfo');
    if nargout == 2
        save([basename, '.' store, '.data.mat'], 'data');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% more processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create lfp file
if contains(store, 'Raw') && nargout == 1 
    LFPfromDat('basepath', basepath, 'cf', [450], 'chunksize', 5e6,...
        'nchans', size(dat, 1), 'fsOut', 1250,...
        'fsIn', blockHead.stores.(store).fs)
end

end

% EOF