function datInfo = catDatMemmap(varargin)

% concatenates specific parts of dat files. 
% below also segment to concat states.
%
% INPUT:
%   datFiles    string array with full path and name of dat files. the 
%               order of in which the files are given will determine how
%               they are concatenated
%   parts       cell array where each cell is a n x 2 matrix describing the
%               parts to concatenate [samples]. for example, parts{1} = [80
%               Inf], parts{2} = [0 15], will concatenate the last 80
%               samples of the first file with the first 15 samples of the
%               second file
%   newpath     string. path where new file should be save. if empty than
%               new file will be in pwd
%   newname     string. name of new file. if empty will be extracted from
%               newpath. if newpath not specified will be named as original
%               file with the extension '_new'
%   precision   char. sample precision {'int16'}
%   nchans      numeric. number of channels in dat files {35}.
%   saveVar     logical. save datInfo {true} or not (false).
%
% OUTPUT
%   datInfo     struct with fields describing original and processed files
%
% CALLS:
%   bz_BasenameFromBasepath
%   class2bytes
%
% TO DO LIST:
%   # allow concate tstamps
%   # batch processing (done)
%   # combind with catDat (system commands)
%   # combine parts with clip in n2chunks
%   # allow user to define x time from end in parts
%
% 08 feb 21 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'datFiles', []);
addOptional(p, 'parts', []);
addOptional(p, 'newpath', '', @ischar);
addOptional(p, 'newname', '', @ischar);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
datFiles = p.Results.datFiles;
parts = p.Results.parts;
newpath = p.Results.newpath;
newname = p.Results.newname;
precision = p.Results.precision;
nchans = p.Results.nchans;
saveVar = p.Results.saveVar;

% size of one data point in bytes
nbytes = class2bytes(precision);

% make sure datFiles are cell
datFiles = cellstr(datFiles);

chunksize = 5e5;

% handle names for new path and new file
if isempty(newpath)
    newpath = fileparts(datFiles{1});
    newpath = fullfile(newpath, 'newdat');
    mkdir(newpath)
end
cd(newpath)
basename = bz_BasenameFromBasepath(newpath);
if isempty(newname)
    newname = [basename '.dat'];
end

% handle parts    
if isempty(parts)
    for ifiles = 1 : length(datFiles)
        parts{ifiles} = [1 Inf];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange files and concatenate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open new dat file
fid = fopen(newname, 'w');
if(fid == -1)
    error('cannot open file');
end

for i = 1 : length(datFiles)
    % memory map to original files
    info = dir(datFiles{i});
    nsamps = floor(info.bytes / nbytes / nchans);
    m = memmapfile(datFiles{i}, 'Format', {'int16' [nchans nsamps] 'mapped'});
    raw = m.data;
    
    for ii = 1 : size(parts{i}, 1)
        if parts{i}(ii, 1) < 1
            parts{i}(ii, 1) = 1;
        end
        if parts{i}(ii, 2) > nsamps
            parts{i}(ii, 2) = nsamps;
        end
        chunks = [parts{i}(ii, 1) : chunksize : parts{i}(ii, 2) - chunksize;...
            parts{i}(ii, 1) + chunksize : chunksize : parts{i}(ii, 2)]';
        chunks(end, 2) = parts{i}(ii, 2);
        nchunks = size(chunks, 1);
        for iii = 1 : nchunks
            d = raw.mapped(:, chunks(iii, 1) : chunks(iii, 2));
            
            % write
            fwrite(fid, d(:), precision);
        end
    end  
    clear m
end

% close file and map
rc = fclose(fid);
fclose('all');
if rc == -1
    error('cannot save new file');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange datInfo and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infoname = fullfile(newpath, [basename, '.datInfo.mat']);
if exist(infoname, 'file')
    load(infoname)
end

datInfo.origFile = datFiles;
datInfo.parts = parts;

if saveVar
    save(infoname, 'datInfo', '-v7.3');
end

fprintf('that took %.2f minutes\n', toc / 60)

return

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cat states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
statesfile = 'D:\Data\lh81\lh81_210207_045300\lh81_210207_045300.AccuSleep_states.mat';
file1 = 'D:\Data\lh81\lh81_210206_190000\lh81_210206_190000.AccuSleep_states.mat';
file2 = 'D:\Data\lh81\lh81_210207_065300\lh81_210207_065300.AccuSleep_states.mat';
tcut = [2 4];

% params
nstates = 4;
boutLen = 2.5;

% get labels
ss1 = load(file1, 'ss');
ss1 = ss1.ss;
ss2 = load(file2);
ss2 = ss2.ss;
ss = ss1;
l1 = length(ss1.labels);
l2 = length(ss2.labels);
idx1 = l1 : -1 : l1 - tcut(1) * 60 * 60 / 2.5;
idx2 = 1 : 1 + tcut(2) * 60 * 60 / 2.5;
ss.labels = [ss1.labels(idx1); ss2.labels(idx2)];

% convert labels to state bouts. 
for i = 1 : nstates
    binaryVec = zeros(length(ss.labels), 1);
    binaryVec(ss.labels == i) = 1;
    stateEpisodes = binary2bouts('vec', binaryVec, 'minDur', [], 'maxDur', [],...
        'interDur', [], 'exclude', false); % these are given as indices are equivalent to seconds
    ss.boutTimes{i} = stateEpisodes * boutLen;
end

save(statesfile, 'ss')
