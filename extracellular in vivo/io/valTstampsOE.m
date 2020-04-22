function datInfo = valTstampsOE(varargin)

% fix for bug in open ephys binary format. description of problem: the
% first recording in most (but not all) experiments have more data points
% than timestamps. This is true for continuous but not event streams. The
% excess data points appear on all channels (including auxiliary) and their
% value is zero. Thus, it looks as if the data is zero-padded. The number
% of zeros varies between experiments, from 241 to 139249, and they always
% occur within the first second of the recording (zeroPaddingData.pdf). To
% clarify, as far as I could tell this only occurs after I start
% acquisition (press play) but not between recordings. solution: find
% inconsistencies in timestamps and remove corresponding samples from dat
% file. 
%
% INPUT:
%   basepath    string. path to .dat file (not including dat file itself)
%   fname       string. name of dat file. if empty and more than one dat in
%               path, will be extracted from basepath
%   chunksize   size of data to load at once [samples]{5e6}. 
%               if empty will load entire file (be careful!).
%               for 35 channels in int16, 5e6 samples = 350 MB.
%   bkup        logical. keep original file (true) or not {false}
%   nchans      numeric. original number of channels in dat file {35}.
%   precision   char. sample precision {'int16'} 
%   saveVar     logical. save datInfo {true} or not (false).
%
% OUTPUT
%   datInfo     struct with fields describing original and processed files
%
% CALLS:
%   bz_BasenameFromBasepath
%   class2bytes
%   n2chunks
%
% 22 apr 20 LH      


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'chunksize', 5e6, @isnumeric);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'mapch', [], @isnumeric);
addOptional(p, 'rmvch', [], @isnumeric);
addOptional(p, 'pli', 0, @isnumeric);
addOptional(p, 'bkup', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
fname = p.Results.fname;
chunksize = p.Results.chunksize;
precision = p.Results.precision;
nchans = p.Results.nchans;
bkup = p.Results.bkup;
saveVar = p.Results.saveVar;

% size of one data point in bytes
nbytes = class2bytes(precision);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle dat file
cd(basepath)
datFiles = dir([basepath filesep '**' filesep '*dat']);
if isempty(datFiles)
    error('no .dat files found in %s', basepath)
end
if isempty(fname)
    if length(datFiles) == 1
        fname = datFiles.name;
    else
        fname = [bz_BasenameFromBasepath(basepath) '.dat'];
        if ~contains({datFiles.name}, fname)
            error('please specify which dat file to process')
        end
    end
end
[~, basename, ~] = fileparts(fname);
tempname = [basename '.tmp.dat'];
    
% load timestamps
tFiles = dir([basepath filesep '**' filesep '*timestamps.npy']);
for i = 1 : length(datFiles)
    idx = find(strcmp({tFiles.folder}, datFiles(i).folder));
    if length(idx) > 1
        warning(['more than one timestamps.npy found in %s\n',...
            'skipping tstamps concatination'], basepath)
        break
    elseif length(idx) < 1
        warning(['timestamps.npy not found in %s\n',...
            'skipping tstamps concatination'], basepath)
        break
    else
        tstamps = readNPY(fullfile(tFiles(idx).folder, tFiles(idx).name));
    end
end

% find inconsistencies in timestamps and arrange chunks accordingly
idx = find(diff(tstamps) > 1);
if isempty(idx)
    return
end
clip = zeros(length(idx), 2);
for i = 1 : length(idx)
    clip(i, :) = [tstamps(idx), tstamps(idx + 1)];
end
info = dir(fname);
nsamps = info.bytes / nbytes / nchans;
chunks = n2chunks('n', nsamps, 'chunksize', chunksize, 'clip', clip);
nchunks = size(chunks, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write dat without clip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open temp dat file
fid = fopen(tempname, 'w');
if(fid == -1)
    error('cannot open file');
end

m = memmapfile(fname, 'Format', {'int16' [nchans nsamps] 'mapped'});

% go over chunks
for i = 1 : nchunks
    fprintf('working on chunk %d / %d\n', i, nchunks)
    d = m.data.mapped(:, chunks(i, 1) : chunks(i, 2));
    fwrite(fid, d(:), precision);
end

% close file and map
clear m
rc = fclose(fid);
fclose('all');
if rc == -1
    error('cannot save new file');
end

% rename new and original files, remove if no bkup requested
movefile(fname, [basename '_orig.dat'])
movefile(tempname, fname)
if ~bkup
    delete([basename '_orig.dat'])
end

fprintf('\nfinished processing %s. \nFile size = %.2f MB\n', fname, info.bytes / 1e6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange datInfo and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infoname = fullfile(basepath, [basename, '.datInfo.mat']);
if exist(infoname, 'file')
    load(infoname)
end

datInfo.clip = clip;

if saveVar
    save(infoname, 'datInfo');
end


end

% EOF