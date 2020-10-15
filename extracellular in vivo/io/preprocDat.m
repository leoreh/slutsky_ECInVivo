function datInfo = preprocDat(varargin)

% pre-process dat files. copies or concatenates to new destination. if no
% new desitation specified than creates backup and works on the new file.
% removes dc, removes pli, remaps and/or removes channels, handles
% acceleration, removes specified samples, creates info file.
%
% INPUT:
%   basepath    string. path to .dat file (not including dat file itself)
%   fname       string. name of dat file. if empty and more than one dat in
%               path, will be extracted from basepath
%   chunksize   size of data to load at once [samples]{5e6}. 
%               if empty will load entire file (be careful!).
%               for 35 channels in int16, 5e6 samples = 350 MB.
%   bkup        logical. save original file (true) or not {false}
%   clip        mat n x 2 indicating samples to diregard from chunks.
%               for example: clip = [0 50; 700 Inf] will remove the first
%               50 samples and all samples between 700 and n
%   nchans      numeric. original number of channels in dat file {35}.
%   mapch       vec. new order of channels {[]}. 1-based.
%   rmvch       vec. channels to remove (according to original order) {[]}.
%               1-based.
%   pli         numeric. channel from which to extract line crossings.
%               if 0 will not remove pli. 
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
% TO DO LIST:
%   check if works on linux
%   conversion ot mV
%   handle xml
%   datInfo
%
% 09 apr 20 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'chunksize', 5e6, @isnumeric);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'clip', [], @isnumeric);
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
clip = p.Results.clip;
mapch = p.Results.mapch;
rmvch = p.Results.rmvch;
pli = p.Results.pli;
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

% rearrange mapch according to rmvch
orig_mapch = mapch;
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
    
% partition into chunks
info = dir(fname);
nsamps = info.bytes / nbytes / nchans;
chunks = n2chunks('n', nsamps, 'chunksize', chunksize, 'clip', clip);
nchunks = size(chunks, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over chunks and process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\npre-processing %s\n', fname)

% open bkup dat file
fid = fopen(tempname, 'w');
if(fid == -1)
    error('cannot open file');
end

% memory map to original file
m = memmapfile(fname, 'Format', {'int16' [nchans nsamps] 'mapped'});
raw = m.data;

% go over chunks
for i = 1 : nchunks
    % print progress
    if i ~= 1
        fprintf(repmat('\b', 1, length(txt)))
    end
    txt = sprintf('working on chunk %d / %d', i, nchunks);
    fprintf(txt)
    
    % load chunk
    d = raw.mapped(:, chunks(i, 1) : chunks(i, 2));
       
    % remove channels
    if ~isempty(rmvch)                     
        % save channel used for line crossing detection
        if pli
            dpli = d(pli, :);
        end
        d(rmvch, :) = [];
    end
    
    % remap channels
    if ~isempty(mapch)                    
        d = d(mapch, :);
    end
          
    % remove dc
    [d, dc(i, :)] = rmDC(d, 'dim', 2);
    
    % remove pli
    if pli                                 
        linet = lineDetect('x', d(pli, :), 'fs', fs, 'graphics', false);
        for j = 1 : size(d, 1)
            d(j, :) = lineRemove(d(j, :), linet, [], [], 0, 1);
        end
    end
    
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

% check new file
info = dir(fname);
nsampsNew = info.bytes / nbytes / (nchans - length(rmvch));
if ~isequal(nsampsNew, nsamps)
    warning('processing failed, dats are of different length')
end
fprintf('\nfinished processing %s. \nFile size = %.2f MB\n', fname, info.bytes / 1e6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange datInfo and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infoname = fullfile(basepath, [basename, '.datInfo.mat']);
if exist(infoname, 'file')
    load(infoname)
end

datInfo.mapCh = orig_mapch;
datInfo.rmvCh = rmvch;
datInfo.pli = pli;
datInfo.dc = mean(dc, 1);

if saveVar
    save(infoname, 'datInfo', '-v7.3');
end

fprintf('that took %.2f minutes\n', toc / 60)

end

% EOF