function datInfo = preprocDat(varargin)

% pre-process dat files. copies or concatenates to new destination. if no
% new desitation specified than creates backup and works on the new file.
% removes dc, removes pli, remaps and/or removes channels, handles
% acceleration, creates info file based on xml.
%
% INPUT:
%   datpath     string. path to .dat file (not including dat file itself)
%   fname       string. name of dat file. if empty and more than one dat in
%               path, will be extracted from datpath
%   chunksize   size of data to load at once [samples]{5e6}. 
%               if empty will load entire file (be careful!).
%               for 35 channels in int16, 5e6 samples = 350 MB.
%   bkup        logical. save original file (true) or not {false}
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
%
% TO DO LIST:
%   remove abbarent samples if exist
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
addOptional(p, 'datpath', pwd);
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
datpath = p.Results.datpath;
fname = p.Results.fname;
chunksize = p.Results.chunksize;
precision = p.Results.precision;
nchans = p.Results.nchans;
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
cd(datpath)
datFiles = dir([datpath filesep '**' filesep '*dat']);
if isempty(datFiles)
    error('no .dat files found in %s', datpath)
end
if isempty(fname)
    if length(datFiles) == 1
        fname = datFiles.name;
    else
        fname = [bz_BasenameFromBasepath(datpath) '.dat'];
        if ~contains({datFiles.name}, fname)
            error('please specify which dat file to process')
        end
    end
end
[~, basename, ~] = fileparts(fname);
tempname = [basename '.tmp.dat'];

% rearrange mapch according to rmvch
datInfo.origCh = mapch;
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
if isempty(chunksize)       % load entire file
    nchunks = 1;
    chunks = [1 nsamps * nchans];
else                        % load file in chunks
    nchunks = ceil(nsamps / chunksize);
    chunks = [1 : chunksize * nchans : chunksize * nchunks * nchans;...
        chunksize * nchans : chunksize * nchans : chunksize * nchunks * nchans]';
    chunks(nchunks, 2) = nsamps * nchans;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over chunks and process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open bkup dat file
fid = fopen(tempname, 'w');
if(fid == -1)
    error('cannot open file');
end

% go over chunks
for i = 1 : nchunks
    fprintf('working on chunk %d / %d\n', i, nchunks)
    bsize = (diff(chunks(i, :)) + 1) / nchans;
    m = memmapfile(fname, 'Format', precision);
    d = reshape(m.data(chunks(i, 1) : chunks(i, 2)), [nchans bsize]);
       
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
    
    % acceleration
    
    fwrite(fid, d(:), precision);
end

% close file and map
clear m
rc = fclose(fid);
fclose('all')
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
    error('processing failed, dats are of different length')
end
fprintf('\nfinished processing %s. \nFile size = %.2f MB\n', fname, info.bytes / 1e6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange datInfo and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datInfo.dc = mean(dc, 1);
if saveVar
    save(fullfile(newpath, 'datInfo.mat', 'datInfo'));
end

toc
end
% info = dir(newname);


% EOF