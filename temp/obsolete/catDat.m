function datInfo = catDat(varargin)

% concatenates all .dat files in basepath via system commands (copy for
% windows). can also be used to copy a dat file from basepath to newpath.
% also concatenates timestamps.npy files (open ephys; assuming there is one
% file in each dat folder) and saves them within datInfo.
% for each file also validates that the length of tstamps and dat is equal,
% if not than removes samples from dat. this is based on a bug found in OE
% binary format.
%
% INPUT:
%   basepath    cell array of strings describing paths to .dat files {pwd}.
%               if multiple dat files exist than concat must be true
%   newpath     string. path where new file should be save. if empty than
%               new file will be save in basepath
%   newname     string. name of new file. if empty will be extracted from
%               newpath. if newpath not specified will be named as original
%               file with the extension '_new'
%   concat      logical. concatenate dat files {true} or not (false)
%   precision   char. sample precision {'int16'}
%   nchans      numeric. number of channels in dat file {35}.
%   saveVar     logical. save datInfo {true} or not (false).
%
% OUTPUT
%   datInfo     struct with fields describing original and processed files
%
% CALLS:
%   bz_BasenameFromBasepath
%   class2bytes
%   datInfo
%   valTstamps
%
% TO DO LIST:
%   # handle multiple dats with copy only (no concatenate)
%   # adapt for linux
%   # change tstamps to binary file (supposed to be much faster read/write)
%   # handle edge effects during concatination
%
% 09 apr 20 LH
% 22 apr 20 LH      valTstampsOE
% 20 may 20 LH      input multiple paths
% 19 aug 20 LH      map to tstamps instead of direct readNPY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'newpath', '', @ischar);
addOptional(p, 'newname', '', @ischar);
addOptional(p, 'concat', true, @islogical);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
newpath = p.Results.newpath;
newname = p.Results.newname;
concat = p.Results.concat;
precision = p.Results.precision;
nchans = p.Results.nchans;
saveVar = p.Results.saveVar;

% size of one data point in bytes
nbytes = class2bytes(precision);

% initialize
datFiles = [];
tFiles = [];
tstamps = [];

% set basepath to cell
basepath = cellstr(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange files and concatenate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle names for new path and new file
if isempty(newpath)
    newpath = basepath{1};
end
if isempty(newname)
    if ~isempty(newpath)
        basename = bz_BasenameFromBasepath(newpath);
        newname = [basename '.dat'];
    else
        [~, basename] = fileparts(datFiles(1).name);
        newname = [basename '.new.dat'];
    end
else
    if ~contains(newname, '.dat')
        newname = [newname '.dat'];
    end
end
destination = fullfile(newpath, newname);

% get .dat and .npy files in basepath
for i = 1 : length(basepath)
    datFiles = [datFiles; dir([basepath{i} filesep '**' filesep '*dat'])];
    tFiles = [tFiles; dir([basepath{i} filesep 'continuous' filesep '**' filesep 'timestamps.npy'])];
end
if length(datFiles) > 1 && ~concat
    error(['multiple dat files found in %s\n',...
        'function cannot handle multiple files without concatination'], basepath)
end

% check .dat files integrity and concat tstamps
for i = 1 : length(datFiles)
    source{i} = fullfile(datFiles(i).folder, datFiles(i).name);
    nsamps(i) = datFiles(i).bytes / nbytes / nchans;
    if ~isequal(nsamps(i), round(nsamps(i)))
        error('incorrect nCh for file')
    end
    
    % load corresponding tstamps file
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
        
        % memory map to tstamps file
        tname = fullfile(tFiles(idx).folder, tFiles(idx).name);
        [arrayShape, dataType, ~, ~, totalHeaderLength, ~] =...
            readNPYheader(tname);
        tmap = memmapfile(tname, 'Format', {dataType, arrayShape(end:-1:1), 'd'},...
            'Offset', totalHeaderLength); 
                
        % fix bug in OE binary format
        if length(tmap.Data.d) < nsamps(i)
            warning(['more samples than timestamps in %s\n'...
                'initializing valTstampsOE'], source{i}) 
            validate_OE_tstamps('basepath', tFiles(idx).folder, 'precision', precision,...
                'chunksize', 5e6, 'bkup', false, 'saveVar', true,...
                'nchans', nchans)
            nsamps(i) = length(tmap.Data.d);
        end       
        tstamps = [tstamps; tmap.Data.d];
        
    end
end

% do the copy / concat
fprintf('\ncreating %s from files\n', newname)
cmd = ['!copy /b ' strjoin(source, ' + '), ' ' destination];
eval(cmd);

% check integrity of new file
info = dir(destination);
nsampsNew = info.bytes / nbytes / nchans;
if ~isequal(nsampsNew, sum(nsamps))
    warning('copying failed, dats are of different length')
end
fprintf('new file size = %.2f MB\n', info.bytes / 1e6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange datInfo and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infoname = fullfile(newpath, [basename, '.datInfo.mat']);
if exist(infoname, 'file')
    load(infoname)
end

datInfo.origFile = source;
datInfo.newFile = destination;
datInfo.nsamps = nsamps;

if saveVar
    save(infoname, 'datInfo', '-v7.3');
    save(fullfile(newpath, [basename, '.tstamps.mat']), 'tstamps', '-v7.3')
end

fprintf('that took %.2f minutes\n', toc / 60)

end

% EOF