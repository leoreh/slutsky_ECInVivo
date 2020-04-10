function datInfo = catDat(varargin)

% concatenates all .dat files in basepath. can also be used to copy a
% datfile from basepath to newpath.
%
% INPUT:
%   datpath     string. path to .dat files {pwd}. 
%               if multiple dat files exist than concat must be true
%   newpath     string. path where new file should be save. if empty than
%               new file will be save in datpath
%   newname     string. name of new file. if empty will be extracted from
%               newpath. in newpath not specified will be named as original
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
%
% TO DO LIST:
%   handle multiple dats with copy only (no concatenate)
%   remove abbarent samples if exist
%   check if works on linux
%   conversion ot mV
%
% 09 apr 20 LH      


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'datpath', pwd);
addOptional(p, 'newpath', '', @ischar);
addOptional(p, 'newname', '', @ischar);
addOptional(p, 'concat', true, @islogical);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
datpath = p.Results.datpath;
newpath = p.Results.newpath;
newname = p.Results.newname;
concat = p.Results.concat;
precision = p.Results.precision;
nchans = p.Results.nchans;
saveVar = p.Results.saveVar;

% size of one data point in bytes
nbytes = class2bytes(precision);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange files and concatenate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get .dat files in datpath and check their integrity
datFiles = dir([datpath filesep '**' filesep '*dat']);
if length(datFiles) > 1 && ~concat
    error(['multiple dat files found in %s\n',...
        'function cannot handle multiple files without concatination'], datPath)
end
for i = 1 : length(datFiles)
    source{i} = fullfile(datFiles(i).folder, datFiles(i).name);
    nsamps(i) = datFiles(i).bytes / nbytes / nchans;
    if ~isequal(nsamps(i), round(nsamps(i)))
        error('incorrect nCh for file')
    end
end

% handle names for new path and new file
if isempty(newpath)
    newpath = datpath;
end
if isempty(newname)
    if ~isempty(newpath)
        newname = [bz_BasenameFromBasepath(newpath) '.dat'];
    else
        [~, newname] = fileparts(datFiles(1).name);
        newname = [newname 'new.dat'];
    end
else
    if ~contains(newname, '.dat')
        newname = [newname '.dat'];
    end
end
destination = fullfile(newpath, newname);
 
% do the copy / concat
fprintf('creating %s from files\n', newname)
cmd = ['!copy /b ' strjoin(source, ' + '), ' ' destination];
eval(cmd); 

% check integrity of new file
info = dir(destination);
nsampsNew = info.bytes / nbytes / nchans;
if ~isequal(nsampsNew, sum(nsamps))
    error('copying failed, dats are of different length')
end
fprintf('\ncreated %s. \nFile size = %.2f MB\n', newname, info.bytes / 1e6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange datInfo and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar
    save(fullfile(newpath, 'datInfo.mat', 'datInfo'));
end

end

% EOF