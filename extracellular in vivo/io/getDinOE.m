function din = getDinOE(varargin)

% arranges digital input from open ephys and saves as .mat in newpath.
% concatenates if multiple files exist. 
%
% INPUT:
%   datpath     string. path to .dat and .npy files {pwd}. 
%               if multiple dat files exist than concat must be true
%   newpath     string. path where new file should be saved. if empty than
%               new file will be save in datpath
%   newname     string. name of new file. if empty will be extracted from
%               newpath. in newpath not specified will be named as original
%               file with the extension '_new'
%   concat      logical. concatenate dat files {true} or not (false)
%   nbytes      number of bytes in one data point in .dat file 
%   nchans      numeric. number of channels in dat file {35}.
%   saveVar     logical. save datInfo {true} or not (false).
%
% OUTPUT
%   din         struct with the following fields:
%       data        timestamps of rising phase
%       origDir     path to original folder w/ json file
%
% CALLS:
%   bz_BasenameFromBasepath
%   load_open_ephys_binary
%
% TO DO LIST:
%   handle multiple channels (e.g. vid and stim)
%   options for rising and/or falling phase
%   compare .npy and .dat timestamps
%
% 09 apr 20 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'datpath', pwd);
addOptional(p, 'newpath', '', @ischar);
addOptional(p, 'newname', '', @ischar);
addOptional(p, 'concat', true, @islogical);
addOptional(p, 'nbytes', 2, @isnumeric);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
datpath = p.Results.datpath;
newpath = p.Results.newpath;
newname = p.Results.newname;
concat = p.Results.concat;
nbytes = p.Results.nbytes;
nchans = p.Results.nchans;
saveVar = p.Results.saveVar;

% handle names for new path and new file
if isempty(newpath)
    newpath = datpath;
end
if isempty(newname)
    if ~isempty(newpath)
        newname = [bz_BasenameFromBasepath(newpath)];
    else
        [~, newname] = fileparts(datFiles(1).name);
    end
end
destination = [newpath, filesep, newname, '.din.mat'];

% find corresponding .dat and din.npy files
jsonFiles = dir([datpath filesep '**' filesep '*oebin']);
datFiles = dir([datpath filesep '**' filesep '*dat']);
if length(jsonFiles) ~= length(datFiles)
    error('number of .oebin and .dat files not equal')
end
if length(jsonFiles) > 1 && ~concat
    error(['multiple .oebin files found in %s\n',...
        'function cannot handle multiple files without concatination'], datPath)
end

% go over files, load and arrange
nsamps = 0;
data = [];
for i = 1 : length(jsonFiles)
    
    % create file maps
    jsonName = fullfile(jsonFiles(i).folder, jsonFiles(i).name);
    mDin = load_open_ephys_binary(jsonName, 'events', 1, 'mmap');
    mDat = load_open_ephys_binary(jsonName, 'continuous', 1, 'mmap');

    % play with data
    data = [data; mDin.Timestamps(mDin.Data > 0)];
    
%     for i = 1 : length(data)
%     idx = find(mDat.Timestamps == data(i));
%     v = mDat.Data.Data.mapped(14, idx - 4000 : idx + 4000);
%     figure; plot(v)
%     end
%     
    % get nsamps of recording from corresponding .dat file
    nsamps(i) = datFiles(i).bytes / nbytes / nchans; 
    nsamps(i) = mDat.Timestamps(end);
    totsamps = cumsum(nsamps);
%     if data(end) > totsamps(end)
%         txt = sprintf(['digital input longer than recording duration\n',...
%             'check if argument nchans is correct']);
%         error(txt)
%     end
end

% arrange struct output and save
din.data = data;
din.origDir = {jsonFiles.folder};
if saveVar
    save(destination, 'din');
end

end

% EOF