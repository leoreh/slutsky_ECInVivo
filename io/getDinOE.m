function din = getDinOE(varargin)

% arranges digital input from open ephys and saves as .mat in newpath.
% concatenates if multiple files exist. 
%
% INPUT:
%   basepath    string. path to .dat and .npy files {pwd}. 
%               if multiple dat files exist than concat must be true
%   newpath     string. path where new file should be saved. if empty than
%               new file will be save in basepath
%   newname     string. name of new file. if empty will be extracted from
%               newpath. in newpath not specified will be named as original
%               file with the extension '_new'
%   concat      logical. concatenate dat files {true} or not (false)
%   precision   char. sample precision {'int16'}
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
%   class2bytes
%
% TO DO LIST:
%   handle multiple channels (e.g. vid and stim)
%   options for rising and/or falling phase
%   compare .npy and .dat timestamps
%
% 09 apr 20 LH      
% 20 may 20 LH      input multiple paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'newpath', '', @ischar);
addOptional(p, 'newname', '', @ischar);
addOptional(p, 'concat', true, @islogical);
addOptional(p, 'precision', 'int16', @isstr);
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
nsamps = 0;
tstamps = cell(8, 1);

% set basepath to cell
basepath = cellstr(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle names for new path and new file
if isempty(newpath)
    newpath = basepath;
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
for i = 1 : length(basepath)
    datFiles(i) = dir([basepath{i} filesep '**' filesep '*dat']);
    jsonFiles(i) = dir([basepath{i} filesep '**' filesep '*oebin']);
end

if length(jsonFiles) ~= length(datFiles)
    warning('number of .oebin and .dat files not equal')
    return
end
if length(jsonFiles) > 1 && ~concat
    error(['multiple .oebin files found in %s\n',...
        'function cannot handle multiple files without concatination'], basepath)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(jsonFiles)
    
    % create file maps
    jsonName = fullfile(jsonFiles(i).folder, jsonFiles(i).name);
    mDin = load_open_ephys_binary(jsonName, 'events', 1, 'mmap');
    mDat = load_open_ephys_binary(jsonName, 'continuous', 1, 'mmap');
       
    % separate by channels and get rising phase
    chans = double(unique(mDin.ChannelIndex));
    idxRise = mDin.Data > 0;        % rising phase
    for ii = 1 : length(chans)
        idxCh = mDin.ChannelIndex == chans(ii);
        tstamps{ii} = [tstamps{ii}; mDin.Timestamps(idxRise & idxCh)];
    end
    
    % get nsamps of recording from corresponding .dat file
    nsamps(i) = datFiles(i).bytes / nbytes / nchans; 
    nsamps(i) = mDat.Timestamps(end);
    totsamps = cumsum(nsamps);
    
    % return if there are no events
    if all(cellfun(@isempty, tstamps))
        warning('no events found, skipping Din')
        return
    end
    
%     if data(end) > totsamps(end)
%         txt = sprintf(['digital input longer than recording duration\n',...
%             'check if argument nchans is correct']);
%         error(txt)
%     end
end

% arrange struct output and save
din.tstamps = tstamps;
din.origDir = {jsonFiles.folder};
din.fs = mDin.Header.sample_rate;
din.chans = chans;

if saveVar
    save(destination, 'din');
end

end

% EOF