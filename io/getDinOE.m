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
addOptional(p, 'nchans', [], @isnumeric);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
newpath = p.Results.newpath;
newname = p.Results.newname;
concat = p.Results.concat;
precision = p.Results.precision;
nchans = p.Results.nchans;
saveVar = p.Results.saveVar;

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
    end
end
destination = [newpath, filesep, newname, '.din.mat'];

% find corresponding .dat and din.npy files
for ifile = 1 : length(basepath)
    tmp = dir([basepath{ifile} filesep '**' filesep '*dat']);
    if ~isempty(tmp)
        datFiles(ifile) = tmp;
    end
    tmp = dir([basepath{ifile} filesep '**' filesep '*oebin']);
    if ~isempty(tmp)
        jsonFiles(ifile) = tmp;
    end
end

if length(jsonFiles) ~= length(datFiles)
    warning('number of .oebin and .dat files not equal')
end
if length(jsonFiles) > 1 && ~concat
    error(['multiple .oebin files found in %s\n',...
        'function cannot handle multiple files without concatination'], basepath)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nCreating %s from files:\n', destination)

for ifile = 1 : length(jsonFiles)    
    
    if isempty(jsonFiles(ifile).name)
        continue
    end
    
    % create file maps
    jsonName = fullfile(jsonFiles(ifile).folder, jsonFiles(ifile).name);
    mDin = load_open_ephys_binary(jsonName, 'events', 1, 'mmap');
    
    fprintf('%s...\n', jsonName)
    
    % separate by channels and get rising phase
    chans = double(unique(mDin.ChannelIndex));
    idxRise = mDin.Data > 0;        % rising phase
    for ich = 1 : length(chans)
        idxCh = mDin.ChannelIndex == chans(ich);
        tstamps{ich} = [tstamps{ich}; mDin.Timestamps(idxRise & idxCh)];
    end
    
    % validate Din length by checking the nsamps of recording from
    % corresponding .dat file and. this can take a long time and is
    % typically not necassary. 
    if ~isempty(nchans)
        mDat = load_open_ephys_binary(jsonName, 'continuous', 1, 'mmap');
        nbytes = class2bytes(precision);
        nsamps(ifile) = datFiles(ifile).bytes / nbytes / nchans;
        nsamps(ifile) = mDat.Timestamps(end);
        totsamps = cumsum(nsamps);
        if data(end) > totsamps(end)
            txt = sprintf(['digital input longer than recording duration\n',...
                'check if argument nchans is correct']);
            error(txt)
        end
    end
    
    % return if there are no events
    if all(cellfun(@isempty, tstamps))
        warning('no events in %s', jsonName)
    end
    
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