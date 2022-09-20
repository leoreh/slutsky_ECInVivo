function cat_OE_tstamps(varargin)

% creates tstamps.dat by concatenating timestamps.npy in OE dirs (assumes
% one file in each dat folder). for each file also validates that the
% length of tstamps and dat is equal. if not then removes samples from dat.
% this is based on a bug found in OE binary format (see validate_OE_tstamps).
%
% INPUT:
%   orig_paths  char / cell of chars. full path where the original timestamps files
%               exist. will take all files in these paths and concatenate
%               them in the same order as they are naturally sorted.
%   new_path    string. path where new file should be save. if empty than
%               new file will be saved in orig_paths{1}
%   new_name    string. name of new file (not including path). if empty
%               will be extracted from new_path. if new_path not specified
%               will be named as original file with the extension '_new'
%   precision   char. sample precision {'int16'}
%   nchans      numeric. number of channels in dat file {35}. needed to
%               determine the number of samples in each datfile
%
% OUTPUT
%   tstamps    
%
% CALLS:
%   class2bytes
%   datInfo
%   validate_OE_tstamps
%
% TO DO LIST:
%   # adapt for linux
%   # change tstamps to binary file
%
% 09 apr 20 LH
% 22 apr 20 LH      valTstampsOE
% 20 may 20 LH      input multiple paths
% 19 aug 20 LH      map to tstamps instead of direct readNPY
% 21 dec 21 LH      separated catDat
% 27 dec 21 LH      changed tstamps to binary file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'orig_paths', pwd);
addOptional(p, 'new_path', '', @ischar);
addOptional(p, 'new_name', '', @ischar);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'nchans', [], @isnumeric);

parse(p, varargin{:})
orig_paths = p.Results.orig_paths;
new_path = p.Results.new_path;
new_name = p.Results.new_name;
precision = p.Results.precision;
nchans = p.Results.nchans;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% size of one data point in bytes
nbytes = class2bytes(precision);

% initialize
datFiles = [];
tFiles = [];

% set orig_paths to cell
orig_paths = cellstr(orig_paths);

% handle names for new path and new file
if isempty(new_path)
    new_path = orig_paths{1};
end
if isempty(new_name)
    [~, basename] = fileparts(new_path);
    new_name = fullfile(new_path, [basename '.tstamps.dat']);
end

% open file
% fidTstamps = fopen(new_name, 'w');
% if(fidTstamps == -1)
%     error('cannot open file');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange files and concatenate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get .dat and .npy files in orig_paths
for ifile = 1 : length(orig_paths)
    datFiles = [datFiles; dir([orig_paths{ifile} filesep '**' filesep '*dat'])];
    tFiles = [tFiles; dir([orig_paths{ifile} filesep 'continuous' filesep '**' filesep 'timestamps.npy'])];
end

% check .dat files integrity and concat tstamps
for ifile = 1 : length(datFiles)
    source{ifile} = fullfile(datFiles(ifile).folder, datFiles(ifile).name);
    nsamps(ifile) = datFiles(ifile).bytes / nbytes / nchans;
    if ~isequal(nsamps(ifile), round(nsamps(ifile)))
        error('incorrect nCh for file')
    end
   
    % load corresponding tstamps file
    idx = find(strcmp({tFiles.folder}, datFiles(ifile).folder));
    if length(idx) > 1
        warning(['more than one timestamps.npy found in %s\n',...
            'skipping tstamps concatination'], orig_paths)
        break
    elseif length(idx) < 1
        warning(['timestamps.npy not found in %s\n',...
            'skipping tstamps concatination'], orig_paths)
        break
    else
       
        % memory map to tstamps file
        tname = fullfile(tFiles(idx).folder, tFiles(idx).name);
        [arrayShape, dataType, ~, ~, totalHeaderLength, ~] =...
            readNPYheader(tname);
        tmap = memmapfile(tname, 'Format', {dataType, arrayShape(end:-1:1), 'd'},...
            'Offset', totalHeaderLength);
        raw = tmap.data;
        nsamps_t = length(raw.d);
        
        % fix bug in OE binary format
        if nsamps_t < nsamps(ifile)
            warning(['more samples than timestamps in %s\n'...
                'initializing valTstampsOE'], source{ifile})
            validate_OE_tstamps('basepath', tFiles(idx).folder, 'precision', precision,...
                'chunksize', 5e6, 'bkup', false, 'saveVar', true,...
                'nchans', nchans)
            nsamps(ifile) = nsamps_t;
        end
        
        % write to file
        %         fwrite(fidTstamps, raw.d, 'int64');
        
        % clear data
        clear raw
        clear tmap
       
    end
end

% rc = fclose(fidTstamps);
% if rc == 0
%     fprintf('that took %.2f minutes\n', toc / 60)
% else
%     fprintf('. Failed to create %s!', new_name)
% end

end

% EOF