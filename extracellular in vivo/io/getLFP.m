function lfp = getLFP(varargin)

% gets lfp from lfp file for each channel. can specify channels and
% intervals. can also average across channels.
%  
% INPUT
%   filename    string. filename of lfp file
%   basepath    string. path to load filename and save output {pwd}
%   force       logical. force reload {false}.
%   fs          numeric. sampling frequency {1250}
%   interval    numeric mat. list of intervals to read from lfp file [s]
%   chans       vec. channels to load
%   saveVar     save variable {1}.
%   chavg       cell. each row contain the lfp channels you want to average
%   
% OUTPUT
%   lfp         structure with the following fields:
%   fs          
%   interval    
%   duration    
%   chans
%   timestamps 
%   data  
% 
% 01 apr 19 LH & RA
% 19 nov 19 LH          load mat if exists  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'filename', '');
addOptional(p, 'force', false, @islogical);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'interval', [0 inf], @isnumeric);
addOptional(p, 'chans', [1 : 16], @isnumeric);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'chavg', {}, @iscell);

parse(p,varargin{:})
basepath = p.Results.basepath;
filename = p.Results.filename;
force = p.Results.force;
fs = p.Results.fs;
interval = p.Results.interval;
chans = p.Results.chans;
saveVar = p.Results.saveVar;
chavg = p.Results.chavg;

nchans = length(chans);

if isempty(interval)
    interval = [0 inf];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if file exists
cd(basepath)
if isempty(filename)
    [~, basename] = fileparts(basepath);
    filename = [basename '.lfp'];
end
if ~exist(filename)
    filename = [basename '.eeg'];
    if ~exist(filename)
        error('file %s does not exist', filename)
    end
end
if exist(fullfile(basepath, [filename '.mat']), 'file') && ~force
    load(fullfile(basepath, [filename '.mat']))
    return
end

lfp.data = bz_LoadBinary(filename, 'duration', diff(interval),...
    'frequency', fs, 'nchannels', nchans, 'start', interval(1),...
    'channels', chans, 'downsample', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signal average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(chavg)
    mlfp = zeros(size(chavg, 1), length(lfp.data));
    for i = 1 : size(chavg, 1)
        mlfp(i, :) = mean(lfp.data(:, chavg{i}), 2);
    end    
    lfp.data = mlfp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange and save struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lfp.timestamps = (interval(1) : 1 / fs : interval(1) + (length(lfp.data) - 1) / fs)';

if interval(2) == inf
    inteval(2) = lfp.timestamps(end);
end
lfp.interval = interval;
lfp.duration = length(lfp.data) / fs;
lfp.chans = chans;
lfp.fs = fs;

% save variable
if saveVar   
    save([basepath, filesep, filename, '.mat'], 'lfp')
end



end
