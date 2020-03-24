function lfp = getStim(varargin)

% this is a wrapper to get lfp signals and stim from tdt
%  
% INPUT
%   basename    string. filename of lfp file. if empty retrieved from
%               basepath. if .lfp should not include extension, if .wcp
%               should include extension
%   basepath    string. path to load filename and save output {pwd}
%   extension   load from {'lfp'} (neurosuite), 'abf', or 'wcp'.
%   forceL      logical. force reload {false}.
%   fs          numeric. requested sampling frequency {1250}
%   mapch       new order of channels {[]}.
%   rmvch       channels to remove (according to original order) {[]}
%   clip        array of mats indicating times to diregard from recording.
%               each cell corresponds to a block. for example:
%               clip{3} = [0 50; 700 Inf] will remove the first 50 s from
%               block-3 and the time between 700 s and the end of Block-3
%   saveVar     save variable {1}.
%   chavg       cell. each row contain the lfp channels you want to average
%   
% OUTPUT
%   lfp         structure with the following fields:
%   fs
%   fs_orig
%   extension
%   interval    
%   duration    
%   chans
%   timestamps 
%   data  
% 
% 01 apr 19 LH & RA
% 19 nov 19 LH          load mat if exists  
% 14 jan 19 LH          adapted for wcp and abf and resampling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'andname', '');
addOptional(p, 'blocks', 1, @isnumeric);
addOptional(p, 'mapch', [1 : 16], @isnumeric);
addOptional(p, 'rmvch', [], @isnumeric);
addOptional(p, 'clip', {}, @iscell);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'saveVar', true, @islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
andname = p.Results.andname;
blocks = p.Results.blocks;
mapch = p.Results.mapch;
rmvch = p.Results.rmvch;
clip = p.Results.clip;
fs = p.Results.fs;
saveVar = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stim
[info] = tdt2dat('basepath', basepath, 'store', 'Stim', 'blocks', blocks,...
    'chunksize', 300, 'mapch', [], 'rmvch', [], 'clip', clip);
lfp.stim = bz_LoadBinary([info.filename '.dat'], 'frequency', info.fs, 'start', 0,...
    'duration', Inf, 'nChannels', 1);
lfp.stim = resample(double(lfp.stim), fs, round(info.fs));

% lfp
[info] = tdt2dat('basepath', basepath, 'store', 'Raw1', 'blocks', blocks,...
    'chunksize', 300, 'mapch', mapch, 'rmvch', rmvch, 'clip', clip);
fs_orig = info.fs;

lfp.data = bz_LoadBinary([info.filename '.dat'], 'frequency', info.fs, 'start', 0,...
    'duration', Inf, 'nChannels', length(mapch));
lfp.data = resample(double(lfp.data), fs, round(info.fs));

lfp.stim = interp1(linspace(0, 1, length(lfp.stim)), double(lfp.stim), linspace(0, 1, length(lfp.data)))';

% arrange struct
lfp.fs = fs;
lfp.fs_orig = fs_orig;
lfp.blocks = blocks;
lfp.blockduration = info.blockduration;
lfp.mapch = mapch;
lfp.rmvch = rmvch;
lfp.filename = [info.filename '.' andname '.mat'];

% save variable
if saveVar   
    save(lfp.filename, 'lfp')
end

end
