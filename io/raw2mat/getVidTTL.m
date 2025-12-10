
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TTL input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% event stream format can be found here:
% https://open-ephys.github.io/gui-docs/User-Manual/Recording-data/Binary-format.html

% files and paths
basepath = 'K:\Data\lh130\lh130_230321_084932';
oepath = 'K:\Data\lh130\2023-03-21_08-49-32';
expIdx = 2;
vidPath = 'K:\Data\Video\lh130';
vidFile = 'Cam_1__23339300__20230321_140854683.mp4';
cd(basepath)
[~, basename] = fileparts(basepath);

% get json file
jsonFile = dir([oepath filesep '**' filesep '*oebin']);
jsonFile = jsonFile(expIdx);
jsonName = fullfile(jsonFile.folder, jsonFile.name);

% map event file
mDin = load_open_ephys_binary(jsonName, 'events', 1, 'mmap');
tstamps_raw = double(mDin.Timestamps);
fs = mDin.Header.sample_rate;

% get indices of relevant channel
vidCh = 2;      % TTL channel in OE events stream
if length(unique(mDin.ChannelIndex)) > 1
    warning('more then one TTL was provided')
end
idxCh = mDin.ChannelIndex == vidCh;
if isempty(idxCh) 
    error('check video TTL channel')
end

% get indices of TTL triggers (rising or falling phase)
trigType = 'rise';
switch trigType
    case 'rise'
        idxTrig = mDin.Data > 0;
    case 'fall'
        idxTrig = mDin.Data < 0;
end

% check that tstamps fits video acquisition params 
fps = 50;                                               % frames per second [hz]
ttlDur = 5;                                             % ttl high duration [ms]
ttlPeriod = round(diff(tstamps_raw) / fs * 1000);           % time between ttls
uPeriod = unique(ttlPeriod);
if length(uPeriod) > 2 || sum(uPeriod) ~= 1 / fps * 1000 ||...
        uPeriod(1) ~= ttlDur
    warning('TTL duration does not fit acquisition params')
end

% get tstamps
tstamps = tstamps_raw(idxCh & idxTrig);
dt = diff(tstamps) / fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% video file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vid = VideoReader(fullfile(vidPath, vidFile));
nframes = vid.NumFrames;

dframes = nframes - length(tstamps);
if dframes ~= 0
    warning('%d mismatch between video frames and TTLs', dframes)
end

% fix dframes
ttl.tstamps = tstamps(1 : end - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust tstamps to entire recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dat info
load([basename, '.datInfo.mat'])

% session info
session = CE_sessionTemplate(basepath, 'viaGUI', false,...
    'forceDef', true, 'forceL', false, 'saveVar', false);      

% add recordign length to tstamps 
ttl.tstamps = ttl.tstamps + datInfo.nsamps(1);


% -------------------------------------------------------------------------
% save

ttl.info.header = mDin.Header;
ttl.info.jsonfile = jsonFile;
ttl.info.trigType = trigType;
ttl.info.vidCh = vidCh;
ttl.info.vid = vid;
ttl.info.dframes = dframes;
ttl.info.trigType = trigType;
ttl.info.remarks = 'missing 3 tstamps. added 358099968';
ttlfile = fullfile(basepath, [basename, '.vidTTL.mat']);

save(ttlfile, 'ttl')


