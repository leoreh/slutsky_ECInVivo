
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TTL input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% event stream format can be found here:
% https://open-ephys.github.io/gui-docs/User-Manual/Recording-data/Binary-format.html

% files and paths
basepath = 'E:\Data\lh129\lh129_230218_091553';
oepath = 'E:\Data\lh129\lh129_230218_091553\2023-02-18_11-15-53';
vidPath = 'H:\Data\Video\lh129';
vidFile = 'Cam_1__23339300__20230218_133616998.mp4';
cd(basepath)
[~, basename] = fileparts(basepath);

% get json file
jsonFile = dir([oepath filesep '**' filesep '*oebin']);
% jsonFile = jsonFile(2);
jsonName = fullfile(jsonFile.folder, jsonFile.name);

% map event file
mDin = load_open_ephys_binary(jsonName, 'events', 1, 'mmap');
tstamps = double(mDin.Timestamps);
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
fps = 50;       % frames per second [hz]
pulseDur = 5;   % ttl high duration [ms]
ttlDur = unique(round(diff(tstamps) / fs * 1000));
if length(ttlDur) > 2 || sum(ttlDur) ~= 1 / fps * 1000 ||...
        ttlDur(1) ~= pulseDur
    warning('TTL duration does not fit acquisition params')
end

% get tstamps
tstamps = tstamps(idxCh & idxTrig);
dt = diff(tstamps) / fs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% video file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vid = VideoReader(fullfile(vidPath, vidFile));
nframes = vid.NumFrames;

if nframes ~= length(tstamps)
    warning('mismatch between video frames and TTLs')
end

tstamps = tstamps(1 : end - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust tstamps to entire recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% dat info
load([basename, '.datInfo.mat'])

% session info
session = CE_sessionTemplate(basepath, 'viaGUI', false,...
    'forceDef', true, 'forceL', false, 'saveVar', false);      

% check that dat files in datInfo fit the recording file
% session.general.nsamps
% sum(datInfo.nsamps)

% add recordign length to tstamps 
recfile = 2;        % idx of relevant recording file
% ttl.tstamps = tstamps + datInfo.nsamps(recfile - 1);

ttl.tstamps = tstamps + 253547008;


109547008 + 144000000

1828790272 - 1684790272


% -------------------------------------------------------------------------
% save

ttl.info.header = mDin.Header;
ttl.info.jsonfile = jsonFile;
ttl.info.trigType = trigType;
ttl.info.vidCh = vidCh;
ttl.info.vid = vid;
ttl.info.trigType = trigType;
ttl.info.remarks = '1 stamp more than frames, removed last. added 253547008 to tstamps';
ttlfile = fullfile(basepath, [basename, '.vidTTL.mat']);

save(ttlfile, 'ttl')


