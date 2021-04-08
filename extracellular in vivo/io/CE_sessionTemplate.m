function session = CE_sessionTemplate(s, varargin)

% creates a session metadata struct in the cell explorer format. based on
% sessionTemplate. s can be a session struct or basepath
%
% INPUT:
%   basepath    string. path to recording folder {pwd}
%   importCh    logical. import skipped / synced channel from xml
%   viaGUI      logical. open CR session GUI {true}
%   force       logical. force default params set here to session even if
%               they already exist {false}
%   saveVar     logical. save variable {true}
%
% CALLS:
%   xmltree
%   catDat
%   preprocDat
%   getDinOE
%
% TO DO LIST:
%
% 09 jul 20 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 's', @(X) (ischar(X) && exist(X, 'dir')) || isstruct(X));
addParameter(p, 'importCh', true, @islogical);
addParameter(p, 'viaGUI', true, @islogical);
addParameter(p, 'force', true, @islogical);
addParameter(p, 'saveVar', true, @islogical);

% Parsing inputs
parse(p, s, varargin{:})
importCh = p.Results.importCh;
viaGUI = p.Results.viaGUI;
force = p.Results.force;
saveVar = p.Results.saveVar;

% Initializing session struct and defining basepath
if ischar(s)
    basepath = s;
    cd(basepath)
elseif isstruct(s)
    session = s;
    if isfield(session.general, 'basePath') &&...
            exist(session.general.basePath, 'dir')
        basepath = session.general.basePath;
        cd(basepath)
    else
        basepath = pwd;
    end
end

% initialize other params
session.channelTags = [];

% load existing basename.session.mat file if exist
[~, basename, ~] = fileparts(basepath);
sessionName = fullfile(basepath, [basename, '.session.mat']);
if ~exist('session', 'var') && exist(sessionName, 'file')
    sprintf('loading %s', sessionName)
    load(sessionName)
elseif ~exist('session', 'var')
    session = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standard params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumes file structure: animal/experiment/
pathPieces = regexp(basepath, filesep, 'split');

% General metadata
session.general.basePath    =  basepath; % Full path
session.general.name        = pathPieces{end}; % Session name / basename
session.general.version     = 5; % Metadata version
session.general.sessionType = 'Chronic'; % Type of recording: Chronic, Acute
session.general.nsamps      = 0;        % see end of script

% limited animal metadata
if ~isfield(session, 'animal') || force
    session.animal.name         = pathPieces{end - 1};
    session.animal.sex          = 'Male';
    session.animal.species      = 'Mouse';
    session.animal.strain       = 'C57BL';
    session.animal.geneticLine  = 'WT';
end

if ~isfield(session.general, 'experimenters') ||...
        isempty(session.general.experimenters) || force
    session.general.experimenters = 'LH';
end
if isfield(session.general, 'notes') || force
    session.general.notes = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extracellular
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default params (will be replaced later by xml params)
if ~isfield(session, 'extracellular') ||...
        (isfield(session, 'extracellular') && (~isfield(session.extracellular, 'sr')) ||...
        isempty(session.extracellular.sr)) || force
    session.extracellular.sr = 20000;           % Sampling rate
    session.extracellular.nChannels = 64;       % number of channels
    session.extracellular.fileName = '';        % (optional) file name of raw data if different from basename.dat
    session.extracellular.electrodeGroups.channels = {[1:session.extracellular.nChannels]}; % default list. change if needed
    session.extracellular.nElectrodeGroups = numel(session.extracellular.electrodeGroups);
    session.extracellular.spikeGroups = session.extracellular.electrodeGroups;
    session.extracellular.nSpikeGroups = session.extracellular.nElectrodeGroups;
end
if ~isfield(session, 'extracellular') ||...
        (isfield(session, 'extracellular') && (~isfield(session.extracellular, 'leastSignificantBit')) ||...
        isempty(session.extracellular.leastSignificantBit)) || force
    session.extracellular.leastSignificantBit = 0.195; % [�V]. Intan = 0.195, Amplipex = 0.3815
end
if ~isfield(session, 'extracellular') ||...
        (isfield(session, 'extracellular') && (~isfield(session.extracellular,'probeDepths')) ||...
        isempty(session.extracellular.probeDepths)) || force
    session.extracellular.probeDepths = 0;
end
if ~isfield(session, 'extracellular') ||...
        (isfield(session, 'extracellular') && (~isfield(session.extracellular, 'precision')) ||...
        isempty(session.extracellular.precision)) || force
    session.extracellular.precision = 'int16';
end

% Spike sorting
if ~isfield(session, 'spikeSorting')
    session.spikeSorting{1}.format = 'Neurosuite'; % Sorting data-format
    session.spikeSorting{1}.method = 'KiloSort'; % Sorting algorithm
    % generate relative path to clustered data (default is Kilosort output
    % folder generated by KiloSortWrapper)
    kiloSortFolder = dir('Kilosort_*');
    dirFlags = [kiloSortFolder.isdir];
    kiloSortFolder = kiloSortFolder(dirFlags);
    if ~isempty(kiloSortFolder)
        session.spikeSorting{1}.relativePath = kiloSortFolder.name;
    else
        session.spikeSorting{1}.relativePath = ''; % Relative path to the clustered data (here assumed to be the basepath)
    end
    session.spikeSorting{1}.channels = [];
    session.spikeSorting{1}.manuallyCurated = 1;
    session.spikeSorting{1}.notes = '';
end

% Brain regions
% Brain regions  must be defined as index 1. Can be specified on a channel
% or electrode group basis (below example for CA1 across all channels)
% session.brainRegions.CA1.channels = 1:128; % Brain region acronyms from
% Allan institute: http://atlas.brain-map.org/atlas?atlas=1)
% session.brainRegions.CA1.electrodeGroups = 1:8;

% Channel tags
% Channel tags must be defined as index 1. Each tag is a
% fieldname with the channels or spike groups as subfields. Below examples
% shows 5 tags (Theta, Ripple, RippleNoise, Cortical, Bad)
% session.channelTags.Theta.channels = 64;              % Theta channel
% session.channelTags.Ripple.channels = 64;             % Ripple channel
% session.channelTags.RippleNoise.channels = 1;         % Ripple Noise reference channel
% session.channelTags.Cortical.electrodeGroups = 3;     % Cortical spike groups
% session.channelTags.Bad.channels = 1;                 % Bad channels
% session.channelTags.Bad.electrodeGroups = 1;          % Bad spike groups (broken shanks)

% Analysis tags
if ~isfield(session, 'analysisTags') ||...
        (isfield(session, 'analysisTags') && (~isfield(session.analysisTags, 'probesLayout')) ||...
        isempty(session.analysisTags.probesLayout)) || force
    session.analysisTags.probesLayout = 'staggered'; % Probe layout: linear,staggered,poly2,edge,poly3,poly5
    session.analysisTags.probesVerticalSpacing = 10; % [�m] Vertical spacing between sites.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmlName = fullfile(session.general.basePath, [session.general.name, '.xml']);
sessionInfoName = fullfile(session.general.basePath, [session.general.name, '.sessionInfo.mat']);

% load from bzcode sessionInfo
if exist(sessionInfoName, 'file')
    load(sessionInfoName, 'sessionInfo')
    if sessionInfo.spikeGroups.nGroups > 0
        session.extracellular.nSpikeGroups = sessionInfo.spikeGroups.nGroups; % Number of spike groups
        session.extracellular.spikeGroups.channels = sessionInfo.spikeGroups.groups; % Spike groups
    else
        warning('No spike groups exist in the xml. Anatomical groups used instead')
        session.extracellular.nSpikeGroups = size(sessionInfo.AnatGrps, 2); % Number of spike groups
        session.extracellular.spikeGroups.channels = {sessionInfo.AnatGrps.Channels}; % Spike groups
    end
    session.extracellular.nElectrodeGroups = size(sessionInfo.AnatGrps,2); % Number of electrode groups
    session.extracellular.electrodeGroups.channels = {sessionInfo.AnatGrps.Channels}; % Electrode groups
    session.extracellular.sr = sessionInfo.rates.wideband; % Sampling rate of dat file
    session.extracellular.srLfp = sessionInfo.rates.lfp; % Sampling rate of lfp file
    session.extracellular.nChannels = sessionInfo.nChannels; % Number of channels
    % Changing index from 0 to 1:
    session.extracellular.electrodeGroups.channels=cellfun(@(x) x+1,session.extracellular.electrodeGroups.channels,'un',0);
    session.extracellular.spikeGroups.channels=cellfun(@(x) x+1,session.extracellular.spikeGroups.channels,'un',0);
    % load from xml
elseif exist('LoadXml.m', 'file') && exist(xmlName, 'file')
    sessionInfo = LoadXml(xmlName);
    if isfield(sessionInfo, 'SpkGrps')
        session.extracellular.nSpikeGroups = length(sessionInfo.SpkGrps); % Number of spike groups
        session.extracellular.spikeGroups.channels = {sessionInfo.SpkGrps.Channels}; % Spike groups
    else
        warning('No spike groups exist in the xml. Anatomical groups used instead')
        session.extracellular.nSpikeGroups = size(sessionInfo.AnatGrps, 2); % Number of spike groups
        session.extracellular.spikeGroups.channels = {sessionInfo.AnatGrps.Channels}; % Spike groups
    end
    session.extracellular.nElectrodeGroups = size(sessionInfo.AnatGrps,2); % Number of electrode groups
    session.extracellular.electrodeGroups.channels = {sessionInfo.AnatGrps.Channels}; % Electrode groups
    session.extracellular.sr = sessionInfo.SampleRate; % Sampling rate of dat file
    session.extracellular.srLfp = sessionInfo.lfpSampleRate; % Sampling rate of lfp file
    session.extracellular.nChannels = sessionInfo.nChannels; % Number of channels
    % Changing index from 0 to 1:
    session.extracellular.electrodeGroups.channels=cellfun(@(x) x+1,session.extracellular.electrodeGroups.channels,'un',0);
    session.extracellular.spikeGroups.channels=cellfun(@(x) x+1,session.extracellular.spikeGroups.channels,'un',0);
else
    warning('No sessionInfo.mat or xml file loaded')
    sessionInfo = [];
end

% if (~isfield(session.general, 'date') ||...
%         isempty(session.general.date)) && isfield(sessionInfo, 'Date') || force
%     if ~isempty(sessionInfo)
%         session.general.date = sessionInfo.Date;
%     end
% end
if isfield(session,'extracellular') &&...
        isfield(session.extracellular,'nChannels') || force
    fullpath = fullfile(session.general.basePath,[session.general.name,'.dat']);
    if exist(fullpath,'file')
        temp2_ = dir(fullpath);
        session.extracellular.nSamples = temp2_.bytes/session.extracellular.nChannels/2;
        session.general.duration = session.extracellular.nSamples/session.extracellular.sr;
    end
end

% Importing channel tags from sessionInfo
if isfield(sessionInfo, 'badchannels')
    if isfield(session.channelTags, 'Bad')
        session.channelTags.Bad.channels = unique([session.channelTags.Bad.channels,sessionInfo.badchannels+1]);
    else
        session.channelTags.Bad.channels = sessionInfo.badchannels+1;
    end
end

if isfield(sessionInfo, 'channelTags')
    tagNames = fieldnames(sessionInfo.channelTags);
    for iTag = 1:length(tagNames)
        if isfield(session, 'channelTags') && isfield(session.channelTags,tagNames{iTag})
            session.channelTags.(tagNames{iTag}).channels = unique([session.channelTags.(tagNames{iTag}).channels,sessionInfo.channelTags.(tagNames{iTag})+1]);
        else
            session.channelTags.(tagNames{iTag}).channels = sessionInfo.channelTags.(tagNames{iTag})+1;
        end
    end
end

% Importing brain regions from sessionInfo
if isfield(sessionInfo, 'region')
    load BrainRegions.mat
    regionNames = unique(cellfun(@num2str,sessionInfo.region,'uni',0));
    regionNames(cellfun('isempty',regionNames)) = [];
    for iRegion = 1:length(regionNames)
        if any(strcmp(regionNames{iRegion},BrainRegions(:,2)))
            session.brainRegions.(regionNames{iRegion}).channels = find(strcmp(regionNames{iRegion},sessionInfo.region));
        elseif strcmp(lower(regionNames{iRegion}),'hpc')
            session.brainRegions.HIP.channels = find(strcmp(regionNames{iRegion},sessionInfo.region));
        else
            disp(['Brain region does not exist in the Allen Brain Atlas: ' regionNames{iRegion}])
            regionName = regexprep(regionNames{iRegion}, {'[%() ]+', '_+$'}, {'_', ''});
            tagName = ['brainRegion_', regionName];
            if ~isfield(session,'channelTags') || all(~strcmp(tagName, fieldnames(session.channelTags)))
                disp(['Creating a channeltag with assigned channels: ' tagName])
                session.channelTags.(tagName).channels = find(strcmp(regionNames{iRegion},sessionInfo.region));
            end
        end
    end
end

% Epochs derived from MergePoints
if exist(fullfile(basepath,[session.general.name,'.MergePoints.events.mat']),'file')
    load(fullfile(basepath,[session.general.name,'.MergePoints.events.mat']),'MergePoints')
    for i = 1:size(MergePoints.foldernames,2)
        session.epochs{i}.name =  MergePoints.foldernames{i};
        session.epochs{i}.startTime =  MergePoints.timestamps(i,1);
        session.epochs{i}.stopTime =  MergePoints.timestamps(i,2);
    end
end

% Importing time series from intan metadatafile info.rhd
% session = loadIntanMetadata(session);

% % Importing Neuroscope xml parameters (skipped channels, dead channels, notes and experimenters)
% if importCh && exist(xmlName, 'file')
%     [sessionInfo, rxml] = LoadXml(xmlName);
%     
%     % Removing dead channels by the skip parameter in the xml file
%     if importSkippedChannels
%         order = [sessionInfo.AnatGrps.Channels];
%         skip = find([sessionInfo.AnatGrps.Skip]);
%         badChannels_skipped = order(skip)+1;
%     else
%         badChannels_skipped = [];
%     end
%     % Removing dead channels by comparing AnatGrps to SpkGrps in the xml file
%     if importSyncedChannels & isfield(sessionInfo,'SpkGrps')
%         skip2 = find(~ismember([sessionInfo.AnatGrps.Channels], [sessionInfo.SpkGrps.Channels])); % finds the indices of the channels that are not part of SpkGrps
%         badChannels_synced = order(skip2)+1;
%     else
%         badChannels_synced = [];
%     end
%     
%     if isfield(session,'channelTags') & isfield(session.channelTags,'Bad')
%         session.channelTags.Bad.channels = unique([session.channelTags.Bad.channels,badChannels_skipped,badChannels_synced]);
%     else
%         session.channelTags.Bad.channels = unique([badChannels_skipped,badChannels_synced]);
%     end
%     % Importing notes
%     try
%         if isfield(session.general,'notes')
%             session.general.notes = [session.general.notes,'   Notes from xml: ',rxml.child(1).child(3).value,'   Description: ' rxml.child(1).child(4).value];
%         else
%             session.general.notes = ['Notes: ',rxml.child(1).child(3).value,'   Description from xml: ' rxml.child(1).child(4).value];
%         end
%     end
%     % Importing experimenters
%     try
%         if ~isfield(session.general,'experimenters') || isempty(session.general.experimenters)
%             session.general.experimenters = rxml.child(1).child(2).value;
%         end
%     end
% end

% nsamps
[~, basename] = fileparts(basepath);
nchans = length([session.extracellular.electrodeGroups.channels{:}]);
nbytes = class2bytes('int16');
fraw = dir([basename '.dat']);
if ~isempty(fraw)
    nsamps = fraw.bytes / nbytes / nchans;
    session.general.nsamps = nsamps;
end

% Finally show GUI if requested by user
if viaGUI
    session = gui_session(session);
end

if saveVar
    save(sessionName, 'session')
end

end

% EOF