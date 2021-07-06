function [EMG, EEG, sigInfo] = as_prepSig(eegData, emgData, varargin)

% prepare signals for accusleep. filters and subsamples eeg and emg. two
% options to input data; ALT 1 is to directly input eeg and emg data as
% vectors (e.g. via getLFP). this requires inputing the signals sampling
% frequency. ALT 2 is to load traces from an lfp / dat file. recommend
% using session info struct (cell explorer format) and inputing the signal
% channels within the file.
%
% INPUT:
%   eegData     ALT 1; eeg numeric data (1 x n1).
%               ALT 2; string. name of file with data. must include
%               extension (e.g. lfp / dat)
%   emgData     ALT 1; emg numeric data (1 x n2).
%               ALT 2; string. name of file with data. must include
%               extension (e.g. lfp / dat)
%   eegFs       numeric. eeg sampling frequency
%   emgFs       numeric. emg sampling frequency
%   eegCh       numeric. channel number of eeg to load from lfp file. can
%               be a vector and then the channels will be averaged
%   emgCh       numeric. channel number of eeg to load from lfp file. for
%               oe recording system
%   eegNchans   numeric. no channels in eeg data file. if empty will be
%               extracted from session info file
%   emgNchans   numeric. no channels in emg data file. if empty will equal
%               to eegNchans
%   basepath    string. path to recording folder {pwd}
%   saveVar     logical. save ss var {true}
%   forceLoad   logical. reload recordings even if mat exists
%   inspectSig  logical. inspect signals via accusleep gui {false}
%
% DEPENDENCIES:
%   rmDC
%   iosr.DSP.SINCFILTER     for low-pass filtering EEG data
%
% TO DO LIST:
%       # implement cleanSig
%       # input nchans for emg / eeg files separately
%
% 19 apr 21 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'eegCh', 1, @isnumeric);
addOptional(p, 'emgCh', 1, @isnumeric);
addOptional(p, 'eegFs', [], @isnumeric);
addOptional(p, 'emgFs', [], @isnumeric);
addOptional(p, 'eegNchans', [], @isnumeric);
addOptional(p, 'emgNchans', [], @isnumeric);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'inspectSig', false, @islogical);
addOptional(p, 'forceLoad', false, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
eegCh           = p.Results.eegCh;
emgCh           = p.Results.emgCh;
eegFs           = p.Results.eegFs;
emgFs           = p.Results.emgFs;
eegNchans       = p.Results.eegNchans;
emgNchans       = p.Results.emgNchans;
saveVar         = p.Results.saveVar;
inspectSig      = p.Results.inspectSig;
forceLoad       = p.Results.forceLoad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analysis params (decided by RA, HB, and LH 20 apr 21)
eegCf = 60;         % cutoff frequency for eeg
emgCf = [10 600];   % cutoff frequency for emg
fs = 1250;          % requested sampling frequency

% file names
cd(basepath)
mousepath = fileparts(basepath);
[~, basename] = fileparts(basepath);
eegfile = [basename '.AccuSleep_EEG.mat'];
emgfile = [basename '.AccuSleep_EMG.mat'];
sigInfofile = [basename '.AccuSleep_sigInfo.mat'];
sessionInfoFile = [basename, '.session.mat'];

% initialize
sigInfo = [];
recDur = [];

% reload data if already exists and return
if exist(emgfile, 'file') && exist(eegfile, 'file') && ~forceLoad
    fprintf('\n%s and %s already exist. loading...\n', emgfile, eegfile)
    load(eegfile, 'EEG')
    load(emgfile, 'EMG')
    return
end
    
% import toolbox for filtering    
import iosr.dsp.*

% session info
if exist(sessionInfoFile, 'file')
    load([basename, '.session.mat'])
    recDur = session.general.duration;
    if isempty(eegNchans)
        eegNchans = session.extracellular.nChannels;
    end
    if isempty(eegFs)
        eegFs = session.extracellular.srLfp;
    end
    if isempty(emgFs)
        emgFs = eegFs;
    end
end
if isempty(emgNchans)
    emgNchans = eegNchans;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eeg
fprintf('\nworking on %s\n', basename)
fprintf('loading data...')
if ischar(eegData) || isempty(eegData)
    if isempty(eegData)
        eegData = [basename '.lfp'];
        if ~exist(eegData, 'file')
            eegData = [basename '.dat'];
            if ~exist(eegData, 'file')
                error('could not fine %s. please specify lfp file or input data directly', eegfile)
            end
        end
    end
        
    % load average across given channels
    eegOrig = double(bz_LoadBinary(eegData, 'duration', Inf,...
        'frequency', eegFs, 'nchannels', eegNchans, 'start', 0,...
        'channels', eegCh, 'downsample', 1));
    if size(eegOrig, 2) > 1
        eegOrig = mean(eegOrig, 2);
    end
else
    eegOrig = eegData;
end

% emg
if ischar(emgData) || isempty(emgData)
    if isempty(emgData)
        emgData = [basename '.lfp'];
        if ~exist(emgData, 'file')
            emgData = [basename '.dat'];
            if ~exist(emgData, 'file')
                error('could not fine %s. please specify lfp file or input data directly', eegfile)
            end
        end
    end
        
    emgOrig = double(bz_LoadBinary(emgData, 'duration', Inf,...
        'frequency', emgFs, 'nchannels', emgNchans, 'start', 0,...
        'channels', emgCh, 'downsample', 1));        
else
    emgOrig = emgData;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter and downsample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low-pass filter eeg to assure nyquist.
% note accusleep only uses spectrogram up to 50 Hz
if ~isempty(eegCf)
    fprintf('\nfiltering EEG, cutoff = %d Hz', eegCf)
    filtRatio = eegCf / (eegFs / 2);
    eegOrig = iosr.dsp.sincFilter(eegOrig, filtRatio);
end

if ~isempty(emgCf)
    fprintf('\nfiltering EMG, cutoff = %d Hz', emgCf)
    filtRatio = emgCf / (emgFs / 2);
    emgOrig = iosr.dsp.sincFilter(emgOrig, filtRatio);
end

% remove DC component
fprintf('\nremoving DC component\n')
eegOrig = rmDC(eegOrig, 'dim', 1);
if numel(emgCf) ~= 2
    emgOrig = rmDC(emgOrig, 'dim', 1);
end

if isempty(recDur)
    recDur = length(emgOrig) / emgFs;
end

% validate recording duration and sampling frequency 
emgDur = length(emgOrig) / emgFs;
eegDur = length(eegOrig) / eegFs;
if abs(emgDur - eegDur) > 2
    warning(['EEG and EMG are of differnet duration (diff = %.2f s).\n',...
        'Check data and sampling frequencies.\n'], abs(emgDur - eegDur))
end
if isempty(recDur)
    recDur = eegDur;
end
tstamps_sig = [1 / fs : 1 / fs : recDur];

% subsample emg and eeg to the same length. assumes both signals span the
% same time interval. interpolation, as opposed to idx subsampling, is
% necassary for cases where the sampling frequency is not a round number
% (tdt). currently, this is done even when fs == round(fs) for added
% consistancy.
fprintf('downsampling to %d Hz\n', fs)
EMG = [interp1([1 : length(emgOrig)] / emgFs, emgOrig, tstamps_sig,...
    'pchip')]';
EEG = [interp1([1 : length(eegOrig)] / eegFs, eegOrig, tstamps_sig,...
    'pchip')]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finilize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigInfo.newFs = fs;
sigInfo.eegFs = eegFs;
sigInfo.emgFs = emgFs;
sigInfo.eegCf = eegCf;
sigInfo.emgCf = emgCf;
sigInfo.emgCh = emgCh;
sigInfo.eegCh = eegCh;
if ischar(emgData)
    sigInfo.emgFile = emgData;
else
    sigInfo.emgFile = 'inputData';
end
if ischar(eegData)
    sigInfo.eegFile = eegData;
else
    sigInfo.eegFile = 'inputData';
end

% save files
if saveVar
    fprintf('saving signals... ')
    save(eegfile, 'EEG')
    save(emgfile, 'EMG')
    save(sigInfofile, 'sigInfo')
    fprintf('done.\nthat took %.2f s\n\n', toc)
end

% inspect signals
if inspectSig
    AccuSleep_viewer(EEG, EMG, fs, 1, [], [])
end

return

% EOF

% -------------------------------------------------------------------------
% interpolate labels to a different fs
oldFs = 1000;
newFs = 1250;

% subsample
labelsNew = labels(1 : newFs / oldFs : end)

% interpolate
labelsNew = interp1([1 : floor(length(EEG) / oldFs)], labels, [1 : floor(length(EEG) / newFs)]);
labelsNew = round(labelsNew);       % necassary only if upsampling
labels = labelsNew;
