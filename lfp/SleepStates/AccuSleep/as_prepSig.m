function sSig = as_prepSig(eegData, emgData, varargin)

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
%               extension (e.g. lfp / dat). If empty will be the same as
%               eegData
%   fs          numeric. requested new sampling frequency
%   eegFs       numeric. eeg sampling frequency
%   emgFs       numeric. emg sampling frequency
%   eegCh       numeric. channel number of eeg to load from lfp file. can
%               be a vector and then the channels will be averaged. 
%   emgCh       numeric. channel number of eeg to load from lfp file. for
%               oe recording system
%   emgCf       numeric. cut off frequency for emg signal [10 200] or acc
%               signal [10 450]. decided by RA, HB, and LH 20 apr 21
%   eegCf       numeric. low-pass frequency for eeg signal [60]
%   eegNchans   numeric. no channels in eeg data file. if empty will be
%               extracted from session info file
%   emgNchans   numeric. no channels in emg data file. if empty will equal
%               to eegNchans
%   sigWin      numeric 1 x 2. time to load in sec {[0 Inf]}
%   basepath    string. path to recording folder {pwd}
%   saveVar     logical. save signals as mat var {true}
%   forceLoad   logical. reload recordings even if mat exists
%   inspectSig  logical. inspect signals via accusleep gui {false}
%
% OUTPUT
%   sSig        struct
% 
% DEPENDENCIES:
%   iosr.DSP.SINCFILTER     for low-pass filtering EEG data
%
% TO DO LIST:
%       # implement cleanSig (irrelavent)
%       # input nchans for emg / eeg files separately (done)
%       # implement tsa_filt
%       # chunk processing
%
% 19 apr 21 LH  updates:
% 06 jan 22     combined signals to struct
% 25 feb 22     specBand instead of createSpectrogram
% 28 apr 22     ftarget starting at 0.5. must also consider rayleigh

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
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'eegNchans', [], @isnumeric);
addOptional(p, 'emgNchans', [], @isnumeric);
addOptional(p, 'emgCf', [10 450], @isnumeric);
addOptional(p, 'eegCf', [450], @isnumeric);
addOptional(p, 'sigWin', [0 Inf], @isnumeric);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'inspectSig', false, @islogical);
addOptional(p, 'forceLoad', false, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
eegCh           = p.Results.eegCh;
emgCh           = p.Results.emgCh;
eegFs           = p.Results.eegFs;
emgFs           = p.Results.emgFs;
fs              = p.Results.fs;
emgCf           = p.Results.emgCf;
eegCf           = p.Results.eegCf;
eegNchans       = p.Results.eegNchans;
emgNchans       = p.Results.emgNchans;
sigWin          = p.Results.sigWin;
saveVar         = p.Results.saveVar;
inspectSig      = p.Results.inspectSig;
forceLoad       = p.Results.forceLoad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sig info
if ischar(emgData)
    sSig.info.emgFile = emgData;
else
    sSig.info.emgFile = 'inputData';
end
if ischar(eegData)
    sSig.info.eegFile = eegData;
else
    sSig.info.eegFile = 'inputData';
end

% assert sigWin
if numel(sigWin) ~= 2 || sigWin(2) <= sigWin(1)
    error('check sigWin')
end

% file names
cd(basepath)
[~, basename] = fileparts(basepath);
sigfile = [basename '.sleep_sig.mat'];
sessionfile = [basename, '.session.mat'];

% reload data if already exists and return
if exist(sigfile, 'file') && ~forceLoad
    fprintf('\n%s already exists. loading...\n', sigfile)
    load(sigfile, 'info')
    return
end
    
% import toolbox for filtering    
import iosr.dsp.*

% session info and signal params
if exist(sessionfile, 'file')
    load(sessionfile)
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
        
    % load channels and average
    eegData = binary_load(eegData, 'duration', diff(sigWin),...
        'fs', eegFs, 'nCh', eegNchans, 'start', sigWin(1),...
        'ch', eegCh, 'downsample', 1);
    if size(eegData, 2) > 1
        eegData = mean(eegData, 2);
    end
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
        
    emgData = binary_load(emgData, 'duration', diff(sigWin),...
        'fs', emgFs, 'nCh', emgNchans, 'start', sigWin(1),...
        'ch', emgCh, 'downsample', 1);        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter and downsample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the eeg and emg signals must be precisely the same length. this becomes
% complicated when recording with tdt because of their non-integer sampling
% rate. if an .lfp file exists, it is best to leave the lfp data untouched
% (interpolation adds a strong slow component to the spectrogram). if an
% .lfp file was not created (e.g. raw eeg from tdt is used), then it is
% best to follow the pipeline in LFPfromDat (filter followed by
% subsampling). this cannot be used for the emg signal, which requires
% interpolation to be the same length as the eeg signal. this is ok if the
% emg is filtered after the interpolation. note that although accusleep
% only uses the spectrogram up to 50 hz, there is no benifit to low-pass
% filtering the eeg if the sampling frequency is maintained at 1250 (which
% is better for the emg rms signal).

% handle eeg data (must apply batch proccessing due to memory)
if eegFs ~= fs
    
    filtRatio = eegCf / (eegFs / 2);
    fsRatio = eegFs / fs;

    % low-pass filter
    if ~isempty(eegCf)
        fprintf('\nfiltering EEG, cutoff = %d Hz', eegCf)
        eegData = iosr.dsp.sincFilter(eegData, filtRatio);
    end
    
    % subsample
    fprintf('\nsubsampling eeg from %.2f to %d Hz', eegFs, fs)
    eegData = eegData(fsRatio : fsRatio : length(eegData));

else
    eegData = eegData;
end

% handle emg data
if emgFs ~= fs || length(emgData) ~= length(eegData)
    tstamps_sig = [1 : length(eegData)] / fs;
    
    % interpolate
    fprintf('\ninterpolating emg from %.2f to %d Hz', emgFs, fs)
    emgData = [interp1([1 : length(emgData)] / emgFs, emgData, tstamps_sig,...
        'spline')]';
end
% filter
if ~isempty(emgCf)
    fprintf('\nfiltering EMG, cutoff = %d Hz', emgCf)
    filtRatio = emgCf / (fs / 2);
    emgData = iosr.dsp.sincFilter(emgData, filtRatio);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare spectrogram and emg rms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spectrogram
fprintf('\ncreating spectrogram...')
spec = calc_spec('sig', eegData,...
    'fs', fs, 'graphics', false, 'saveVar', false, 'force', true,...
    'padfft', 0, 'winstep', 1, 'logfreq', false,...
    'ftarget', [0.5 : 0.2 : 20, 20.5 : 0.5 : 50]);
sSig.spec = spec.s;
sSig.spec_tstamps = spec.tstamps;
sSig.spec_freq = spec.freq;

% emg rms
sSig.emg_rms = processEMG(emgData, fs, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finilize and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% info
sSig.emg = emgData(:);
clear emgData
sSig.eeg = eegData(:);
clear eegData
sSig.fs = fs;
sSig.info.eegFs = eegFs;
sSig.info.emgFs = emgFs;
sSig.info.eegCf = eegCf;
sSig.info.emgCf = emgCf;
sSig.info.emgCh = emgCh;
sSig.info.eegCh = eegCh;
sSig.info.emgNchans = emgNchans;
sSig.info.eegNchans = eegNchans;
sSig.info.runtime = datetime("now");

% save files
if saveVar
    fprintf('\nsaving signals... ')
    save(sigfile, '-struct', 'sSig')
    fprintf('done.\nthat took %.2f s\n\n', toc)
end

% inspect signals
if inspectSig
    AccuSleep_viewer(sSig, [], [])
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
labelsNew = interp1([1 : floor(length(eegData) / oldFs)], labels, [1 : floor(length(eegData) / newFs)]);
labelsNew = round(labelsNew);       % necassary only if upsampling
labels = labelsNew;

% -------------------------------------------------------------------------
% call for lfp and emg from same file
sSig = as_prepSig([basename, '.lfp'], [],...
    'eegCh', [29 : 32], 'emgCh', 33, 'saveVar', true, 'eegNchans', nchans,...
    'inspectSig', false, 'forceLoad', true, 'eegFs', 1250, 'emgFs', [],...
    'emgCf', [10 450]);

% call for emg dat file
sSig = as_prepSig([basename, '.lfp'], [basename, '.emg.dat'],...
    'eegCh', [5 : 8], 'emgCh', 1, 'saveVar', true, 'emgNchans', 1, 'eegNchans', nchans,...
    'inspectSig', false, 'forceLoad', true, 'eegFs', 1250, 'emgFs', 1250,...
    'emgCf', [10 450]);

% call for accelerometer
sSig = as_prepSig([basename, '.lfp'], acc.mag,...
    'eegCh', [20 : 22], 'emgCh', [], 'saveVar', true, 'emgNchans', [],...
    'eegNchans', nchans, 'inspectSig', false, 'forceLoad', true,...
    'eegFs', 1250, 'emgFs', 1250, 'eegCf', [], 'emgCf', [10 450], 'fs', 1250);

% -------------------------------------------------------------------------
% fix special case where emg changes in the middle of a recording, such as
% the sudden appearance of strong HR artifacts. Thus, although wake and
% nrem are still noticeably different, the difference between them is not
% homegenous and the same calibration matrix cannot be used for the entire
% recording. so, classify the recording twice (or more) and each time
% during the calibration use only the manual (calibration) labels from the
% parts of the recording that look similar. then, combine the labels from
% the different classifications. The transitions will probably need to be
% manually labeled. 
