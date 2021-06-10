function as_inspectSeparation(EEG, EMG, labels, varargin)

% plot the separation of states based on emg (rms) and eeg (spectrogram)
% signals
%
% INPUT:
%   EMG             numeric. emg data (1 x n)
%   EEG             numeric. eeg data (1 x n)
%   labels          numeric. 
%   fs              numeric. sampling frequency
%   saveFig         logical. save figure {true}
%
% DEPENDENCIES
%   AccuSleep (modified in slutskycode)
%   IOSR.DSP.SINCFILTER     for filtering data
% 
% TO DO LIST
%   calc fft power through spectrogram
%
% 08 jun 21 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
fs              = p.Results.fs;
saveFig         = p.Results.saveFig;

basepath = pwd;
[~, basename] = fileparts(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constants
epochLen = 1;        
minBoutLen = epochLen;

% get params from configuration file
[cfg_colors, cfg_names, ~] = as_loadConfig([]);
nstates = length(cfg_names);

% validate data
if length(EEG) ~= length(EMG)
    error('EEG and EMG must be the same length')
end

if isempty(fs)
    fs = 1250;
end
newFs = 128;

tic;
fprintf('\nPreparing signals, this may take some time...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep signals (similar pipeline from AccuSleep_classify)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% standardize
eeg = standardizeSR(EEG, fs, newFs);
emg = standardizeSR(EMG, fs, newFs);

% calibration 
calibrationData = createCalibrationData(standardizeSR(EEG, fs, 128),...
    standardizeSR(EMG, fs, 128), labels, 128, epochLen);

% select relavent labels
labelsIdx = find(labels < nstates);

% create spectrogram
[s, ~, f] = createSpectrogram(eeg, newFs, epochLen);

% frequency indices
[~, f1idx] = min(abs(f - 1));
[~, f4idx] = min(abs(f - 4));
[~, f6idx] = min(abs(f - 6));
[~, f12idx] = min(abs(f - 12));
[~, f20idx] = min(abs(f - 20)); % index in f of 20Hz
[~, f50idx] = min(abs(f - 50)); % index in f of 50Hz

% select frequencies up to 50 Hz, and downsample between 20 and 50 Hz
s = s(:, [1 : (f20idx - 1), f20idx : 2 : f50idx]);

% take log
s = log(s);

% scale the spectrogram
for j = 1:size(s, 2)
    s(:, j) = (s(:, j) - calibrationData(j, 1)) ./ calibrationData(j, 2);
    s(:, j) = (s(:, j) + 4.5) ./ 9; % clip z scores
end

% clip
s(s < 0) = 0;
s(s > 1) = 1;

sDelta = sum(s(:, f1idx : f4idx), 2);
sTheta = sum(s(:, f6idx : f12idx), 2);
sRatio = sDelta ./ sTheta;

% calculate log rms for each EMG bin
processedEMG = processEMG(emg, newFs, epochLen);

% scale the EMG
processedEMG = (processedEMG - calibrationData(end, 1)) ./ calibrationData(end, 2);
processedEMG = (processedEMG + 4.5) ./ 9;

% clip
processedEMG(processedEMG < 0) = 0;
processedEMG(processedEMG > 1) = 1;

fprintf('\ndone in %.1f sec\n\n', toc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert labels to state epochs. 
for istate = 1 : nstates
    binaryVec = zeros(length(labels), 1);
    binaryVec(labels == istate) = 1;
    stateEpisodes = binary2epochs('vec', binaryVec, 'minDur', [], 'maxDur', [],...
        'interDur', [], 'exclude', false); % these are given as indices and are equivalent to seconds
    stateEpochs{istate} = stateEpisodes * epochLen;
end

% calculate power for each epoch separatly
win = hann(2 ^ (nextpow2(2 * newFs) - 1));
noverlap = 0.25 * newFs;
freqs = [0 : 0.2 : 50]; % freqencies for power analysis
fftSum = zeros(nstates - 1, length(freqs));
for istate = 1 : nstates - 1
    for iepoch = 1 : length(stateEpochs{istate})
        % get eeg data in epoch
        dataIdx = stateEpochs{istate}(iepoch, 1) * newFs :...
            stateEpochs{istate}(iepoch, 2) * newFs - 1;
        
        % calc power
        [pow, ~] = pwelch(eeg(dataIdx), win, noverlap, freqs, newFs);
        fftSum(istate, :) = fftSum(istate, :) + pow; 
               
        % eeg ratio high / low
        eegRatio{istate}(iepoch) = sum(pow(freqs > 25)) / sum(pow(freqs < 5));        
               
        % emg rms 
        emgRMS{istate}(iepoch) = rms(emg(dataIdx));
    end
    fftSum(istate, :) = fftSum(istate, :) / length(stateEpochs{istate});    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(false)
fh = figure;

% spectrogram vs. emg
subplot(2, 2, 1)
gscatter(sRatio(labelsIdx), processedEMG(labelsIdx),...
    labels(labelsIdx), [cell2nanmat(cfg_colors')]')
ylabel('Norm. EMG RMS')
xlabel('Delta / Theta Ratio')
legend(cfg_names)

% emg vs. time
subplot(2, 2, 2)
hold on
for istate = 1 : nstates - 1
    stateLabels = find(labels == istate);
    scatter(stateLabels / epochLen / 60 / 60,...
        processedEMG(stateLabels),...
        3, cfg_colors{istate})
end
axis tight
ylim([0 1])
xlabel('Time [h]')
ylabel('Norm. EMG RMS')

subplot(2, 2, 3)
ph = plot(freqs, fftSum, 'LineWidth', 3);
set(ph, {'color'}, cfg_colors(1 : nstates - 1))
legend(cfg_names)
xlim([0 20])
xlabel('Frequency [Hz]')
ylabel('Power [dB]')

subplot(2, 2, 4)
hold on
for istate = [1, 2, 4, 5]
    histogram(log10(eegRatio{istate}), 30, 'Normalization', 'probability')
end
ylabel('Norm. Counts')
xlabel('EEG High/Low Ratio')

if saveFig
    figpath = fullfile('graphics', 'sleepState');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_stateSeparation', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF
