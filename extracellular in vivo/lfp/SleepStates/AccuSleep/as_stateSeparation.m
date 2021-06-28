function as_stateSeparation(EEG, EMG, labels, varargin)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[~, basename] = fileparts(basepath);

% constants
epochLen = 1;        
minBoutLen = epochLen;

% get params from configuration file
[cfg_colors, cfg_names, ~] = as_loadConfig([]);
nstates = length(cfg_names);
sstates = 1 : nstates - 1;      % selected states (ignore bin)

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
% prep signals (similar pipeline to AccuSleep_classify)
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
[spec, taxis, f] = createSpectrogram(eeg, newFs, epochLen);

% frequency indices
[~, f1idx] = min(abs(f - 1));
[~, f4idx] = min(abs(f - 4));
[~, f6idx] = min(abs(f - 6));
[~, f12idx] = min(abs(f - 12));
[~, f20idx] = min(abs(f - 20)); % index in f of 20Hz
[~, f50idx] = min(abs(f - 50)); % index in f of 50Hz

% select frequencies up to 50 Hz, and downsample between 20 and 50 Hz
spec = spec(:, [1 : (f20idx - 1), f20idx : 2 : f50idx]);

% take log
spec = log(spec);

% scale the spectrogram
for j = 1:size(spec, 2)
    spec(:, j) = (spec(:, j) - calibrationData(j, 1)) ./ calibrationData(j, 2);
    spec(:, j) = (spec(:, j) + 4.5) ./ 9; % clip z scores
end

% clip
spec(spec < 0) = 0;
spec(spec > 1) = 1;

sDelta = sum(spec(:, f1idx : f4idx), 2);
sTheta = sum(spec(:, f6idx : f12idx), 2);
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
    epLen{istate} = [diff(stateEpochs{istate}')]';
end

% calculate power for each epoch separatly
win = hann(2 ^ (nextpow2(2 * newFs) - 1));
noverlap = 0.25 * newFs;
faxis = [0 : 0.2 : 50]; % freqencies for power analysis
fftSum = zeros(nstates - 1, length(faxis));
for istate = sstates
    for iepoch = 1 : length(stateEpochs{istate})
        % skip epochs shorter than one bin (may happen at end of recording)
        if epLen{istate}(iepoch) < 1 
            continue
        end
        
        % get eeg data in epoch
        dataIdx = stateEpochs{istate}(iepoch, 1) * newFs :...
            stateEpochs{istate}(iepoch, 2) * newFs - 1;
        
        % calc power
        [pow, ~] = pwelch(eeg(dataIdx), win, noverlap, faxis, newFs);
        fftSum(istate, :) = fftSum(istate, :) + pow; 
               
        % eeg ratio high / low
        eegRatio{istate}(iepoch) = sum(pow(faxis > 25)) / sum(pow(faxis < 5));        
               
        % emg rms 
        emgRMS{istate}(iepoch) = rms(emg(dataIdx));
    end
    fftSum(istate, :) = fftSum(istate, :) / length(stateEpochs{istate});    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(false)
fh = figure('Color', 'w', 'Position', [0.08, 0.04, 0.83, 0.88]);
sb1 = axes('Position' ,[0.08, 0.85, 0.83, 0.12]);   % emg vs time
sb2 = axes('Position' ,[0.08, 0.70, 0.83, 0.12]);
sb3 = axes('Position' ,[0.08, 0.42, 0.35, 0.24]);
sb4 = axes('Position' ,[0.55, 0.42, 0.35, 0.24]);
sb5 = axes('Position' ,[0.08, 0.08, 0.35, 0.24]);
sb6 = axes('Position' ,[0.55, 0.08, 0.35, 0.24]);
linkaxes([sb1, sb2], 'x');      % link spectrogram and emg rms
set([sb1, sb2, sb3, sb4, sb5, sb6], 'box', 'off', 'TickLength', [0 0])

% emg vs time colored by state
fh.CurrentAxes = sb1;
hold on
for istate = sstates
    stateLabels = find(labels == istate);
    scatter(stateLabels / epochLen / 60 / 60,...
        processedEMG(stateLabels),...
        3, cfg_colors{istate})
end
axis tight
ylim([min(processedEMG) 1])
ylabel('Norm. EMG RMS')
set(gca, 'XTick', [], 'YTick', [])

% spectrogram
fh.CurrentAxes = sb2;
% time axis in hours
nbins = length(taxis);
tSpec = ((1 : nbins) * epochLen - epochLen / 2) / 3600; % spectrogram time axis, in seconds
showFreqs = find(faxis <= 15);  % choose freqs to display
% take a sample of the spectrogram to help initialize the colormap
sampleBins = randperm(nbins, round(nbins / 10));
specSample = reshape(spec(sampleBins, showFreqs), 1, length(sampleBins) * length(showFreqs));
caxis = prctile(specSample, [6 99.5]);
% plot 
imagesc(tSpec, faxis(showFreqs), spec(:, showFreqs)', caxis);
colormap(AccuSleep_colormap());
axis('xy')
ylabel('Frequency [Hz]')
xlabel('Time [h]')
set(gca, 'YTick', [])

% spectral power per state
fh.CurrentAxes = sb3;
ph = plot(faxis, fftSum, 'LineWidth', 3);
set(ph, {'color'}, cfg_colors(sstates))
xlim([0 20])
xlabel('Frequency [Hz]')
ylabel('Power [dB]')     
set(gca, 'YTick', [])

% scatter spectrogram vs. emg (could not use axes handle w/ gscatter)
fh.CurrentAxes = sb4;
hold on
for istate = sstates
scatter(sRatio(labels == istate), processedEMG(labels == istate)',...
    2, cfg_colors{istate}, 'filled')
end
ylabel('Norm. EMG RMS')
xlabel('Delta / Theta Ratio')
set(gca, 'YTick', [])

% histogram of epoch lengs for NREM, REM, and WAKE
fh.CurrentAxes = sb5;
hold on
for istate = [1, 4, 5]
    hh = histogram(log10(epLen{istate}), 30, 'Normalization', 'probability');
    set(hh, 'EdgeColor', 'none', 'FaceColor', cfg_colors{istate}, 'FaceAlpha', 0.4)
end
ylabel('Norm. Counts')
xlabel('Epoch Length [log(s)]')
set(gca, 'YTick', [])

% percent time in state
fh.CurrentAxes = sb6;
pie(sum(cell2nanmat(epLen(sstates)), 1, 'omitnan'), ones(1, length(sstates)));
hold on
ph = findobj(sb6, 'Type', 'Patch');
set(ph, {'FaceColor'}, flipud(cfg_colors(sstates)))
legend(cfg_names, 'Units', 'normalized', 'Position', [0.81 0.11 0.10 0.20]);

% save figure
if saveFig
    figpath = fullfile('graphics', 'sleepState');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_stateSeparation', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF
