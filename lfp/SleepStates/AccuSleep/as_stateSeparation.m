function as_stateSeparation(sSig, ss, varargin)

% plot the separation of states based on emg (rms) and eeg (spectrogram)
% signals
%
% INPUT:
%   sSig            struct. see as_prepSig.m
%   ss          	struct. see as_classify.m 
%   saveFig         logical. save figure {true}
%
% DEPENDENCIES
%   AccuSleep (modified in slutskycode)
%   IOSR.DSP.SINCFILTER     for filtering data
% 
% TO DO LIST
%   calc fft power through spectrogram
%
% 08 jun 21 LH      updates:
% 12 jan 22 LH      ss as input instead of labels


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
saveFig         = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[~, basename] = fileparts(basepath);

% get params from configuration file
cfg = as_loadConfig();
epochLen = cfg.epochLen;
nstates = cfg.nstates;
sstates = 1 : nstates;      

% validate data
if length(sSig.eeg) ~= length(sSig.emg)
    error('eeg and emg must be the same length')
end

labels = ss.labels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep signals (similar pipeline to AccuSleep_classify)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calibration 
calData = ss.info.calibrationData;

% scale the emg
sSig.emg_rms = (sSig.emg_rms - calData(end, 1)) ./ calData(end, 2);
sSig.emg_rms = (sSig.emg_rms + 4.5) ./ 9;

% clip
sSig.emg_rms(sSig.emg_rms < 0) = 0;
sSig.emg_rms(sSig.emg_rms > 1) = 1;

% -------------------------------------------------------------------------
% normalize spectrogram for power in specific bands
specDisp = sSig.spec;
taxis = sSig.spec_tstamps;

% frequency indices
[~, f1idx] = min(abs(sSig.spec_freq - 1));
[~, f4idx] = min(abs(sSig.spec_freq - 4));
[~, f6idx] = min(abs(sSig.spec_freq - 6));
[~, f12idx] = min(abs(sSig.spec_freq - 12));
[~, f20idx] = min(abs(sSig.spec_freq - 20)); % index in f of 20Hz
[~, f50idx] = min(abs(sSig.spec_freq - 50)); % index in f of 50Hz

% select frequencies up to 50 Hz, and downsample between 20 and 50 Hz
specNorm = specDisp(:, [1 : (f20idx - 1), f20idx : 2 : f50idx]);
% take log
specNorm = log(specNorm);
% scale the spectrogram
for j = 1:size(specNorm, 2)
    specNorm(:, j) = (specNorm(:, j) - calData(j, 1)) ./ calData(j, 2);
    specNorm(:, j) = (specNorm(:, j) + 4.5) ./ 9; % clip z scores
end
% clip
specNorm(specNorm < 0) = 0;
specNorm(specNorm > 1) = 1;

sDelta = sum(specNorm(:, f1idx : f4idx), 2);
sTheta = sum(specNorm(:, f6idx : f12idx), 2);
sRatio = sDelta ./ sTheta;

% calculate psd according to states
[psd, faxis] = calc_psd('sig', sSig.eeg,...
    'fs', sSig.fs, 'graphics', false, 'winCalc', ss.stateEpochs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(false)
set(groot, 'DefaultAxesFontSize', 12)
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
for istate = 1 : nstates
    stateLabels = find(labels == istate);
    scatter(stateLabels / epochLen / 60 / 60,...
        sSig.emg_rms(stateLabels),...
        3, cfg.colors{istate})
end
axis tight
ylim([min(sSig.emg_rms) 1])
ylabel('Norm. emg RMS')
set(gca, 'XTick', [], 'YTick', [])

% spectrogram
fh.CurrentAxes = sb2;
% time axis in hours
nbins = length(taxis);
tSpec = ((1 : nbins) * epochLen - epochLen / 2) / 3600; % spectrogram time axis, in seconds
showFreqs = find(faxis <= 15);  % choose freqs to display
% take a sample of the spectrogram to help initialize the colormap
sampleBins = randperm(nbins, round(nbins / 10));
specSample = reshape(specDisp(sampleBins, showFreqs), 1, length(sampleBins) * length(showFreqs));
caxis1 = prctile(specSample, [6 98]);
% plot 
imagesc(tSpec, faxis(showFreqs), specDisp(:, showFreqs)', caxis1);
colormap(AccuSleep_colormap());
axis('xy')
ylabel('Frequency [Hz]')
xlabel('Time [h]')
% set(gca, 'YTick', [])

% spectral power per state
fh.CurrentAxes = sb3;
ph = plot(faxis, psd ./ sum(psd, 2), 'LineWidth', 3);
set(ph, {'color'}, cfg.colors(sstates))
xlim([0 30])
xlabel('Frequency [Hz]')
ylabel('Norm PSD')     
% set(gca, 'YTick', [])

% scatter spectrogram vs. emg (could not use axes handle w/ gscatter)
fh.CurrentAxes = sb4;
hold on
for istate = sstates
scatter(sRatio(labels == istate), sSig.emg_rms(labels == istate)',...
    2, cfg.colors{istate}, 'filled')
end
ylabel('Norm. emg RMS')
xlabel('Delta / Theta Ratio')
set(gca, 'YTick', [])

% histogram of epoch lengths for NREM, REM, and WAKE
fh.CurrentAxes = sb5;
hold on
epMat = cell2nanmat(ss.epLen(sstates));
plot([1 : size(epMat, 2)], mean(epMat, 1, 'omitnan'),...
    'kd', 'markerfacecolor', 'k')
boxplot(epMat, 'PlotStyle', 'traditional', 'Whisker', 6);
bh = findobj(sb5, 'Tag', 'Box');
bh = flipud(bh);
for ibox = 1 : length(bh)
    patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
        cfg.colors{sstates(ibox)}, 'FaceAlpha', 0.5)
end
xticklabels(cfg.names(sstates))
xtickangle(45)
set(sb5, 'YScale', 'log')
ylabel('Epoch Length [log(s)]')
ylim([0 ceil(prctile(epMat(:), 99.99))])
    
% percent time in state
fh.CurrentAxes = sb6;
pie(sum(cell2nanmat(ss.epLen(sstates), 2), 1, 'omitnan'), ones(1, length(sstates)));
hold on
ph = findobj(sb6, 'Type', 'Patch');
set(ph, {'FaceColor'}, flipud(cfg.colors(sstates)))
legend(cfg.names, 'Units', 'normalized', 'Position', [0.81 0.11 0.10 0.20]);

% save figure
if saveFig
    figpath = fullfile('graphics', 'sleepState');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_stateSeparation', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF
