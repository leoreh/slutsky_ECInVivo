function [EMG, EEG, info] = as_cleanSig(eegOrig, emgOrig, lfpfile, varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean recording from artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% two alternatives - semi-automatically or manually (see below)
% -------------------------------------------------------------------------
% a different approach is to mark artifacts during the calibration process
% as undefined. than, must insert the calibration labels in place of the
% final results or use the gui and select the option only overwrite
% undefined epochs. however, this will not prevent the artifacts from
% influincing the spectrogram if it is normalized.
% -------------------------------------------------------------------------
% currently uses a temporal resolution that corresponds to epoch length but
% should be improved to < 1 s. just need to find a way to restore deleted
% sections afterwards.

cleanRes = 2.5;       % size of bins for cleaning the data [s]
artifactThr = 2;      % threshold for removing sig in alt 2 [z scores]
if exist(artifactsfile) && ~forceLoad
    load(artifactsfile)
else
    if cleanRec == 1
        % ALT 1: manual mark bad times using AccuSleep gui (use NREM
        % as a bad epoch). when done save labels.
        AccuSleep_viewer(EEG, EMG, fs, cleanRes, [], artifactsfile)
        uiwait
        load(artifactsfile)
        ignoretimes(labels == 3) = 1;
        ignoretimes(labels == 4) = 0;
        
        % ALT 2: semi-automatically find bad epochs from spectrogram
    elseif cleanRec == 2
        [s, ~, ~] = createSpectrogram(standardizefs(EEG, fs, 128), 128, cleanRes);
        ignoretimes = zeros(floor(recDur / cleanRes), 1);
        ignoretimes(zscore(mean(zscore(s), 2)) > artifactThr) = 1;
    else
        ignoretimes = [];
    end
end
save(artifactsfile, 'ignoretimes')

% interpolate badtimes to length of the signals
if ~isempty(ignoretimes)
    rmtimes = [interp1([1 : length(ignoretimes)] * cleanRes, ignoretimes, tstamps_sig,...
        'linear', 'extrap')]';
else
    rmtimes = zeros(length(EEG), 1);
end
rmtimes = round(rmtimes);
EEG = EEG(~rmtimes);
EMG = EMG(~rmtimes);
tRemoved = sum(rmtimes) / fs; % total time removed from data [s]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restore bad times to labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tstamps_labels = [0 : epochLen : recDur];
% labelsLen = floor(length(EMG_newFs) / SR / epochLen);
% x = [interp1([1 : length(badtimes)] / cleanRes, badtimes,...
%     [1 : labelsLen] * epochLen, 'linear', 'extrap')]';
% x = round(x);

ignoreEpochs = binary2epochs('vec', ignoretimes, 'minDur', [], 'maxDur', [],...
    'interDur', [], 'exclude', false); 
for i = 1 : size(ignoreEpochs, 1)
    labels = [labels(1 : ignoreEpochs(i, 1) - 1);...
        ones(diff(ignoreEpochs(i, :)), 1) * 4;...
    labels(ignoreEpochs(i, 1) : end)];

    labels_calibration = [labels_calibration(1 : ignoreEpochs(i, 1) - 1);...
        ones(diff(ignoreEpochs(i, :)), 1) * 4;...
    labels_calibration(ignoreEpochs(i, 1) : end)];
end


end