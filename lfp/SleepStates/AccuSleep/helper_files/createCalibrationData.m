function [calibrationData] = createCalibrationData(spec, freq, emg_rms, labels)
% CREATECALIBRATIONDATA  Calculate parameters for mixture z-scoring
% Zeke Barger 021321
%
%   Background:
%   Before classifying EEG/EMG data from a given subject recorded on a given
%   recording setup, it is necessary to scale the features of the data so
%   that they resemble the classifier's training data. See the accompanying
%   paper for details.
%
%   Output:
%   calibrationData - a set of parameters for mixture z-scoring. These are
%   required by AccuSleep_GUI or AccuSleep_classify when classifying sleep
%   stages automatically.

% load config data
cfg = as_loadConfig();
nstates = cfg.nstates;
boutLen = cfg.boutLen;

% set the fixed mixture weights
weights = cfg.weights(:); % rem, wake, nrem, ...

% check if there are at least a few labeled bouts for each state
count_by_state = zeros(1, nstates);
for i = 1:nstates
    count_by_state(i) = sum(labels==i);
end
if any(count_by_state < 3)
    calibrationData = [];
    disp('At least a few bouts of each state must be labeled');
    return
end

% select frequencies up to 50 Hz, and downsample between 20 and 50 Hz
% [~,f20idx] = min(abs(freq - 20)); % index in f of 20Hz
% [~,f50idx] = min(abs(freq - 50)); % index in f of 50Hz
% spec = spec(:, [1:(f20idx-1), f20idx:2:f50idx]);
% take log of the spectrogram
spec = log(spec);

% check if labels has the correct length
diffLD = length(freq) - length(labels);
if diffLD > 0
    if diffLD > 2
        calibrationData = [];
        disp('Labels are not the proper length for this recording');
        return
    else
        spec = spec(1 : length(labels), :);
        emg_rms = emg_rms(1 : length(labels));
    end
end

% make an image for the entire recording
spec = [spec, emg_rms']; % we only need one of the sSig.emg columns

m = zeros(size(spec, 2), nstates);
v = zeros(size(spec, 2), nstates);
% for each sleep stage
for i = 1:nstates
    % calculate mean and var for all features
    m(:,i) = mean(spec(labels==i,:));
    v(:,i) = var(spec(labels==i,:));
end

% mixture means are just weighted combinations of state means
calibrationData = zeros(size(spec, 2), 2);
calibrationData(:,1) = m * weights;

% mixture variance is given by law of total variance
% sqrt to get the standard deviation
calibrationData(:,2) = sqrt(v * weights +...
    ((m - repmat(calibrationData(:,1),1,nstates)).^2) * weights);
