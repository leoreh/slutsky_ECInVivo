function [labels, netScores] = AccuSleep_classify(spec, freq, emg_rms, calibrationData, net)
% AccuSleep_classify  Classify brain states
% Zeke Barger 021321
%
%   Arguments: (note that these are data, not paths to files)
%   sSig.eeg - a 1-D matrix of sSig.eeg data
%   sSig.emg - a 1-D matrix of sSig.emg data
%   net - a trained network. To create one, use AccuSleep_train
%   SR - sampling rate (Hz)
%   boutLen - length of each bout (sec). THIS MUST MATCH THE EPOCH LENGTH
%              THAT WAS USED WHEN TRAINING THE NETWORK
%   calibrationData - matrix specifying how to scale features of the
%       sSig.eeg/sSig.emg. It is recommended to create one of these for each
%       combination of mouse and recording equipment. This is the output of
%       createCalibrationData.m
%   minBoutLen (optional) - minimum length (sec) of a bout of any brain state

% load config data
cfg = as_loadConfig();
nstates = cfg.nstates;
boutLen = cfg.boutLen;
minBoutLen = cfg.minBoutLen;
netScores = [];

% select frequencies up to 50 Hz, and downsample between 20 and 50 Hz.
% redundant after changes to as_prepSig, May2022
% [~,f20idx] = min(abs(freq - 20)); % index in f of 20Hz
% [~,f50idx] = min(abs(freq - 50)); % index in f of 50Hz
% spec = spec(:, [1:(f20idx-1), f20idx:2:f50idx]);

% if the spectrogram isn't the same height as the network, that's a problem
if size(spec,2) ~= (net.Layers(1,1).InputSize(2) - 9)
   disp('Error: frequency axes for network and data do not match')
   labels = [];
   return
end

% take log of the spectrogram
spec = log(spec);

% scale the spectrogram
for j = 1:size(spec,2)
    spec(:,j) = (spec(:,j) - calibrationData(j,1)) ./ calibrationData(j,2);
    spec(:,j) = (spec(:,j) + 4.5)./9; % clip z scores
end
% scale the sSig.emg
emg_rms = (emg_rms - calibrationData(end,1)) ./ calibrationData(end,2);
emg_rms = (emg_rms + 4.5)./9;

% clip them
emg_rms(emg_rms < 0)=0;
emg_rms(emg_rms > 1)=1;
spec(spec < 0)=0;
spec(spec > 1)=1;

% find how much time on either side of central timepoint to include
% this is based on the size of the trained network.
pad = round((net.Layers(1,1).InputSize(1) - 1)/2);

% pad the spectrogram and sSig.emg
spec = [repmat(spec(1,:), pad, 1); spec; repmat(spec(end,:), pad, 1)];
emg_rms = [repmat(emg_rms(1), 1, pad), emg_rms, repmat(emg_rms(end), 1, pad)];

% preallocate the image stack
X = zeros((pad*2+1),net.Layers(1,1).InputSize(2),1,length(emg_rms)-pad*2);

% create an image for each time step
for i = (pad+1):(length(emg_rms)-pad)
    X(:,:,1,(i-pad)) = [spec((i-pad):(i+pad),:), repmat(emg_rms((i-pad):(i+pad))',1,9)];
end

% classify
X = uint8(X.*255);
[labels, netScores] = classify(net, X);
labels = double(labels);

% put labels in correct orientation
if isrow(labels)
    labels = labels';
end

% remove bouts that are too short
if minBoutLen > boutLen
    labels = enforceMinDuration(labels, ones(1,nstates) * ceil(minBoutLen / boutLen), ...
        1:nstates, 0);
end
