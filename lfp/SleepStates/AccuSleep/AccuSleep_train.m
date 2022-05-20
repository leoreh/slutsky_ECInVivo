function [net, netInfo] = AccuSleep_train(basepaths, epochs, imageLocation)
% AccuSleep_train  Train a network for classifying brain states
% Zeke Barger, 021321
%
%   Arguments:
%   basepaths - a cell array of strings depicting the full path to the data
%   folder. each folder should contain a [basename].sleep_sig and
%   [basename].labelsMan.
%   epochs - the number of epochs for the network to consider at once when
%       scoring each individual epoch. More epochs provide more context,
%       but using fewer epochs is more efficient.
%       ***Must be an odd number greater than or equal to 9***
%       For 2.5s epochs, a reasonable value is 13.
%   imageLocation (optional) - training images will be stored here. If
%       no location is specified, a temporary folder will be created in the
%       current directory and deleted after the network is trained.   
%
%   Output:
%   net - a trained network (SeriesNetwork object) that can be used by
%   AccuSleep_GUI or AccuSleep_classify for classification.
%   Note - Unfortunately, it seems that networks trained using MATLAB 
%   2019a or later are not readable by versions 2018b or earlier. 
%   However, networks trained using MATLAB 2018b or earlier seem to be 
%   forward compatible.


%% Check the inputs and prepare the image folder
net = [];
switch nargin
    case {0, 1}
        error('Not enough arguments')
    case 2
        imageLocation = [char(cd),filesep,'training_images_',...    
            char(datetime(now,'ConvertFrom','datenum',...
            'Format','yyyy-MM-dd_HH-mm-ss'))];
        deleteImages = 1;
    case 3
        deleteImages = 0;
end

% load config data
cfg = as_loadConfig();
nstates = cfg.nstates;
epochLen = cfg.epochLen;

% make directory to hold training images
% remove slash at the end of the path if it's there already
if strcmp(imageLocation(end),'\') || strcmp(imageLocation(end),'/')
    imageLocation = imageLocation(1:end-1);
end
% put a folder inside the target location
imageLocation = [imageLocation,filesep,'training_images_',...    
            char(datetime(now,'ConvertFrom','datenum',...
            'Format','yyyy-MM-dd_HH-mm-ss'))];
mkdir(imageLocation);
% make folders for each class
for i = 1:nstates
    mkdir([imageLocation,filesep,num2str(i)])
end

% calculate how many epochs on either side of central epoch to include in
% the images
pad = round((epochs - 1) / 2);

%%  Process the recordings
nFiles = size(basepaths, 1);
disp(['Processing ',num2str(nFiles),' files (this may take a while)'])
for i = 1:nFiles
    % load the files
    basepath = basepaths{i};
    [~, basename] = fileparts(basepath);
    sigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
    labelmanfile = fullfile(basepath, [basename, '.sleep_labelsMan.mat']);
    
    SR = load(sigfile, 'fs');
    data = struct;
%     data.a = load(basepaths{i,1});
%     data.b = load(basepaths{i,2});
    data.c = load(labelmanfile, 'labels');
    s = load(sigfile, 'spec');
    s = s.spec;
    f = load(sigfile, 'spec_freq');
    f = f.spec_freq;
    processedEMG = load(sigfile, 'emg_rms');
    processedEMG = processedEMG.emg_rms;
    
    % make sure the files look ok
%     switch checkFiles(data, SR, epochLen, i, nstates)
%         % user opts to skip
%         case 0
%             continue
%             % user opts to quit
%         case -1
%             % delete the images
%             if deleteImages
%                 rmdir(imageLocation);
%             end
%             return
%     end
    
    % create the spectrogram
%     [s, ~, f] = createSpectrogram(data.a.EEG, SR, epochLen);
    % select frequencies up to 50 Hz, and downsample between 20 and 50 Hz.
    % redundant after changes to as_prepSig, May2022
%     [~,f20idx] = min(abs(f - 20)); % index in f of 20Hz
%     [~,f50idx] = min(abs(f - 50)); % index in f of 50Hz
%     s = s(:, [1:(f20idx-1), f20idx:2:f50idx]);
    % take log of the spectrogram
    s = log(s);
    % calculate log rms EMG for each epoch
%     processedEMG = processEMG(data.b.EMG, SR, epochLen);
    % make sure labels are the right length
    if length(data.c.labels) > length(processedEMG)
        data.c.labels = data.c.labels(1:length(processedEMG));
    end
    if length(data.c.labels) < length(processedEMG)
        processedEMG = processedEMG(1:length(data.c.labels));
        s = s(1:length(data.c.labels),:);
    end
    
    % create an image for the entire recording
    im = [s, repmat(processedEMG',1,9)];
    
    % calibrate the features
    weights = cfg.weights(:); % rem wake nrem ...
    m = zeros(size(im, 2), nstates);
    v = zeros(size(im, 2), nstates);
    % for each brain state
    for j = 1:nstates
        % calculate mean and var for all features
        m(:,j) = mean(im(data.c.labels==j,:));
        v(:,j) = var(im(data.c.labels==j,:));
    end
    % mixture means are just weighted combinations of state means
    calibrationData = zeros(size(im, 2), 2);
    calibrationData(:,1) = m * weights;
    % mixture variance is given by law of total variance
    % sqrt to get the standard deviation
    calibrationData(:,2) = sqrt(v * weights +...
        ((m - repmat(calibrationData(:,1),1,nstates)).^2) * weights);
    
    % perform scaling
    for j = 1:size(im,2)
        im(:,j) = (im(:,j) - calibrationData(j,1)) ./ calibrationData(j,2);
        im(:,j) = (im(:,j) + 4.5)./9; % clip z scores
    end
    
    % clip the image
    im(im < 0)=0;
    im(im > 1)=1;
    % generate the images
    disp(['Creating images for recording ',num2str(i)])
    % take all datapoints with sufficient epochs on either side
    for j = (pad+1):(length(processedEMG)-pad)
        % only take timepoints with labels in range 1:n_states
        if data.c.labels(j) > 0 && data.c.labels(j) < (nstates+1)
            thisImg = im((j-pad):(j+pad),:);
            % create filename
            fName = [imageLocation,filesep,num2str(data.c.labels(j)),...
                filesep,'rec',num2str(i),'t',num2str(j-pad,'%05.f'),'.png'];
            % write to file 
            imwrite(thisImg, fName);
        end
    end
end

%% Oversample training data to achieve balanced classes
disp('Balancing classes')
% count examples of each class
counts = zeros(1, nstates);
for i = 1:nstates
    d = dir([imageLocation,filesep,num2str(i)]);
    counts(i) = length(d)-2;
end

% make more if necessary
for i = 1:nstates
   if counts(i) < max(counts)
      numSamples = max(counts) - counts(i);
      d = dir([imageLocation,filesep,num2str(i)]);
      for j = 1:numSamples
          k = ceil(rand*counts(i))+2;
           copyfile([d(k).folder,filesep,d(k).name], [d(k).folder,filesep,num2str(j),'__',d(k).name,'.png']);
          
      end
   end
end

%% Train the network
% create image datastore
imds = imageDatastore(imageLocation, 'IncludeSubfolders',true,...
    'LabelSource','foldernames');

% get image size
img = readimage(imds,1);
imsize = size(img);
% get number of images per label
labelCount = countEachLabel(imds);
% and decide how much data to use for training
trainFraction = 0.8;
% also decide how often to perform validation
vf = .012*labelCount.Count(1) - 6.4;
% make sure it's at least 30 but less than 200
vf = round(max([min([vf, 200]), 30]));

% split images into training and validation sets with 80%  used for training
[imdsTrain,imdsValidation] = splitEachLabel(imds,trainFraction,'randomize');

% set training options
options = trainingOptions('sgdm', ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.85, ...
    'LearnRateDropPeriod',1, ...
    'InitialLearnRate',0.015, ...
    'MaxEpochs',10, ...
    'Shuffle','every-epoch', ...
    'ValidationData',imdsValidation, ...
    'ValidationFrequency',vf, ...
    'ValidationPatience',15,...
    'Verbose',false, ...
    'ExecutionEnvironment','auto',...
    'MiniBatchSize',256,...
    'Plots','none');

layers = [
    imageInputLayer([imsize(1) imsize(2) 1]);%,'AverageImage',ai)
    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(nstates)
    softmaxLayer
    classificationLayer];

% train
disp('Training network')
[net, netInfo.train] = trainNetwork(imdsTrain,layers,options);

disp('Training complete: Final validation accuracy:')
disp([num2str(netInfo.train.ValidationAccuracy(end)),'%'])

% delete the images
if deleteImages
    rmdir(imageLocation,'s');
end

%% Other functions

% make sure each set of eeg, emg, and label files looks correct
% offer the option to skip or quit on encountering a problem
function [looksGood] = checkFiles(data, SR, epoch_length, i, n_states)
looksGood = 1;

% need to have all the correct variables
if ~isfield(data.a,'EEG') || ~isfield(data.b,'EMG') || ~isfield(data.c,'labels')
    looksGood = showError(['Recording ',num2str(i),' has improperly formatted files. ']);
    return
end
% eeg and emg must be the same length (we'll trust that the content is correct)
if length(data.a.EEG)~=length(data.b.EMG)
    looksGood = showError(['Recording ',num2str(i),' has EEG and EMG files of different lengths. ']);
    return
end
% labels are numeric
if ~isa(data.c.labels,'numeric')
    looksGood = showError(['Recording ',num2str(i),' has non-numeric labels. ']);
    return
end
% should be a vector
if length(size(data.c.labels)) > 2 || min(size(data.c.labels)) >  1
    looksGood = showError(['Recording ',num2str(i),' has labels that are not a vector. ']);
    return
end
% in the range 1:n_states (at least some of them)
if ~any(data.c.labels > 0 & data.c.labels < (n_states+1))
    looksGood = showError(['Recording ',num2str(i),' has no labels in the range 1:',...
        num2str(n_states),'. ']);
    return
end
% labels must be approximately the same length as EEG / EMG
if abs(length(data.c.labels)*epoch_length*SR - length(data.a.EEG)) / length(data.a.EEG) > 0.05
    looksGood = showError(['Recording ',num2str(i),' has EEG and label files of different lengths. ']);
    return
end

% display an error message, and interpret the response
function [q] = showError(message)
q = 0;
answer = questdlg([message,'Skip this recording, or quit?'], 'Error','Skip','Quit','Skip');
if strcmp(answer,'Quit')
    q = -1;
end
