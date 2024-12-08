function [stateEpochs, epochStats] = as_epochs(varargin)

% recieves a labels vector and returns the state epochs. basically a
% wrapper for binary2epochs with params adjusted for sleep scoring. can
% separate the recording to timebins and calculate epochs for each time
% bin, can apply duration thresholds, and remove spec outliers. if flgEmg,
% will create a labels vector from the bimodal distribution of emg rms
%
% INPUT:
%   labels          numeric. see as_classify 
%   minDur          numeric. minimum duration of an epoch. 
%                   if length(minDur) == nstates than a different minimum
%                   duration will be applied to each state
%   interDur        numeric. combine epochs separated by <= interDur
%   sstates         numeric. values of labels for which epochs will be
%                   calculated. default is 1 : nstates
%   timebins        n x 2 numeric of time windows (indices to labels). 
%                   epochs of each state will be calculated for each time
%                   window separately. 
%   nbins           numeric. if timebins not specified, will divide labels
%                   to nbins timebins
%   confMarg        2 x 1 numeric. confidance margin. for example, [1 2]
%                   will change the epoch [789 810] to [790 808]. 
%                   this is to assure that the previous and next state do
%                   not influence the current state
%   rmOtl           logical. remove spectrogram outliers. see get_specOutliers.m
%   emgFlg          logical
%   graphics        logical
% 
% OUTPUT
%
% DEPENDENCIES
% 
% TO DO LIST
%
% 12 jan 22 LH      updates:
% 26 mar 24             rmOtl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'labels', [], @isnumeric);
addOptional(p, 'minDur', [10, 5, 5, 10, 5, 5], @isnumeric);
addOptional(p, 'interDur', 4, @isnumeric);
addOptional(p, 'sstates', [], @isnumeric);
addOptional(p, 'timebins', [], @isnumeric);
addOptional(p, 'nbins', 1, @isnumeric);
addOptional(p, 'confMarg', [], @isnumeric);
addOptional(p, 'rmOtl', false, @islogical);
addOptional(p, 'flgEmg', false, @islogical);
addOptional(p, 'graphics', false, @islogical);

parse(p, varargin{:})
labels          = p.Results.labels;
minDur          = p.Results.minDur;
interDur        = p.Results.interDur;
sstates         = p.Results.sstates;
timebins        = p.Results.timebins;
nbins           = p.Results.nbins;
confMarg        = p.Results.confMarg;
rmOtl           = p.Results.rmOtl;
flgEmg          = p.Results.flgEmg;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file
basepath = pwd;
[~, basename] = fileparts(pwd);
sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

% state params
if flgEmg
    namePrefix = 'epochsEmg';
    clr = {[150, 70, 55] / 255, [200, 170, 100] / 255};
    snames = {'High-EMG', 'Low-EMG'};
    sstates = [1, 2];

    % emg data
    emg = load(sleepfile, 'emg_rms');
    emg = emg.emg_rms;

    % find threshold to separate the bimodal distribution of emg
    [~, cents] = kmeans(emg(:), 2);
    emgThr = mean(cents);

    % create EMG labels
    labels = double(emg > emgThr);
    labels(emg < emgThr) = 2;

else
    
    cfg = as_loadConfig();
    clr = cfg.colors(sstates);
    snames = cfg.names(sstates);

end

% selected states
if isempty(sstates)
    nstates = cfg.nstates;
    sstates = 1 : nstates;
else
    nstates = length(sstates);
end

% validate duration thresholds
if length(minDur) == 1
    minDur = repmat(minDur, nstates, 1);
elseif length(minDur) ~= nstates
    error('minDur length is different than the number of states')
end
if length(interDur) == 1
    interDur = repmat(interDur, nstates, 1);
elseif length(interDur) ~= nstates
    error('interDur length is different than the number of states')
end

% timebins 
if isempty(timebins)
    timebins = n2chunks('n', length(labels), 'nchunks', nbins);
end
nbins = size(timebins, 1);
binLen = diff(timebins');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create state epochs
for ibin = 1 : nbins
    
    labels_bin = labels(timebins(ibin, 1) : timebins(ibin, 2));
    
    for istate = 1 : nstates
        
        binaryVec = zeros(length(labels_bin), 1);
        binaryVec(labels_bin == sstates(istate)) = 1;

        stateEpochs{ibin, istate} = binary2epochs('vec', binaryVec,...
            'minDur', minDur(istate), 'maxDur', [],...
            'interDur', interDur(istate), 'exclude', false, 'printFlag', false);

        stateEpochs{ibin, istate} = stateEpochs{ibin, istate} + timebins(ibin, 1) - 1;
    end
end

% remove spec outliers
origEpochs = stateEpochs;
otl = [];
if rmOtl
    
    % get outliers
    otl = get_specOutliers('basepath', pwd, 'saveVar', true,...
        'flgCalc', false, 'flgForce', false, 'graphics', false);
    
    % subtract outlier epochs from state epochs 
    for ibin = 1 : nbins
        for istate = 1 : length(origEpochs)
            if ~isempty(origEpochs{istate})
                stateEpochs{ibin, istate} = SubtractIntervals(origEpochs{istate}, otl.epochs);

                % re-apply min duration criterion
                epLen = stateEpochs{ibin, istate}(:, 2) - stateEpochs{ibin, istate}(:, 1);
                badIdx = epLen < minDur(istate);
                stateEpochs{ibin, istate}(badIdx, :) = [];
            end
        end
    end

end

% remove confindance margins
if ~isempty(confMarg)
    funh = @(x) [x(:, 1) + confMarg(1) x(:, 2) - confMarg(2)];
    stateEpochs = cellfun(funh, stateEpochs, 'uni', false);
end

% epoch stats
epochStats.epLen = cellfun(@(x) (diff(x')'), stateEpochs, 'uni', false);
epochStats.nepochs = cellfun(@length, epochStats.epLen);
epochStats.totDur = cellfun(@sum, epochStats.epLen);
epochStats.prctDur = epochStats.totDur * 100 ./ binLen';
epochStats.confMarg = confMarg;
epochStats.otl = otl;
epochStats.origEpochs = origEpochs;
epochStats.timebins = timebins;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~graphics
    return
end

% open figure
fh = figure;
setMatlabGraphics(true)
set(fh, 'WindowState','maximized');
tlayout = [4, nstates];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, basename, 'interpreter', 'none')

% hypnogram
axh = nexttile(th, 1, [1, nstates]); hold on; cla
for istate = 1 : nstates
    tmpEpochs{istate} = cell2padmat(stateEpochs(:, istate), 1);
end
plot_hypnogram('stateEpochs', tmpEpochs, 'clr', clr, 'axh', axh,...
    'sstates', [1 : nstates])
axis tight
yLimit = ylim;
if size(timebins, 1) > 1
    plot([timebins(2 : end, 1), timebins(2 : end, 1)], yLimit, '--k')
end
xval = 3600 : 3600 : length(labels);
xticks(xval)
xticklabels(string(xval / 3600))
xlabel('Time (hr)')

% percent duration
tilebias = nstates;
for istate = 1 : nstates
    axh = nexttile(th, istate + tilebias, [1, 1]); cla
    dataMat = squeeze(epochStats.prctDur(:, istate))';
    plot_boxMean('dataMat', dataMat, 'clr', clr{istate},...
        'plotType', 'bar', 'axh', axh, 'allpnts', true)
    ylabel('Duration (%)')
end

% no. epochs
tilebias = nstates * 2;
for istate = 1 : nstates
    axh = nexttile(th, istate + tilebias, [1, 1]); cla
    dataMat = squeeze(epochStats.nepochs(:, istate, :))';
    plot_boxMean('dataMat', dataMat, 'clr', clr{istate},...
        'plotType', 'bar', 'axh', axh, 'allpnts', true)
    ylabel('No. Epochs')
end

% epoch length
tilebias = nstates * 3;
for istate = 1 : nstates
    axh = nexttile(th, istate + tilebias, [1, 1]); cla
    dataMat = cell2padmat(epochStats.epLen(:, istate, :), 2);
    plot_boxMean('dataMat', dataMat, 'clr', clr{istate},...
        'plotType', 'bar', 'axh', axh, 'allpnts', false)
    ylabel('Mean Epoch Length (s)')
end


end

% EOF