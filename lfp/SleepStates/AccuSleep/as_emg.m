function ssEmg = as_emg(varargin)

% recieves (or loads) an emg_rms signal (from sSig) and creates a labelsOrig
% vector from the bimodal distribution. Then creates bouts and saves
% in a new ss struct. Can manually review the labels in accusleep.
%
% INPUT:
%   basepath        char. fullpath to recording folder {pwd}
%   minDur          numeric. minimum duration of a bout. 
%                   if length(minDur) == nstates than a different minimum
%                   duration will be applied to each state
%   interDur        numeric. combine bouts separated by <= interDur
%   saveVar         logical. save ss var {true}
%   flgInspct       logical. manually review classification 
% 
% OUTPUT
%
% DEPENDENCIES
%   as_bouts
% 
% TO DO LIST
%
% 25 dec 24 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addOptional(p, 'minDur', 5, @isnumeric);
addOptional(p, 'interDur', 4, @isnumeric);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'flgInspct', false, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
minDur          = p.Results.minDur;
interDur        = p.Results.interDur;
saveVar         = p.Results.saveVar;
flgInspct       = p.Results.flgInspct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
[~, basename] = fileparts(basepath);
sigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
statesfile = [basename '.sleep_statesEmg.mat'];

% state config
cfg = as_loadConfig('flgEmg', true);
snames = cfg.names;
clr = cfg.colors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data and create labelsOrig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% emg data
load(sigfile, 'emg_rms');

% find threshold to separate the bimodal distribution of emg
[~, cents] = kmeans(emg_rms(:), 2);
emgThr = mean(cents);

% create EMG labelsOrig
labels = double(emg_rms > emgThr);
labels(emg_rms < emgThr) = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create state bouts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bouts = as_bouts('labels', labels,...
    'minDur', minDur, 'interDur', interDur, 'flgOtl', true,...
    'sstates', [1, 2], 'graphics', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create and save struct 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ssEmg.info.names = cfg.names;
ssEmg.info.colors = cfg.colors;
ssEmg.info.boutLen = cfg.boutLen;
ssEmg.info.minBoutLen = cfg.minBoutLen;
ssEmg.info.emgThr = emgThr;
ssEmg.info.minDur = minDur;
ssEmg.info.interDur = interDur;
ssEmg.info.runtime = datetime("now");
ssEmg.labels = bouts.labels;
ssEmg.bouts = bouts;

if saveVar
    save(statesfile, 'ssEmg')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manual inspection 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flgInspct
    sSig = load(sigfile);
    fname_labels = [basename, '.sleep_labelsManEmg.mat'];
    
    labels = bouts.labelsOrig;
    AccuSleep_viewer(sSig, labels, fname_labels, cfg)
 
end


end

% EOF