function [emg, spec, boutTimes] = recRep_prep(varargin)

% loads and prepares data for plotting a recording representative (see
% recRep_plot). this was separated so that the figure can be plotted
% without reloading the data each time. the data loaded here does not need
% to be updated with changing the time windows
%
% INPUT:
%   basepath    char. path to session folder {pwd}
%   panels2plot string array. determines which plots to include and in what
%               order. can be "spec", "emg", "hypnogram", "raster", "raw"
%   fs_eeg      numeric. sampling frequency for downsampling eeg
%
% OUTPUT
%
% CALLS
%
% TO DO LIST
%
% 05 mar 24 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'panels2plot', []);
addParameter(p, 'fs_eeg', []);

parse(p, varargin{:})
basepath        = p.Results.basepath;
panels2plot     = p.Results.panels2plot;
fs_eeg          = p.Results.fs_eeg;



grph_save('fh', gcf, 'fname', 'Representative_lh96_Emg', 'frmt', {'ai', 'jpg', 'fig'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(basepath)
[~, basename] = fileparts(basepath);

% load data
varsFile = ["session"; "units"; "sleep_statesEmg"];
v = basepaths2vars('basepaths', {basepath}, 'vars', varsFile);
if isfield(v, 'sleep_statesEmg'), [v.ssEmg] = v.sleep_statesEmg; v = rmfield(v, 'sleep_statesEmg'); end
if isfield(v, 'ssEmg'), v.ss = v.ssEmg; end

% sleep signals
filename = fullfile(basepath, [basename, '.sleep_sig.mat']);
load(filename, 'emg_rms');
emg = movmean(emg_rms, 33);
load(filename, 'eeg');

% re-calc state epochs for better visualization.
% LSLEEP => NREM; N/REM => REM
minDur = [10, 10, 5, 4, 4, 4];
% minDur = 4;
interDur = 20;
labels = v.ss.labels;
labels(labels == 6) = 5;
labels(labels == 3) = 4;
bouts = as_bouts('labels', labels,...
    'minDur', minDur, 'interDur', interDur);
boutTimes = bouts.times;

% subsample eeg and re-calc spec
eeg = eeg(1 : 1250 / fs_eeg : end);
spec = calc_spec('sig', eeg, 'fs', fs_eeg, 'graphics', false, 'saveVar', false,...
    'padfft', -1, 'winstep', 1, 'logfreq', true, 'ftarget', [],...
    'ch', [{1}], 'force', true);

end


% EOF
