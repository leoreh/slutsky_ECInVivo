function bands = frStates_sessions(mname, varargin)

% arranges 
%
% INPUT:
%   mname           string. mouse name to analyze
%   idxBsl          numeric. indices to baseline sessions for normalization
%   flgNormTime     logical. normalize each session to baseline (1st
%                   session) {true}. only applied to bands
%   flgNormBand     logical. normalize psd to broadband {true}
%   flgEmg          logical. load [basename].psdEmg.mat {true}
%   flgAnalyze      logical. re-analyze psd {false}
%   saveVar         logical. save ss var {true}
%   graphics        logical. plot confusion chart and state separation {true}
%
% DEPENDENCIES:
%   calc_psd        
%
% TO DO LIST:
%
% 20 oct 22 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'flgEmg', true, @islogical);
addOptional(p, 'flgAnalyze', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
flgEmg          = p.Results.flgEmg;
flgAnalyze      = p.Results.flgAnalyze;
saveVar         = p.Results.saveVar;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load or analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
prct = 70;
flgEmg = true;
mname = 'lh96';
ntiles = 2;

% load session vars
varsFile = ["sleep_states"; "datInfo"; "session"; "psdEmg"; "spikes";...
    "fr"; "units"];
varsName = ["ss"; "datInfo"; "session"; "psd"; "spikes";...
    "fr"; "units"];
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName);
nfiles = length(basepaths);

% organize cell array of state ratio (nfiles x unitType x ntiles)
dataMat = cell(nfiles, 2, ntiles);
for ifile = 1 : nfiles
    
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

    % get window from psd
    if ~isempty(v(ifile).psd)
        wins = v(ifile).psd.info.wins;
    else
        wins = [0, Inf];
    end

    % load emg data
    emg = load(sleepfile, 'emg_rms');
    emg = emg.emg_rms;

    % get indices to high- and low-emg
    labels = double(emg > prctile(emg, prct));
    labels(emg < prctile(emg, 100 - prct)) = 2;

    % limit indices to time window and get "state" epochs
    labels = labels(wins(:, 1) : wins(:, 2));
    [stateEpochs, ~] = as_epochs('labels', labels, 'minDur', 10, 'interDur', 4);
    stateEpochs = stateEpochs([1 : 2]);

    % calc firing rate
    fr = calc_fr(v(ifile).spikes.times, 'basepath', basepath,...
        'graphics', false, 'binsize', 60, 'saveVar', false,...
        'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf],...
        'stateEpochs', stateEpochs);   
    stateRat = squeeze(fr.states.ratio(2, 1, :));

    % organize cell array of state ratio
    for iunit = 1 : 2
        units = v(ifile).units.clean(iunit, :);
        
        % get a cell array of indices according to mfr percentiles per unit
        % type
        [~, tileIdx] = vec2tileMat(fr.mfr, ntiles, units);

        % organize state ratio according to percentiles
        dataMat(ifile, iunit, :) = cellfun(@(x) stateRat(x), tileIdx, 'uni', false);
    end
end

% graphics and to prism
iunit = 1;
unitClr = {'b', 'r'};
unitType = {'RS', 'FS'};
stateNames = {'High-EMG', 'Low-EMG'};

fh = figure;
% fh.Position = [0.1 0.1 0.8 0.8];
th = tiledlayout(1, ntiles, 'TileSpacing', 'Compact');
for itile = 1 : ntiles
    axh = nexttile;
    ydata = cell2nanmat(dataMat(:, iunit, itile), 2);
    plot_boxMean('dataMat', ydata, 'clr', unitClr{iunit}, 'allPnts', true)
    hold on
        plot(xlim, [0, 0], '--k')
    ylabel({sprintf('%s - %s /', stateNames{2}, stateNames{1}),...
        sprintf('%s + %s', stateNames{2}, stateNames{1})})
    xlabel('Session')
    title(sprintf('Tile #%d', itile))

end
title(th, mname)