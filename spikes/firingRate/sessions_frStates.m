function frTiles = sessions_frStates(mname, varargin)

% arranges 
%
% INPUT:
%   mname           string. mouse name to analyze
%   ntiles          numeric. number of tiles to divide mfr
%   flgEmg          logical. load [basename].psdEmg.mat {true}
%   flgAnalyze      logical. re-analyze psd {false}
%   saveVar         logical. save ss var {true}
%   graphics        logical. plot confusion chart and state separation {true}
%
% DEPENDENCIES:
%
% TO DO LIST:
%
% 20 oct 22 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'ntiles', 2, @isnumeric);
addOptional(p, 'flgEmg', true, @islogical);
addOptional(p, 'flgAnalyze', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
ntiles          = p.Results.ntiles;
flgEmg          = p.Results.flgEmg;
flgAnalyze      = p.Results.flgAnalyze;
saveVar         = p.Results.saveVar;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and prepare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load session vars
if flgEmg
    varsFile = ["sleep_states"; "datInfo"; "session"; "psdEmg"; "spikes";...
        "fr"; "units"];
else
    varsFile = ["sleep_states"; "datInfo"; "session"; "psd"; "spikes";...
        "fr"; "units"];
end
varsName = ["ss"; "datInfo"; "session"; "psd"; "spikes";...
    "fr"; "units"];
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName);
nfiles = length(basepaths);

% get params for psd file or use defaults. this is for emg states only
if ~isempty(v(1).psd)
    prct = v(1).psd.info.input.prct;
    sstates = v(1).psd.info.input.sstates;
else
    prct = 70;
    sstates = [1, 4];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze and organize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize cell array of state ratio (nfiles x unitType x ntiles)
dataMat = cell(nfiles, 2, ntiles);
for ifile = 1 : nfiles
    
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    
    if flgAnalyze
        
        % get window from psd
        if ~isempty(v(ifile).psd)
            wins = v(ifile).psd.info.wins;
        else
            wins = [0, Inf];
        end
        
        if flgEmg

            % load emg data
            sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
            emg = load(sleepfile, 'emg_rms');
            emg = emg.emg_rms;

            % get indices to high- and low-emg
            labels = double(emg > prctile(emg, prct));
            labels(emg < prctile(emg, 100 - prct)) = 4;

            % limit indices to time window and get "state" epochs
            labels = labels(wins(:, 1) : wins(:, 2));
            [stateEpochs, ~] = as_epochs('labels', labels, 'minDur', 10, 'interDur', 4);
        else
            stateEpochs = v(ifile).ss.stateEpochs;
        end

        % calc firing rate
        v(ifile).fr = calc_fr(v(ifile).spikes.times, 'basepath', basepath,...
            'graphics', false, 'binsize', 60, 'saveVar', true,...
            'smet', 'none', 'winBL', wins, 'winCalc', [0, Inf],...
            'stateEpochs', stateEpochs);
        stateRat = squeeze(v(ifile).fr.states.ratio(4, 1, :));

    else
        stateRat = squeeze(v(ifile).fr.states.ratio(4, 1, :));
    end

    % organize cell array of state ratio
    for iunit = 1 : 2
        units = v(ifile).units.clean(iunit, :);
        
        % get a cell array of indices according to mfr percentiles per unit
        % type
        [~, tileIdx] = vec2tileMat(v(ifile).fr.mfr, ntiles, units);

        % organize state ratio according to percentiles
        frTiles(ifile, iunit, :) = cellfun(@(x) stateRat(x), tileIdx, 'uni', false);
    end
end

% save
if saveVar
    mousepath = fileparts(basepaths{1});
    datafile = fullfile(mousepath, [mname, '_frTiles.mat']);
    save(datafile, 'frTiles')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % graphics and to prism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics

    iunit = 1;
    unitClr = {'b', 'r'};
    unitType = {'RS', 'FS'};
    if flgEmg
        stateNames = {'High-EMG', 'Low-EMG'};
    else
        stateNames = v(1).ss.info.names(sstates);
    end

    fh = figure;
    fh.Position = [0.1 0.1 0.8 0.8];
    th = tiledlayout(1, ntiles, 'TileSpacing', 'Compact');
    for itile = 1 : ntiles
        axh = nexttile;
        ydata = cell2nanmat(frTiles(:, iunit, itile), 2);
        plot_boxMean('dataMat', ydata, 'clr', unitClr{iunit}, 'allPnts', true)
        hold on
        plot(xlim, [0, 0], '--k')
        ylabel({sprintf('%s - %s /', stateNames{2}, stateNames{1}),...
            sprintf('%s + %s', stateNames{2}, stateNames{1})})
        xlabel('Session')
        title(sprintf('Tile #%d', itile))

    end
    title(th, mname)

end

end

% EOF