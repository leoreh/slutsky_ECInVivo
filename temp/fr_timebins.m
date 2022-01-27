function frBins = fr_timebins(varargin)

% INPUT:
%   basepath        string. path to recording folder {pwd}
%   timebins        n x 2 numeric. timebins to calc fr [s]
%   sstates         numeric. selected states (indexed according to as_cfg)
%   saveVar         logical {true}
%   tbins_txt       cell of chars for xticklabels
%   forceA          logical. reanalyze recording even if fr_bins var exists
%                   {false}
%   graphics        logical {true}
%   saveFig         logical {true}
%
% DEPENDENCIES
%   firingRate
%
% TO DO LIST
%
% 07 jan 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'timebins', [], @isnumeric);
addOptional(p, 'sstates', [], @isnumeric);
addOptional(p, 'tbins_txt', []);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
timebins        = p.Results.timebins;
sstates         = p.Results.sstates;
tbins_txt       = p.Results.tbins_txt;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;
saveFig         = p.Results.saveFig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(tbins_txt) & size(timebins, 1) == 8
    tbins_txt = {'0-3ZT', '3-6ZT', '6-9ZT', '9-12ZT',...
        '12-15ZT', '15-18ZT', '18-21ZT', '21-24ZT'};
else
    tbins_txt = split(num2str(1 : size(timebins, 1)));
end

unitchar = {'RS', 'FS'};

% load vars from each session
varsFile = ["fr"; "sr"; "spikes"; "cell_metrics"; "datInfo"; "session"];
varsName = ["fr"; "sr"; "spikes"; "cm"; "datInfo"; "session"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% file
cd(basepath)
[~, basename] = fileparts(basepath);
frfile = fullfile(basepath, [basename, '.fr_bins.mat']);

% timebins
if isempty(timebins)
    if isfield(v.session.general, 'timebins')
        timebins = v.session.general.timebins;
    else
        error('no timebins specified')
    end
end
nwin = size(timebins, 1);

% states
cfg = as_loadConfig();
nstates = cfg.nstates;
if isempty(sstates)
    sstates = 1 : nstates;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc firing rate in time bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if already analyzed
if exist(frfile, 'file') && ~forceA
    load(frfile)
else
    for iwin = 1 : nwin        
        frBins(iwin) = firingRate(v.spikes.times,...
            'basepath', basepath, 'graphics', false,...
            'binsize', 60, 'saveVar', false, 'smet', 'GK', 'winBL',...
            [0, Inf], 'winCalc', timebins(iwin, :), 'forceA', true);
        
    end
end

% re-organize vars of interest in 3d mat of state x win x unit
for iwin = 1 : nwin
    for istate = 1 : length(sstates)
        stateMfr(istate, iwin, :) = mean(frBins(iwin).states.fr{istate}, 2, 'omitnan');
        stateRat(istate, iwin, :) = squeeze(frBins(iwin).states.ratio(1, istate, :));
    end
    stateGain(:, iwin, :) = frBins(iwin).states.gain;
end
units = selectUnits('basepath', basepath);
units = units.idx;

if saveVar
    save(frfile, 'frBins')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
       
    for iunit = 1 : 2
        fh = figure;
        yLimitMfr = [min(stateMfr(:, :, units(iunit, :)), [], 'all'),...
            max(stateMfr(:, :, units(iunit, :)), [], 'all')];
        yLimitRatio = [min(stateRat(:, :, units(iunit, :)), [], 'all'),...
            max(stateRat(:, :, units(iunit, :)), [], 'all')];
    
        for istate = 1 : length(sstates)
            
            % mfr in state
            subplot(2, length(sstates), istate)
            dataMat = squeeze(stateMfr(istate, :, units(iunit, :)));
            plot_boxMean('dataMat', dataMat', 'clr', cfg.colors{sstates(istate)})
            ylabel(sprintf('MFR %s', cfg.names{sstates(istate)}))
            ylim(yLimitMfr)
            xticklabels(tbins_txt)
            xtickangle(45)
            
            % state ratio
            if sstates(istate) ~= 1
                subplot(2, length(sstates), istate + length(sstates))
                dataMat = squeeze(stateRat(istate, :, units(iunit, :)));
                plot_boxMean('dataMat', dataMat', 'clr', cfg.colors{sstates(istate)})
                ylabel({sprintf('%s - %s /', cfg.names{1}, cfg.names{sstates(istate)}),...
                    sprintf('%s + %s', cfg.names{1}, cfg.names{sstates(istate)})})
                ylim(yLimitRatio)
                xticklabels(tbins_txt)
                xtickangle(45)
            end                   
        end       
        sgtitle([basename, '_', unitchar{iunit}])
        
        if saveFig
            figpath = fullfile(basepath, 'graphics', 'sleepState');
            mkdir(figpath)
            figname = fullfile(figpath, [basename, '_frBins_', unitchar{iunit}]);
            export_fig(figname, '-tif', '-transparent', '-r300')
        end
    end
end

end

% EOF


