function psdBins = psd_states_timebins(varargin)

% calculating the psd per state in different time bins of a
% session.
%
% INPUT:
%   basepath        string. path to recording folder {pwd}
%   fs              numeric. sampling frequency of lfp data. if empty will
%                   be extracted from session struct. also determines the
%                   new sampling frequency of the eeg signal.
%   fsEeg           numeric. sampling frequency of eeg data
%   chEeg           numeric. channels of eeg signal. if greater than 2
%                   will assume it is recorded in same lfp file. otherwise
%                   will assume it is recorded in emg.dat file. if empty
%                   will assume no eeg recording exists. 
%   nchansEeg       numeric. number of channels in [basename].emg.dat file
%   nbins           numeric. how many time bins to separate the recording
%   timebins        numeric. points of interest. chunks will be created
%                   in relation to these poitns [sec].  
%   saveVar         logical. save ss var {true}
%   forceA          logical. reanalyze recordings even if ss struct
%                   exists (false)
%   graphics        logical. plot confusion chart and state separation {true}
%
% DEPENDENCIES:
%
% TO DO LIST:
%
% 07 jan 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fsEeg', [], @isnumeric);
addOptional(p, 'chEeg', 1, @isnumeric);
addOptional(p, 'nchansEeg', [], @isnumeric);
addOptional(p, 'timebins', [], @isnumeric);
addOptional(p, 'tbins_txt', []);
addOptional(p, 'sstates', [1, 4, 5], @isnumeric);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
fsEeg           = p.Results.fsEeg;
chEeg           = p.Results.chEeg;
nchansEeg       = p.Results.nchansEeg;
timebins        = p.Results.timebins;
tbins_txt       = p.Results.tbins_txt;
sstates         = p.Results.sstates;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tbins_txt = {'0-3ZT', '3-6ZT', '6-9ZT', '9-12ZT',...
    %             '12-15ZT', '15-18ZT', '18-21ZT', '21-24ZT'};

% file
cd(basepath)
[~, basename] = fileparts(basepath);
psdfile = fullfile(basepath, [basename, '.psd_bins.mat']);

% load session vars
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% recording params
nchans = v.session.extracellular.nChannels;
fsDat = v.session.extracellular.sr;
fsLfp = v.session.extracellular.srLfp;
if isempty(fsEeg)
    fsEeg = v.ss.info.sSig.emgFs;
end
if isempty(nchansEeg)
    nchansEeg = nchans;
end
chLfp = v.ss.info.sSig.eegCh;

% timebins
% timebins(isinf(timebins)) = floor(v.session.general.duration);
nbins = size(timebins, 1);

if isempty(tbins_txt) && nbins == 4
    tbins_txt = {'0-6 ZT', '6-12 ZT', '12-18 ZT', '18-24 ZT'};
elseif isempty(tbins_txt)
    tbins_txt = split(num2str(1 : nbins));
end

% params for filtering eeg
import iosr.dsp.*
filtRatio = 450 / (fsEeg / 2);
fsRatio = (fsEeg / fsLfp);

% state params
faxis = 0.5 : 0.2 : 100;

% initialize
psdLfp = nan(nbins, length(sstates), length(faxis));
psdEeg = nan(nbins, length(sstates), length(faxis));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if already analyzed
if exist(psdfile, 'file') && ~forceA
    load(psdfile)
else
    for iwin = 1 : nbins
                        
        % get state epochs
        labels = v.ss.labels(timebins(iwin, 1) : timebins(iwin, 2));
        [stateEpochs, epStats(iwin)] = as_epochs('labels', labels);

        % ---------------------------------------------------------------------
        % hippocampal lfp
        sig = double(bz_LoadBinary([basename, '.lfp'],...
            'duration', diff(timebins(iwin, :)) + 1,...
            'frequency', fsLfp, 'nchannels', nchans, 'start', timebins(iwin, 1),...
            'channels', chLfp, 'downsample', 1));
        sig = mean(sig, 2);

        psdLfp(iwin, :, :) = psd_timebins('sig', sig, 'winCalc', stateEpochs(sstates),...
            'fs', fsLfp, 'faxis', faxis, 'graphics', false);  
        
        % ---------------------------------------------------------------------
        % frontal eeg
        if ~isempty(chEeg)
            if chEeg < 3
                sig = double(bz_LoadBinary([basename, '.emg.dat'],...
                    'duration', diff(timebins(iwin, :)) + 1,...
                    'frequency', fsEeg, 'nchannels', 2, 'start', timebins(iwin, 1),...
                    'channels', chEeg, 'downsample', 1));
                if fsEeg ~= fsLfp
                    sig = [iosr.dsp.sincFilter(sig, filtRatio)]';
                    sig = sig(fsRatio : fsRatio : length(sig));
                end
            else
                sig = double(bz_LoadBinary([basename, '.lfp'],...
                    'duration', diff(timebins(iwin, :)) + 1,...
                    'frequency', fsLfp, 'nchannels', nchansEeg, 'start', timebins(iwin, 1),...
                    'channels', chEeg, 'downsample', 1));
            end
            psdEeg(iwin, :, :) = psd_timebins('sig', sig, 'winCalc', stateEpochs(sstates),...
                'fs', fsLfp, 'faxis', faxis, 'graphics', false);

        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % arrange in struct and save
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % organize epoch lengths according to time bins
    for iwin = 1 : nbins
        for istate = 1 : length(sstates)
            epLen_temp{iwin, istate} = epStats(iwin).epLen{sstates(istate)};
        end
        totDur(iwin, :) = epStats(iwin).totDur(sstates);
    end
    for istate = 1 : length(sstates)
        epLen{istate} = cell2nanmat(epLen_temp(:, istate), 2);
    end
    
    psdBins.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
    psdBins.info.chLfp = chLfp;
    psdBins.info.chEeg = chEeg;
    psdBins.info.fsLfp = fsLfp;
    psdBins.info.fsEeg = fsEeg;
    psdBins.info.freq = faxis;
    psdBins.timebins = timebins;
    psdBins.psdLfp = psdLfp;
    psdBins.psdEeg = psdEeg;
    psdBins.epStats = epStats;
    psdBins.epLen_organized = epLen;
    psdBins.totDur_organized = totDur;
    
    if saveVar
        save(psdfile, 'psdBins')
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    if all(isnan(psdBins.psdEeg), 'all')
        nrows = 2;
    else
        nrows = 3;
    end

    setMatlabGraphics(false)
    alphaIdx = linspace(0.5, 1, nbins);
    cfg = as_loadConfig();
    lim_fAxis = faxis >= 1;
    smf = 17;
    gk = gausswin(smf);
    gk = gk / sum(gk);
    
    fh = figure;
    fh.Position = [0.1 0.1 0.8 0.8];

    for istate = 1 : length(sstates)
        
        % lfp
        subplot(nrows, length(sstates), istate)
        psdMat = squeeze(psdBins.psdLfp(:, istate, lim_fAxis));
%         for iwin = 1 : nbins
%             psdMat(iwin, :) = conv(psdMat(iwin, :), gk, 'same');
%         end
%         psdMat = psdMat ./ sum(psdMat(:, 21 : end), 2);
        ph = plot(faxis(lim_fAxis), psdMat', 'LineWidth', 2);
        for iwin = 1 : length(ph)
            ph(iwin).Color(istate) = cfg.colors{sstates(istate)}(istate) - iwin * 0.003;
            ph(iwin).Color(4) = alphaIdx(iwin);
        end
        title(cfg.names{sstates(istate)})
        xlabel('Frequency [Hz]')
        ylabel('LFP PSD [mV^2/Hz]')
        ylabel('Norm. LFP PSD')
        legend(tbins_txt, 'Location', 'Southwest', 'NumColumns', 2, 'FontSize', 9)
        set(gca, 'YScale', 'log', 'XScale', 'log')
        xlim([faxis(find(lim_fAxis, 1)), faxis(find(lim_fAxis, 1, 'last'))])
        
        % eeg
        if nrows == 3
            subplot(nrows, length(sstates), istate + length(sstates))
            psdMat = squeeze(psdBins.psdEeg(:, istate, lim_fAxis));
            %     for iwin = 1 : nbins
            %         psdMat(iwin, :) = conv(psdMat(iwin, :), gk, 'same');
            %     end
            psdMat = psdMat ./ sum(psdMat, 2);
            ph = plot(faxis(lim_fAxis), psdMat', 'LineWidth', 2);
            for iwin = 1 : length(ph)
                ph(iwin).Color(istate) = cfg.colors{sstates(istate)}(istate) - iwin * 0.003;
                ph(iwin).Color(4) = alphaIdx(iwin);
            end
            xlabel('Frequency [Hz]')
            ylabel('EEG PSD [mV^2/Hz]')
            ylabel('Norm. EEG PSD')
            set(gca, 'YScale', 'log', 'XScale', 'log')
            xlim([faxis(find(lim_fAxis, 1)), faxis(find(lim_fAxis, 1, 'last'))])
        end
        
        % state duration
        subplot(nrows, length(sstates), istate + (nrows - 1) * length(sstates))
        epMat = psdBins.epLen_organized{istate};
        boxplot(epMat, 'PlotStyle', 'traditional', 'Whisker', 100);
        bh = findobj(gca, 'Tag', 'Box');
        bh = flipud(bh);
        for ibox = 1 : length(bh)
            patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
                cfg.colors{sstates(istate)}, 'FaceAlpha', alphaIdx(ibox))
        end
        ylabel('Epoch Length [log(s)]')
        set(gca, 'YScale', 'log')
        ylim([0 ceil(prctile(epMat(:), 99.99))])
        yyaxis right
        plot([1 : size(epMat, 2)], psdBins.totDur_organized(:, istate) / 60,...
            'kd', 'markerfacecolor', 'k')
        ylabel('State duration [min]')
        ax = gca;
        set(ax.YAxis(1), 'color', cfg.colors{sstates(istate)})
        set(ax.YAxis(2), 'color', 'k')
        xticklabels(tbins_txt)
        xtickangle(45)
        
    end
    sgtitle(basename)
    
    saveFig = true;
    if saveFig
        figpath = fullfile(basepath, 'graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, [basename, '_psdBins']);
        export_fig(figname, '-jpg', '-transparent', '-r300')
    end
    
    
end

