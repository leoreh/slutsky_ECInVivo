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
%   saveVar         logical. save ss var {true}
%   forceA          logical. reanalyze recordings even if ss struct
%                   exists (false)
%   graphics        logical. plot confusion chart and state separation {true}
%
% DEPENDENCIES:
%   AccuSleep (modified in slutskycode)
%   IOSR.DSP.SINCFILTER     for filtering data
%   binary2epochs
%
% TO DO LIST:
%       # filter before resampling to assure nyquist (done)
%       # implement forceLoad (done)
%       # uigetf for net (done)
%       # graphics (done)
%       # implement tsa filter
%       # batch processing
%       # load sig instead of input (memory)
%
% 07 jan 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fsEeg', 6103.515625, @isnumeric);
addOptional(p, 'chEeg', 1, @isnumeric);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
fsEeg           = p.Results.fsEeg;
chEeg           = p.Results.chEeg;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file
cd(basepath)
[~, basename] = fileparts(basepath);
psdfile = fullfile(basepath, [basename, '.psdBins.mat']);

% load session vars
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% recording params
fileinfo = dir([basename, '.dat']);
recLen = floor(fileinfo.bytes / 2 / nchans / fsDat);
fsDat = v.session.extracellular.sr;
fsLfp = v.session.extracellular.srLfp;
nchans = v.session.extracellular.nChannels;
chLfp = v.ss.info.sSig.eegCh;

% params for filtering eeg
import iosr.dsp.*
filtRatio = 450 / (fsEeg / 2);
fsRatio = (fsEeg / fsLfp);

% state params
faxis = 0.2 : 0.2 : 120;
sstates = [1, 4, 5];

% separate recording to timebins
nbins = 4;
csec = floor(cumsum(v.datInfo.nsamps) / fsDat);
[~, injIdx] = min(abs(csec - 6 * 60 * 60));
injTime = csec(injIdx);
timebins = n2chunks('n', recLen, 'chunksize', ceil(recLen / nbins));
timebins(1, 2) = floor(injTime);
timebins(2, 1) = ceil(injTime) + 1;
timebins(end, end) = length(v.ss.labels);

tbins_txt = {'0-6 ZT', '6-12 ZT', '12-18 ZT', '18-24 ZT'};

% initialize
psdLfp = nan(nbins, length(sstates), length(faxis));
psdEeg = nan(nbins, length(sstates), length(faxis));

% check if already analyzed
if exists(psdfile) && ~forceA
    load(psdfile, 'psdBins')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iwin = 1 : nbins
    
    labels = v.ss.labels(timebins(iwin, 1) : timebins(iwin, 2));
    
    % ---------------------------------------------------------------------
    % hippocampal lfp
    sig = double(bz_LoadBinary([basename, '.lfp'],...
        'duration', diff(timebins(iwin, :)) + 1,...
        'frequency', fsLfp, 'nchannels', nchans, 'start', timebins(iwin, 1),...
        'channels', chLfp, 'downsample', 1));
    sig = mean(sig, 2);

    [psdLfp(iwin, :, :), ~, epStats(iwin)] = psd_states('eeg', sig,...
        'labels', labels, 'fs', fsLfp, 'faxis', faxis,...
        'graphics', false, 'sstates', sstates);

    % ---------------------------------------------------------------------
    % frontal eeg
    if ~isempty(chEeg)
        if chEeg < 3
            sig = double(bz_LoadBinary([basename, '.emg.dat'],...
                'duration', diff(timebins(iwin, :)) + 1,...
                'frequency', fsEeg, 'nchannels', 2, 'start', timebins(iwin, 1),...
                'channels', chEeg, 'downsample', 1));
            sig = [iosr.dsp.sincFilter(sig, filtRatio)]';
            sig = sig(fsRatio : fsRatio : length(sig));
        else
            sig = double(bz_LoadBinary([basename, '.lfp'],...
                'duration', diff(timebins(iwin, :)) + 1,...
                'frequency', fsLfp, 'nchannels', nchans, 'start', timebins(iwin, 1),...
                'channels', chEeg, 'downsample', 1));
        end
        [psdEeg(iwin, :, :), ~, epStats(iwin)] = psd_states('eeg', sig,...
            'labels', labels, 'fs', fsLfp, 'faxis', faxis,...
            'graphics', false, 'sstates', sstates);
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
    epLen{istate} = cell2nanmat(epLen_temp(:, istate));
end

psdBins.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
psdBins.info.chLfp = chLfp;
psdBins.info.chEeg = chEeg;
psdBins.info.fsLfp = fsLfp;
psdBins.info.fsEeg = fsEeg;
psdBins.timebins = timebins;
psdBins.psdLfp = psdLfp;
psdBins.psdEeg = psdEeg;
psdBins.epStats = epStats;
psdBins.epLen_organized = epLen;
psdBins.totDur_organized = totDur;

if saveVar
    save(psdfile, 'psdBins')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    setMatlabGraphics(false)
    alphaIdx = linspace(0.5, 1, nbins);
    cfg = as_loadConfig();
    lim_fAxis = faxis > 1;
    smf = 7;
    gk = gausswin(smf);
    gk = gk / sum(gk);
    
    fh = figure;
    for istate = 1 : length(sstates)
        
        % lfp
        subplot(3, 3, istate)
        psdMat = squeeze(psdBins.psdLfp(:, istate, lim_fAxis));
        %     for iwin = 1 : nbins
        %         psdMat(iwin, :) = conv(psdMat(iwin, :), gk, 'same');
        %     end
        %     psdMat = psdMat ./ sum(psdMat, 2);
        ph = plot(faxis(lim_fAxis), psdMat', 'LineWidth', 2);
        for iwin = 1 : length(ph)
            ph(iwin).Color(istate) = cfg.colors{sstates(istate)}(istate) - iwin * 0.003;
            ph(iwin).Color(4) = alphaIdx(iwin);
        end
        title(cfg.names{sstates(istate)})
        xlabel('Frequency [Hz]')
        ylabel('LFP PSD [mV^2/Hz]')
        legend(tbins_txt, 'Location', 'Southwest')
        set(gca, 'YScale', 'log', 'XScale', 'log')
        xlim([faxis(find(lim_fAxis, 1)), faxis(find(lim_fAxis, 1, 'last'))])
        
        % eeg
        subplot(3, 3, istate + 3)
        psdMat = squeeze(psdBins.psdEeg(:, istate, lim_fAxis));
        %     for iwin = 1 : nbins
        %         psdMat(iwin, :) = conv(psdMat(iwin, :), gk, 'same');
        %     end
        %     psdMat = psdMat ./ sum(psdMat, 2);
        ph = plot(faxis(lim_fAxis), psdMat', 'LineWidth', 2);
        for iwin = 1 : length(ph)
            ph(iwin).Color(istate) = cfg.colors{sstates(istate)}(istate) - iwin * 0.003;
            ph(iwin).Color(4) = alphaIdx(iwin);
        end
        xlabel('Frequency [Hz]')
        ylabel('EEG PSD [mV^2/Hz]')
        set(gca, 'YScale', 'log', 'XScale', 'log')
        xlim([faxis(find(lim_fAxis, 1)), faxis(find(lim_fAxis, 1, 'last'))])
        
        % state duration
        subplot(3, 3, istate + 6)
        epMat = epLen{istate};
        boxplot(epMat, 'PlotStyle', 'traditional', 'Whisker', 6);
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
        plot([1 : size(epMat, 2)], totDur(:, istate) / 60,...
            'kd', 'markerfacecolor', 'k')
        ylabel('State duration [min]')
        ax = gca;
        set(ax.YAxis(1), 'color', cfg.colors{sstates(istate)})
        set(ax.YAxis(2), 'color', 'k')
        xticklabels(tbins_txt)
        xtickangle(45)
        
    end
    
    saveFig = true;
    if saveFig
        figpath = fullfile(basepath, 'graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, [basename, '_psdBins']);
        export_fig(figname, '-jpg', '-transparent', '-r300')
    end
    
    
end

