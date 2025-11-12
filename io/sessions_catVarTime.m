function [expData, info] = sessions_catVarTime(varargin)

% concatenates a variable from different sessions. assumes sessions are
% contineous. concatenates according to the time of day extracted from
% basenames. 
% 
% INPUT
%   basepaths       cell of chars to recording sessions
%   mname           char of mouse name. if basepaths is empty will get
%                   basepaths from the sessionList.xlsx file
%   xTicksBinsize   numeric. x-axis ticks every binsize hr
%   graphics        logical {true}
%   saveFig         logical {true}
%   markRecTrans    logical {false}. add vertical lines on the transition
%                   between recordings
%   dataPreset      string or cell of string depicting the variable to cat. 
%                   can be any combination of 'sr', 'spec', 'fr', 'ripp',
%                   'emg_rms', 'bands', 'srsu', 'hypnogram',
%                   'hypnogram_emg', or 'spec_eeg'
%   dataAlt         numeric. data alternative. relavent for some presets:
%                   for 'spec'  represents channel
%                   for 'bands' represents channel
%                   for 'srsu'  represents unittype (1 = rs; 2 = fs)
%   hAx             handle to plot axis      
% 
% EXAMPLE
% mname = 'lh96';
% [srData, tidx, tidxLabels] = sessions_catVarTime('mname', mname, 'dataPreset', 'both', 'graphics', false);
% 
% TO DO LIST
% 
% UPDATES
%   12 feb 22       added ripples
%   18 apr 24       organized time info

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'basepaths', {}, @iscell);
addParameter(p, 'mname', '', @ischar);
addParameter(p, 'dataPreset', 'sr');
addParameter(p, 'dataAlt', 1, @isnumeric);
addParameter(p, 'xTicksBinsize', 12, @isnumeric);
addParameter(p, 'graphics', true, @islogical);
addParameter(p, 'saveFig', true, @islogical);
addParameter(p, 'markRecTrans', true, @islogical);
addParameter(p, 'hAx', []);

parse(p, varargin{:})
basepaths       = p.Results.basepaths;
mname           = p.Results.mname;
dataPreset      = p.Results.dataPreset;
dataAlt         = p.Results.dataAlt;
xTicksBinsize   = p.Results.xTicksBinsize;
graphics        = p.Results.graphics;
saveFig         = p.Results.saveFig;
markRecTrans    = p.Results.markRecTrans;
hAx             = p.Results.hAx;

zt0 = guessDateTime('0900');    % lights on at 09:00 AM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run recursively through the function if several data vars are requested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iscell(dataPreset)
    hFig = figure;
    set(hFig, 'WindowState', 'maximized');
    th = tiledlayout(length(dataPreset), 1);
    th.TileSpacing = 'tight';
    th.Padding = 'none';
    set(hFig, 'DefaultAxesFontSize', 16);
    title(th, mname, 'interpreter', 'none', 'FontSize', 20)

    for idata = 1 : length(dataPreset)
        sb(idata) = nexttile;
        [expData{idata}, info{idata}] = sessions_catVarTime('mname', mname,...
            'basepaths', basepaths, 'dataPreset', dataPreset{idata},...
            'graphics', graphics, 'hAx', sb(idata), 'dataAlt', dataAlt,...
            'xTicksBinsize', xTicksBinsize, 'markRecTrans', markRecTrans);
        
        if idata < length(dataPreset)
            set(sb(idata), 'xticklabels', {[]})
            xlabel('')
        end
    end

    linkaxes(sb, 'x')
    axis tight

    % concatenate info struct
    info = catfields([info{:}], 2);

    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varsFile = ["units"; "datInfo"; "session"; string(dataPreset)];
varsName = ["units"; "datInfo"; "session"; string(dataPreset)];

if isempty(basepaths)
    [v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"]);
else
    [v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"]);
end

[~, basenames] = cellfun(@fileparts, basepaths, 'uni', false);
nsessions = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data to concatenate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
fs = v(1).session.extracellular.sr;

switch dataPreset
    case 'spec'       
        for isession = 1 : nsessions
            if ndims(v(isession).spec.s) == 3
                v(isession).data = squeeze(v(isession).spec.s(:, :, dataAlt))';
            else
                v(isession).data = v(isession).spec.s';
            end
        end
        ts = v(isession).spec.info.winstep;     % sampling period [s]            
        faxis = v(isession).spec.freq;
    
    case 'spec_eeg'       
        for isession = 1 : nsessions
            if ndims(v(isession).spec_eeg.s) == 3
                v(isession).data = squeeze(v(isession).spec_eeg.s(:, :, 1))';
            else
                v(isession).data = v(isession).spec_eeg.s';
            end
        end
        ts = v(isession).spec_eeg.info.winstep;     % sampling period [s]            
        faxis = v(isession).spec_eeg.freq;

    case 'bands'
        for isession = 1 : nsessions
            filename = fullfile(basepaths{isession},...
                [basenames{isession}, '.spec.mat']);
            load(filename, 'spec');
            yval = movmean(squeeze(spec.bands.db(:, :, dataAlt)), 100,  2);
            v(isession).data = yval;
        end
        ts = spec.info.winstep;                          % sampling period [s]
    
    case 'hypnogram'
        for isession = 1 : nsessions
            filename = fullfile(basepaths{isession},...
                [basenames{isession}, '.sleep_states.mat']);
            load(filename, 'ss');
            v(isession).data = ss.labels';
        end
        ts = ss.info.boutLen;                           % sampling period [s]
        cfg = as_loadConfig('flgEmg', false);
        sstates = [];

    case 'hypnogram_emg'
        for isession = 1 : nsessions
            filename = fullfile(basepaths{isession},...
                [basenames{isession}, '.sleep_statesEmg.mat']);
            load(filename, 'ssEmg');
            v(isession).data = ssEmg.labels';
        end
        ts = ssEmg.info.boutLen;                           % sampling period [s]
        cfg = as_loadConfig('flgEmg', true);
        sstates = [1, 2];

    case 'emg_rms'
        for isession = 1 : nsessions
            filename = fullfile(basepaths{isession},...
                [basenames{isession}, '.sleep_sig.mat']);
            load(filename, 'emg_rms');
            v(isession).data = emg_rms;
        end
        cfg = as_loadConfig();
        ts = cfg.boutLen;                              % sampling period [s]        
    
    case 'sr'
        for isession = 1 : nsessions
            v(isession).data = v(isession).sr.strd;
        end
        ts = v(1).sr.info.binsize;                      % sampling period [s]
        grp = 1 : v(1).session.extracellular.nSpikeGroups;
        grp = dataAlt;
  
    case 'srsu'
        for isession = 1 : nsessions
            v(isession).data = v(isession).srsu(dataAlt).strd;
        end
        ts = v(1).srsu(dataAlt).info.binsize;          % sampling period [s]
        grp = 1 : v(1).session.extracellular.nSpikeGroups;

    case 'fr'
        % nan pad each session to the max number of units
        nunits = [];
        for isession = 1 : nsessions
            datasz(isession, :) = size(v(isession).fr.strd);
            units = v(isession).units.clean;
            nunits(isession, :) = [sum(units(1, :)), sum(units(2, :))];
            v(isession).rs = v(isession).fr.strd(units(1, :), :);
            v(isession).fs = v(isession).fr.strd(units(2, :), :);
        end
        for isession = 1 : nsessions
            v(isession).data = [v(isession).fr.strd;...
                nan(max(datasz(:, 1)) - datasz(isession, 1),...
                datasz(isession, 2))];
            v(isession).rs = [v(isession).rs;...
                nan(max(nunits(:, 1)) - nunits(isession, 1),...
                datasz(isession, 2))];
            v(isession).fs = [v(isession).fs;...
                nan(max(nunits(:, 2)) - nunits(isession, 2),...
                datasz(isession, 2))];
        end
        ts = v(1).fr.info.binsize;          % sampling period [s]        
        
    case 'ripp'
         for isession = 1 : nsessions
            v(isession).data = v(isession).ripp.rate.rate';
        end
        ts = mode(diff(v(isession).ripp.rate.binedges{1}));                 
        
    otherwise
        sprintf('\nno such data preset\n')
end

ncol = size(v(1).data, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize data mat for all sessions based on the time from the first to
% the last session, and the sampling frequency of the variable. assumes the
% recordings are contineous. 
t_recStart = cellfun(@guessDateTime, basenames, 'uni', true);
t_expStart = t_recStart(1) - max([0, diff(timeofday([zt0, t_recStart(1)]))]);
t_expEnd = guessDateTime(basenames{end});
lastDur = v(end).session.general.duration;
if isstr(lastDur)
    lastDur = str2num(lastDur);
end
t_expEnd = t_expEnd + seconds(lastDur);
expLen = ceil(seconds(t_expEnd - t_expStart) / ts); 
expData = nan(expLen, ncol);
expRs = nan(expLen, ncol);
expFs = nan(expLen, ncol);

% initialize
x_blockTrans = [];
for isession = 1 : nsessions
    
    % find index to dataMat according to recording start
    recIdx = round(max([1, seconds(t_recStart(isession) - t_expStart) / ts]));
    recData = v(isession).data;
    recLen = length(recData);
    
    % insert recording data to experiment data
    expData(recIdx : recIdx + recLen - 1, :) = recData';
    
    if strcmp(dataPreset, 'fr')
        expRs(recIdx : recIdx + recLen - 1, 1 : max(nunits(:, 1))) = v(isession).rs';
        expFs(recIdx : recIdx + recLen - 1, 1 : max(nunits(:, 2))) = v(isession).fs';
    end
    
    % cat block transitions. The new line stores only the transitions
    % between files (eg, different spike sorting) and not between blocks
    % if ~isempty(v(isession).datInfo) && isfield(v(isession).datInfo, 'nsamps')
    %     x_blockTrans = [x_blockTrans, cumsum(v(isession).datInfo.nsamps) / fs / ts + recIdx];
    % else
    %     x_blockTrans = 0;
    % end
    
    x_blockTrans = [x_blockTrans, recIdx];

end

% create timstamps and labels for x-axis
xTicksSamples = xTicksBinsize * 60 * 60;   % [seconds]
zt0idx = round(seconds(diff(timeofday([t_recStart(1), zt0])))) / ts;
tidx = round([zt0idx : xTicksSamples / ts : expLen]);
tStartLabel = datetime(t_expStart.Year, t_expStart.Month, t_expStart.Day,...
    zt0.Hour, zt0.Minute, zt0.Second);
x_labels = string(datestr(datenum(tStartLabel : hours(xTicksBinsize) : t_expEnd),...
    'HH:MM', 2000));
x_ticks = tidx * ts;

% add date to x labels once a day
x_labels(1 : 24 / xTicksBinsize : end) =...
    string(datestr(datenum(tStartLabel : hours(24) : t_expEnd),...
    'yymmdd', 2000));

% data for x-axis
x_data = [1 : ts : ceil(seconds(t_expEnd - t_expStart))];

% assign fr to output
if strcmp(dataPreset, 'fr')
    if dataAlt == 1
        expData = expRs;
    else
        expData = expFs;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ogranize time info struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.x_data = x_data';
info.x_ticks = x_ticks';
info.x_labels = x_labels;
info.x_blockTrans = x_blockTrans';
info.t_blockTrans = tstamp2time('tstamp', x_blockTrans * ts, 'dtstr', t_expStart)';
info.t_recStart = t_recStart';
info.t_expStart = t_expStart;
info.t_expEnd = t_expEnd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    if isempty(hAx)
        hFig = figure;
        set(hFig, 'WindowState','maximized');
        hAx = gca;
    else
        set(gcf, 'CurrentAxes', hAx)
        set(gcf, 'WindowState','maximized');
    end
    switch dataPreset
        case {'spec', 'spec_eeg'}
            
            % take a sample of the spectrogram to help initialize the colormap
            sampleBins = randperm(expLen, round(expLen / 10));
            specSample = reshape(expData(sampleBins, :), 1, length(sampleBins) * length(faxis));
            caxis1 = prctile(specSample, [2 94]);

            % plot
            imagesc(x_data, faxis', expData', caxis1);
            colormap(AccuSleep_colormap());
            axis('xy')
            ylabel('Freq. [Hz]')
            set(gca, 'yscale', 'log')
            ylim([max([0.2, faxis(1)]), faxis(end)])
            yticks([1, 4, 16, 64])

        case 'bands'
            expData = movmean(expData, 33, 1);
            plot(x_data, expData, 'lineWidth', 2.5);
            ylabel('Spectral Power [dB]')
            
            for iband = 1 : length(spec.bands.bandFreqs)
                lgd{iband} = sprintf('%s [%d-%d Hz]',...
                    spec.bands.bandNames{iband}, floor(spec.bands.bandFreqs(iband, 1)),...
                    spec.bands.bandFreqs(iband, 2));
            end
            legend(lgd{1 : end}, 'location', 'northeast')
        
        case {'hypnogram', 'hypnogram_emg'}
            plot_hypnogram('labels', expData, 'hAx', hAx,...
                'clr', cfg.colors, 'sstates', sstates)
            pbaspect([30, 1, 1])

        case 'emg_rms'            
            % plot
            expData = movmean(expData, 33, 1);
            plot(x_data, expData);
            ylabel('EMG RMS')
            
        case {'sr', 'srsu'}
            expData = movmean(expData, 13, 1);
            plot(x_data, expData(:, grp))
            ylabel('MU FR [Hz]')
            legend(split(num2str(grp)))
            axis tight
            
        case 'fr'
            hold on
            plotData = movmean(expRs, 13, 1);
            plot_stdShade('dataMat', plotData, 'xVal', x_data,...
                'hAx', hAx, 'clr', [0, 0, 1]);

            plotData = movmean(expFs, 13, 1);
            plot_stdShade('dataMat', plotData, 'xVal', x_data,...
                'hAx', hAx, 'clr', [1, 0, 0]);

            legend({sprintf('RS <= %d', max(nunits(:, 1))),...
                sprintf('FS <= %d', max(nunits(:, 2)))})
            ylabel('MFR [Hz]')
            axis tight


        case 'ripp'
            plot(x_data, expData)
            ylabel('Ripple Rate [Hz]')
            axis tight
    end
    
    % plot block transitions
    axis tight
    xticks(x_ticks)
    xticklabels(x_labels)
    xtickangle(45)
    xlabel('Time [h]')
    if markRecTrans
        hold on
        plot([x_blockTrans; x_blockTrans] * ts, ylim, '--k', 'HandleVisibility', 'off')
    end
    
    if saveFig
        mousepath = fileparts(basepaths{1});
        [~, mname] = fileparts(mousepath);
        figpath = fullfile(mousepath, 'graphics');
        figname = sprintf('%s_%s_sessions', mname, dataPreset);
        mkdir(figpath)
        figname = fullfile(figpath, figname);
        savefig(gcf, figname, 'compact')
    end
    
end

end

% EOF
