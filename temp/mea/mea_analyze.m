function mea_analyze(varargin)

% recieves a mea struct (see mea_orgNex.m) and runs various analyses on
% spike waveform and spike times.
%
% INPUT:
%   mea         struct. if empty will be loaded from basepath.  
%   basepath    char. path to where mea struct exists. will be used for
%               saveing all other files.
%   winBL       2-element vector describing the window for baseline
%               calculations [sec]. {[0 Inf]}.
%   saveVar     logical. save variables {true}.
%   graphics    logical. plot graphics or not {true}. 
%   forceA      logical. reanalyze even if struct file exists
%
% DEPENDENCIES
%
% TO DO lIST
%
% 23 dec 21 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'mea', []);
addOptional(p, 'basepath', pwd, @ischar);
addOptional(p, 'winBL', [0 Inf], validate_win);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'forceA', false, @islogical);

parse(p, varargin{:})
mea         = p.Results.mea;
basepath    = p.Results.basepath;
winBL       = p.Results.winBL;
saveVar     = p.Results.saveVar;
graphics    = p.Results.graphics;
forceA      = p.Results.forceA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load struct
cd(basepath)
[~, basename] = fileparts(basepath);
if isempty(mea)
    load(fullfile(basepath, [basename, '.mea.mat']))
end

% gen params
fs = mea.info.fs;

% spike timing metrics
st = spktimesMetrics('basepath', basepath, 'spktimes', mea.spktimes,...
    'winCalc', [], 'saveVar', saveVar, 'fs', fs, 'forceA', false);

% spike waveform metrics
swv = spkwvMetrics('wv', mea.wv, 'basepath', basepath, 'fs', fs,...
    'saveVar', saveVar, 'forceA', false);

% firing rate
binsize = 60;
fr = firingRate(mea.spktimes, 'basepath', basepath,...
    'graphics', false, 'binsize', binsize, 'saveVar', saveVar,...
    'smet', 'GK', 'winBL', winBL, 'winCalc', [0, Inf], 'forceA', forceA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classify
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CE thresholds
nunits = length(fr.mfr);
cm.putativeCellType = repmat({'Pyramidal Cell'}, 1, nunits);
cm.putativeCellType(swv.tp < 0.425) = {'Narrow Interneuron'};
cm.putativeCellType(st.tau_rise > 6 & swv.tp < 0.425) = {'Wide Interneuron'};
cell_metrics = cm;
save(fullfile(basepath, [basename, '.cell_metrics.cellinfo.mat']), 'cell_metrics')

% select units based on criterion 
units(1, :) = selectUnits([], cm, fr, 0, [], [], 'pyr');
units(2, :) = selectUnits([], cm, fr, 0, [], [], 'int');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    
    % ---------------------------------------------------------------------
    % classification    
    tp = swv.tp;
    spkw = swv.spkw;
    royer = st.royer;
    lidor = st.lidor;
    tau_rise = st.tau_rise;
    mfr = normalize(fr.mfr, 'range', [0.1 1]);    
    
    fh = figure;
    subplot(1, 2, 1)
    sh = scatter(tp(units(1, :)), royer(units(1, :)),...
        mfr(units(1, :)) * 3000, 'b', '.');
    hold on
    sh = scatter(tp(units(2, :)), royer(units(2, :)),...
        mfr(units(2, :)) * 3000, 'r', '.');
    set(gca, 'yscale', 'log')
    xlabel('Trough to Peak [ms]')
    ylabel('Burstiness (royer)')
    
    subplot(1, 2, 2)
    sh = scatter(spkw(units(1, :)), tau_rise(units(1, :)),...
        mfr(units(1, :)) * 3000, 'b', '.');
    hold on
    sh = scatter(spkw(units(2, :)), tau_rise(units(2, :)),...
        mfr(units(2, :)) * 3000, 'r', '.');
    set(gca, 'yscale', 'log')
    xlabel('Spike Width [ms]')
    ylabel('Tau Rise')        
    
    % save
    figname = fullfile(figpath, [basename, '_cellClass']);
    export_fig(figname, '-jpg', '-transparent', '-r300')
    
    % ---------------------------------------------------------------------
    % figure per unit
    sunits = randperm(nunits, 20);
    for iunit = sunits
        
        fh = figure('Visible', 'off');
        txt = [basename, '_unit#', num2str(iunit)];
        title(txt)
        
        % waveform
        subplot(2, 3, 1)
        x_val = [1 : size(mea.wv, 2)] / fs * 1000;
        wv = mea.wv(iunit, :);
        wv_std = mea.wv_std(iunit, :);
        if strcmp(cm.putativeCellType{iunit}, 'Pyramidal Cell')
            clr = 'b';
        elseif strcmp(cm.putativeCellType{iunit}, 'Narrow Interneuron')
            clr = 'r';
        elseif strcmp(cm.putativeCellType{iunit}, 'Wide Interneuron')
            clr = 'g';
        end            
        plot(x_val, mea.wv(iunit, :), clr, 'LineWidth', 2)
        patch([x_val, flip(x_val)], [wv + wv_std, flip(wv - wv_std)],...
            clr, 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
        xlabel('Time [ms]')
        ylabel('Voltage [mV]')
        subtitle('Waveform')

        % fr vs. time
        subplot(2, 3, 2)
        plot(fr.tstamps / 60 / 60 * 6, fr.strd(iunit, :), 'k')
        xlabel('Time [h]')
        ylabel('Firing Rate [Hz]')
        box off
        axis tight
        subtitle('FR vs. time')

        % fr distribution
        subplot(2, 3, 3)
        hh = histogram(fr.mfr, 49, 'Normalization', 'probability');
        hh.FaceColor = 'k';
        hh.EdgeColor = 'none';
        xlabel('Firing Rate [Hz]')
        ylabel('Probability')
        subtitle('FR distribution')

        % acg narrow
        subplot(2, 3, 4)
        bh = bar(st.info.acg_narrow_tstamps,...
            squeeze(st.acg_narrow(:, 1, iunit)), 'BarWidth', 1);
        ylabel('Rate [Hz]')
        xlabel('Time [ms]')
        bh.FaceColor = 'k';
        bh.EdgeColor = 'none';
        box off
        axis tight
        subtitle('ACG narrow')

        % acg wide
        subplot(2, 3, 5)
        bh = bar(st.info.acg_wide_tstamps,...
            squeeze(st.acg_wide(:, 1, iunit)), 'BarWidth', 1);
        ylabel('Rate [Hz]')
        xlabel('Time [ms]')
        bh.FaceColor = 'k';
        bh.EdgeColor = 'none';
        box off
        axis tight
        subtitle('ACG wide')

        % isi histogram
        subplot(2, 3, 6)
        x_max = 0.1;
        isi = diff(mea.spktimes{iunit});
        hh = histogram(isi(isi < x_max), 100, 'Normalization', 'probability');
        hh.FaceColor = 'k';
        hh.EdgeColor = 'none';
        xlim([0 x_max])
        xlabel('Time [ms]')
        ylabel('Probability')
        subtitle('ISI histogram')
        
        % save
        figname = fullfile(figpath, txt);
        export_fig(figname, '-jpg', '-transparent', '-r300')
        
    end
    
    
    % ---------------------------------------------------------------------
    % firing rate vs. time
    fh = figure;
    sb1 = subplot(3, 1, 1);
    title(basename)
    hold on
    
    xidx = fr.tstamps / 60 / 60 * 6;
    
    % rs
    ph = plot(xidx, fr.strd(units(1, :), :), 'b', 'LineWidth', 1);
    alphaIdx = linspace(1, 0.2, length(ph));
    clrIdx = linspace(0.2, 0.6, length(ph));
    [~, mfr_order] = sort(fr.mfr(units(1, :)));
    for iunit = 1 : length(ph)
        ph(iunit).Color(1) = clrIdx(mfr_order(iunit));
        ph(iunit).Color(2) = alphaIdx(mfr_order(iunit));
        ph(iunit).Color(4) = alphaIdx(mfr_order(iunit));
    end
    set(gca, 'YScale', 'log')
    axis tight
    xlabel('Time [h]')
    ylabel('Firing rate [log(Hz)]')
    set(gca, 'box', 'off')
    legend(sprintf('RS = %d su', sum(units(1, :))))
    
    % fs
    sb2 = subplot(3, 1, 2);
    hold on
    ph = plot(xidx, fr.strd(units(2, :), :), 'r', 'LineWidth', 1);
    alphaIdx = linspace(1, 0.2, length(ph));
    clrIdx = linspace(0.2, 0.6, length(ph));
    [~, mfr_order] = sort(fr.mfr(units(2, :)));
    for iunit = 1 : length(ph)
        ph(iunit).Color(3) = clrIdx(mfr_order(iunit));
        ph(iunit).Color(2) = alphaIdx(mfr_order(iunit));
        ph(iunit).Color(4) = alphaIdx(mfr_order(iunit));
    end
    set(gca, 'YScale', 'log')
    axis tight
    xlabel('Time [h]')
    ylabel('Firing rate [log(Hz)]')
    set(gca, 'box', 'off')
    legend(sprintf('FS = %d su', sum(units(2, :))));
    
    % mean per cell class on a linear scale
    sb3 = subplot(3, 1, 3);
    yLimit = [0 ceil(max(mean(fr.strd(units(2, :), :), 'omitnan')))];
    hold on
    plot(xidx, mean(fr.strd(units(1, :), :), 'omitnan'), 'b', 'LineWidth', 2)
    plot(xidx, mean(fr.strd(units(2, :), :), 'omitnan'), 'r', 'LineWidth', 2)
    axis tight
    xlabel('Time [h]')
    ylabel('Firing rate [Hz]')
    set(gca, 'box', 'off')
    linkaxes([sb1, sb2, sb3], 'x')
    figname = 'fr_time';
    
    % save
    figname = fullfile(figpath, figname);
    export_fig(figname, '-jpg', '-transparent', '-r300')
    
end