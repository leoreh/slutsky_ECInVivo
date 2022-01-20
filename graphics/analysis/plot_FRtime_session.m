function plot_FRtime_session(varargin)

% plots firing rate of mu (sr) and su (pyr and int) across time for a
% single session (according to basepath).

% INPUT
%   basepath        recording session {pwd}
%   saveFig         logical {true}
%   grp             numeric. spike groups to plot
%   dataType        char. plot 'strd' or 'norm'.
%                   boundries specified by each row. 1st row RS 2nd row FS
%   muFlag          logical. plot multi unit (sr) activity even if fr
%                   exists {false}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'grp', [], @isnumeric);
addOptional(p, 'dataType', 'strd', @ischar);
addOptional(p, 'muFlag', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
saveFig = p.Results.saveFig;
grp = p.Results.grp;
dataType = p.Results.dataType;
muFlag = p.Results.muFlag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, basename] = fileparts(basepath);
cd(basepath)

% load data
varsFile = ["session"; "cell_metrics"; "spikes.cellinfo";...
    "fr"; "datInfo"; "sr";];
varsName = ["session"; "cm"; "spikes"; "fr"; "datInfo"; "sr"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% distribute vars
session = v.session; cm = v.cm; spikes = v.spikes; fr = v.fr;
datInfo = v.datInfo; sr = v.sr;

if ~isempty(session)
    fs = session.extracellular.sr;
end
if isempty(grp)
    if ~isempty(session)
        spkgrp = session.extracellular.spikeGroups.channels;
        grp = 1 : length(spkgrp);
    end
end

% x axis in hr
if isempty(fr) || muFlag
    ts = sr.info.binsize;
    xidx = sr.tstamps / 60 / 60;
else
    ts = fr.info.binsize;
    xidx = fr.tstamps / 60 / 60;
end

% idx of block tranisition (dashed lines)
if ~isempty(datInfo) && isfield(datInfo, 'nsamps')
    csum = cumsum(datInfo.nsamps) / fs / 60 / 60;
    tidx = csum(:);
else
    tidx = 0;
end

% units
if ~isempty(fr)
    units = selectUnits('basepath', basepath, 'grp', grp, 'saveVar', false,...
        'forceA', true);
end

switch dataType
    case 'norm'
        ytxt = 'Norm Firing Rate';
    case 'strd'
        ytxt = 'Firing Rate [Hz]';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(false)

fh = figure;
title(basename)

if isempty(fr) | muFlag
    yLimit = ceil([0 max(max(sr.(dataType)(grp, :)))]);
    hold on
    plot(xidx, sr.strd(grp, :)', 'LineWidth', 1)
    plot([tidx tidx], yLimit, '--k', 'LineWidth', 1)
    axis tight
    set(gca, 'box', 'off');
    ylabel(['Multi-unit ', ytxt])
    xlabel('Time [h]')
    lgh = legend(split(num2str(grp)));
    figname = 'sr_time';
    
else
    
    % ---------------------------------------------------------------------
    % individual cells on a log scale
    sb1 = subplot(3, 1, 1);
    if sum(units(1, :)) > 0
        title(basename)
        hold on
        % rs
        ph = plot(xidx, fr.(dataType)(units(1, :), :), 'b', 'LineWidth', 1);
        alphaIdx = linspace(1, 0.2, length(ph));
        clrIdx = linspace(0.2, 0.6, length(ph));
        [~, mfr_order] = sort(fr.mfr(units(1, :)));
        for iunit = 1 : length(ph)
            ph(iunit).Color(1) = clrIdx(mfr_order(iunit));
            ph(iunit).Color(2) = alphaIdx(mfr_order(iunit));
            ph(iunit).Color(4) = alphaIdx(mfr_order(iunit));
        end
        set(gca, 'YScale', 'log')
        plot([tidx tidx], ylim, '--k')
        axis tight
        xlabel('Time [h]')
        ylabel(ytxt)
        set(gca, 'box', 'off')
        legend(sprintf('RS = %d su', sum(units(1, :))))
    end
    
    % fs
    sb2 = subplot(3, 1, 2);
    if sum(units(2, :)) > 0
        hold on
        ph = plot(xidx, fr.(dataType)(units(2, :), :), 'r', 'LineWidth', 1);
        alphaIdx = linspace(1, 0.2, length(ph));
        clrIdx = linspace(0.2, 0.6, length(ph));
        [~, mfr_order] = sort(fr.mfr(units(2, :)));
        for iunit = 1 : length(ph)
            ph(iunit).Color(3) = clrIdx(mfr_order(iunit));
            ph(iunit).Color(2) = alphaIdx(mfr_order(iunit));
            ph(iunit).Color(4) = alphaIdx(mfr_order(iunit));
        end
        set(gca, 'YScale', 'log')
        plot([tidx tidx], ylim, '--k')
        axis tight
        xlabel('Time [h]')
        ylabel(ytxt)
        set(gca, 'box', 'off')
        legend(sprintf('FS = %d su', sum(units(2, :))));
    end
    
    % ---------------------------------------------------------------------
    % mean per cell class on a linear scale
    sb3 = subplot(3, 1, 3);
    hold on
    plot(xidx, mean(fr.(dataType)(units(1, :), :), 1, 'omitnan'), 'b', 'LineWidth', 2)
    ylabel(['RS ' ytxt])
    axis tight
    yyaxis right
    plot(xidx, mean(fr.(dataType)(units(2, :), :), 1, 'omitnan'), 'r', 'LineWidth', 2)
    axis tight
    ylabel(['FS ' ytxt])
    yLimit = ylim;
    plot([tidx tidx], yLimit, '--k')
    axis tight
    ax = gca;
    set(ax.YAxis(1), 'color', 'b')
    set(ax.YAxis(2), 'color', 'r')
    xlabel('Time [h]')
    set(gca, 'box', 'off')
    linkaxes([sb1, sb2, sb3], 'x')
    figname = 'fr_time';
    
end

if saveFig
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, [basename, '_', figname]);
    export_fig(figname, '-jpg', '-transparent', '-r300')
end

end

% EOF