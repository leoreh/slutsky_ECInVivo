
function plot_FRtime_session(varargin)

% plots firing rate of mu (sr) and su (pyr and int) across time for a
% single session (according to basepath).

% INPUT
%   basepath        recording session {pwd}
%   saveFig         logical {true}
%   grp             numeric. spike groups to plot
%   frBoundries     2 x 2 mat. include only units with mfr within the
%                   boundries specified by each row. 1st row RS 2nd row FS
%   muFlag          logical. plot multi unit (sr) activity even if fr
%                   exists {false}

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'grp', [1 : 4], @isnumeric);          
addOptional(p, 'frBoundries', [0.2 Inf; 0.2 Inf], @isnumeric);          
addOptional(p, 'muFlag', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
saveFig = p.Results.saveFig;
grp = p.Results.grp;
frBoundries = p.Results.frBoundries;
muFlag = p.Results.muFlag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mousepath, basename] = fileparts(basepath);
cd(basepath)
guessDateTime(basename)

suFlag = 1;                     % plot only su or all units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data

varArray = getSessionVars('dirnames', {basename}, 'mousepath', mousepath,...
    'sortDir', false);
assignVars(varArray, 1)

fs = session.extracellular.sr;

% x axis in hr
ts = sr.binsize;
xidx = [1 : length(sr.strd)] / ts;

% idx of block tranisition (dashed lines)
if ~isempty(datInfo)
    csum = cumsum(datInfo.nsamps) / fs / 60 / 60;
    tidx = csum(:);
else
    tidx = 0;
end

% units
if ~isempty(fr)
    if length(fr.strd) ~= length(sr.strd)
        warning('check length of mu and su firing rate')
    end
    
    clear units
    units(1, :) = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'pyr');
    units(2, :) = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, 'int');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(false)

fh = figure;
title(basename)

if isempty(fr) | muFlag
    yLimit = ceil([0 max(max(sr.strd(grp, :)))]);
    hold on
    plot(xidx, sr.strd(grp, :)')
    plot([tidx tidx], yLimit, '--k')
    axis tight
    set(gca, 'box', 'off');
    ylabel('Multi-unit firing rate [Hz]')
    xlabel('Time [h]')
    lgh = legend(split(num2str(grp)));
    figname = 'sr_time';
    
else
    
    % ---------------------------------------------------------------------
    % individual cells on a log scale    
    sb1 = subplot(3, 1, 1);
    title(basename)
    hold on
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
    plot([tidx tidx], ylim, '--k')
    axis tight
    xlabel('Time [h]')
    ylabel('RS firing rate [Hz]')
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
    plot([tidx tidx], ylim, '--k')
    axis tight
    xlabel('Time [h]')
    ylabel('FS firing rate [Hz]')
    set(gca, 'box', 'off')
    legend(sprintf('FS = %d su', sum(units(2, :))));

    % ---------------------------------------------------------------------
    % mean per cell class on a linear scale 
    sb3 = subplot(3, 1, 3);
    yLimit = [0 ceil(max(mean(fr.strd(units(2, :), :), 'omitnan')))];
    hold on
    plot(xidx, mean(fr.strd(units(1, :), :), 'omitnan'), 'b', 'LineWidth', 2)
    plot(xidx, mean(fr.strd(units(2, :), :), 'omitnan'), 'r', 'LineWidth', 2)
    plot([tidx tidx], yLimit, '--k')
    axis tight
    xlabel('Time [h]')
    ylabel('Single unit firing rate [Hz]')
    set(gca, 'box', 'off')
    linkaxes([sb1, sb2, sb3], 'x')    
    figname = 'fr_time';

end

if saveFig
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, [basename, '_', figname]);
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF