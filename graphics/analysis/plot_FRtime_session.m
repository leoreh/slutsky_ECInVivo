function plot_FRtime_session(varargin)

% plots firing rate of mu (sr) and su (pyr and int) across time for a
% single session (according to basepath).

% INPUT
%   basepath        recording session {pwd}
%   units           numeric 2 x n of indices of units to plot
%   saveFig         logical {true}
%   dataType        char. plot 'strd' or 'norm'.
%                   boundries specified by each row. 1st row RS 2nd row FS
%   muFlag          logical. plot multi unit (sr) activity even if fr
%                   exists {false}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'basepath', pwd);
addParameter(p, 'units', []);
addParameter(p, 'saveFig', true, @islogical);
addParameter(p, 'grp', [], @islogical);
addParameter(p, 'dataType', 'strd', @ischar);
addParameter(p, 'muFlag', false, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
units       = p.Results.units;
saveFig     = p.Results.saveFig;
grp         = p.Results.grp;
dataType    = p.Results.dataType;
muFlag      = p.Results.muFlag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, basename] = fileparts(basepath);
cd(basepath)

% load data
varsFile = ["session"; "fr.mat"; "datInfo"; "sr.mat";];
varsName = ["session"; "fr"; "datInfo"; "sr"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% distribute vars
% assignVars = structvars(v);
session = v.session;     
fr = v.fr;               
datInfo = v.datInfo;     
sr = v.sr; 

if ~isempty(session)
    fs = session.extracellular.sr;
end

if isempty(grp)
    grp = 1 : session.extracellular.nSpikeGroups;
end

% x axis in hr
if isempty(fr) || muFlag
    xidx = sr.tstamps / 60 / 60;
else
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
    if isempty(units)
        units = selectUnits('basepath', basepath, 'grp', [], 'saveVar', false,...
            'forceA', false);
        units = units.idx;
    end
end

switch dataType
    case 'norm'
        ytxt = 'Norm Firing Rate';
    case 'strd'
        ytxt = 'Firing Rate [Hz]';
end

unitChar = {'RS', 'FS'};
unitClr = {'b', 'r'};

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
    for iunit = 1 : size(units, 1)
        sb{iunit} = subplot(3, 1, iunit);
        if sum(units(1, :)) < 0
            continue
        end
        title(basename)
        hold on
        ph = plot(xidx, fr.(dataType)(units(iunit, :), :),...
            unitClr{iunit}, 'LineWidth', 1);
        alphaIdx = linspace(1, 0.2, length(ph));
        clrIdx = linspace(0.2, 0.6, length(ph));
        [~, mfr_order] = sort(fr.mfr(units(iunit, :)));
        for icell = 1 : length(ph)
            ph(icell).Color(iunit) = clrIdx(mfr_order(icell));
            ph(icell).Color(2) = alphaIdx(mfr_order(icell));
            ph(icell).Color(4) = alphaIdx(mfr_order(icell));
        end
        set(gca, 'YScale', 'log')
        plot([tidx tidx], ylim, '--k')
        axis tight
        xlabel('Time [h]')
        ylabel(ytxt)
        set(gca, 'box', 'off')
        legend(sprintf('%s = %d su',...
            unitChar{iunit}, sum(units(iunit, :))))
    end    
  
    % ---------------------------------------------------------------------
    % mean per cell class on a linear scale
    sb3 = subplot(3, 1, 3);
    hold on
    data = fr.(dataType)(units(1, :), :);
    data(~isfinite(data)) = nan;
    data = mean(data, 1, 'omitnan');
    plot(xidx, data, 'b', 'LineWidth', 2)
    axis tight
    
    % create 2 axes if raw fr
    if strcmp(dataType, 'strd')
        ylabel(['RS ' ytxt])
        yyaxis right
        ylabel(['FS ' ytxt])
        ax = gca;
        set(ax.YAxis(1), 'color', 'b')
        set(ax.YAxis(2), 'color', 'r')
    end
    
    data = fr.(dataType)(units(2, :), :);
    data(~isfinite(data)) = nan;
    data = mean(data, 1, 'omitnan');
    plot(xidx, data, 'r', 'LineWidth', 2)
    axis tight
    
    yLimit = ylim;   
    plot([tidx tidx], yLimit, '--k')
    xlabel('Time [h]')
    set(gca, 'box', 'off')
    linkaxes([sb{1}, sb{2}, sb3], 'x')
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