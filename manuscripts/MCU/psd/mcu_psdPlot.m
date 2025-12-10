function mcu_psdPlot(stats, psd_data, psd_cfg, varargin)

% plots psd data and statistical results from mcu_psdCBP
%
% INPUT
%   stats       output from psd_cbp
%   psd_data    cells of psd mats from mcu_psdOrg
%   psd_cfg     config from mcu_psdOrg
%   clr         rgb triplet
%   axh         handle to axis for plotting {[]}
%   ptype       'heat' or 'freq' {freq}
%   flg_log     logical, plot frequency on log scale
%   faxis       numeric frequency of psd
%
% CALLS
%   plot_stdShade
%
% 05 Jan 25 LH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'axh', [], @(x) isempty(x) || ishandle(x));
addOptional(p, 'ptype', 'freq', @(x) ismember(x, {'freq', 'heat'}));
addOptional(p, 'clr', [0.3 0.3 0.3], @isnumeric);
addOptional(p, 'flg_log', true, @islogical); 
addOptional(p, 'faxis', [], @isnumeric);

parse(p, varargin{:});
axh         = p.Results.axh;
faxis       = p.Results.faxis;
ptype       = p.Results.ptype;
clr         = p.Results.clr;
flg_log     = p.Results.flg_log;

if isempty(faxis)
    faxis = psd_cfg.faxis;
end

if iscell(psd_data)
    ngrps = length(psd_data);
else
    ngrps = 1;
end

if isempty(stats)
    stats.mask = [];
end

if isempty(psd_cfg)
    psd_cfg.grps = split(num2str(1 : ngrps));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figure if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(axh)
    fh = figure;
    set(fh, 'WindowState', 'maximized');
    set(fh, 'DefaultAxesFontSize', 16);
    axh = gca;
end
hold(axh, 'on')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(ptype, 'heat')

    % plot t-statistic (better than p-value as it shows direction of effect)
    tmap = log10(reshape(stats.stat, [length(stats.freq), size(psd_data{1}, 3)]));
    tmap = (reshape(stats.prob, [length(stats.freq), size(psd_data{1}, 3)]));

    X = 1 : size(tmap, 2);
    Y = stats.freq;

    if flg_log

        % Interpolate data onto log-spaced grid
        Y_log = logspace(log10(min(Y)), log10(max(Y)), length(Y)*10);  % 10x more points for smoother interp
        [X_orig, Y_mesh] = meshgrid(X, Y);
        [X_new, Y_new] = meshgrid(X, Y_log);
        tmap_log = interp2(X_orig, Y_mesh, tmap, X_new, Y_new);

        % Plot interpolated data
        imagesc(axh, X, Y_log, tmap_log);
        axis xy
        set(axh, 'YScale', 'log')

    else
        imagesc(axh, X, Y, tmap);
        axis xy
    end
    
    % mark significant clusters
    if ~isempty(stats.mask)
        sigMask = reshape(stats.mask, size(tmap));
        contour(axh, X, Y, sigMask, [0.5 0.5], 'k', 'LineWidth', 2);
    end

    % use diverging colormap centered at zero
    colormap(axh, brewermap([], '*RdBu'));    % or flipud(colormap('jet'))
    c = colorbar;
    c.Label.String = 't-statistic';

    % formatting
    xlabel('Time [days]')
    ylabel('Frequency [Hz]')

else    % frequency plot

    for igrp = 1 : ngrps

        clr = (igrp == 1) * clr + (igrp == 2) * clr * 0.5;

        % plot mean +/- SEM for each group
        plot_stdShade('dataMat', psd_data{igrp}, 'xVal', faxis,...
            'axh', axh, 'clr', clr, 'alpha', 0.5)
    end

    % formatting
    if flg_log
        set(axh, 'YScale', 'log')
    end
    set(axh, 'XScale', 'log')
    xlabel('Frequency [Hz]')
    ylabel('Power')

    % set y limits
    allData = cat(2, psd_data{:});
    ymax = max(mean(allData, 2, 'omitnan'));
    ymin = min(mean(allData, 2, 'omitnan'));
    ylimit = [10^floor(log10(ymin)), 10^ceil(log10(ymax))];
    set(axh, 'YLim', ylimit);

    % mark significant clusters
    if ~isempty(stats.mask)
        sigFreq = stats.freq(stats.mask);
        if ~isempty(sigFreq)
            plot(axh, sigFreq, ylimit(2) * ones(size(sigFreq)), 'k*')
        end
    end

    legend(psd_cfg.grps, 'Location', 'best')
end

end

% EOF