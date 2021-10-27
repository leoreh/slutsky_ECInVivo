
% compares distribution of state-dependent firing rate during different
% times of the same recording or during different sessions. can compare
% bins of mu activity or su firing rate. uses box plots with diamonds
% marking the mean. relies on vars loaded from fr_sessions, session params,
% as_configFile, selecUnits.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general
saveFig = false;
grp = [1 : 4];
stateIdx = [1, 4, 5];
dataType = 'su';        % can be 'mu' (bins of sr) or 'su' (fr)
unitClass = 'int';      % plot 'int', 'pyr', or 'all'
suFlag = 1;             % plot only su or all units
frBoundries = [0 Inf];  % include only units with mean fr in these boundries

% initialize
maxY = 0;
[cfg_colors, cfg_names, ~] = as_loadConfig([]);

% select sessions. if sessionidx > 1 will ingore tstamps
sessionidx = [3];
% sessionidx = 1 : nsessions;

% select of tstamps from the same recording for comparison
assignVars(varArray, sessionidx(1))    
if length(sessionidx) == 1   
    csamps = cumsum(datInfo.nsamps) / fs;
    tstamps = [1, floor(csamps(1));...
        ceil(csamps(1)), Inf];    
    
    tstamps = [1, 6 * 60 * 60 - 1;...
        6 * 60 * 60 Inf];
    
    xLabel = ["Before"; "After"];
    figpath = session.general.basePath;
else
    figpath = fileparts(session.general.basePath);
end
cd(figpath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear frcell units
if length(sessionidx) > 1
    % compare different sessions
    for istate = 1 : length(stateIdx)
        for isession = 1 : length(sessionidx)
            assignVars(varArray, sessionidx(isession))
            switch dataType
                case 'mu'
                    frcell{istate, isession} = sr.states.fr{stateIdx(istate)}(grp, :);
                    frcell{istate, isession} = frcell{istate, isession}(:);
                case 'su'
                    units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, unitClass);
                    nunits = sum(units);
                    frcell{istate, isession} = mean(fr.states.fr{stateIdx(istate)}(units, :), 2, 'omitnan');
            end
            maxY = max([prctile(frcell{istate, isession}, 80), maxY]);
        end
    end
    
    % get x labels for dirnames
    for isession = 1 : length(sessionidx)
        [dt, ~] = guessDateTime(dirnames(sessionidx(isession)));
        xLabel{isession} = datestr(datenum(dt), 'dd/mm_HH:MM');
    end
else
    
    % compare different times from the same session
    for istate = 1 : length(stateIdx)
        for itstamps = 1 : size(tstamps, 1)
            tidx = fr.states.tstamps{stateIdx(istate)} > tstamps(itstamps, 1) &...
                fr.states.tstamps{stateIdx(istate)} < tstamps(itstamps, 2);
            switch dataType
                case 'mu'
                    frcell{istate, itstamps} = sr.states.fr{stateIdx(istate)}(grp, tidx);
                    frcell{istate, itstamps} = frcell{istate, itstamps}(:);
                case 'su'
                    units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, unitClass);
                    nunits = sum(units);
                    frcell{istate, itstamps} = mean(fr.states.fr{stateIdx(istate)}(units, tidx), 2, 'omitnan');
            end
            maxY = max([prctile(frcell{istate, itstamps}, 90), maxY]);
        end
    end
end

% gain factor between two states
clear gainfactor
stateRatio = [1, 4];    % two states for comparison
if ~isempty(stateRatio)
    
    for itime = 1 : size(frcell, 2)
        gainfactor{itime} = frcell{1, itime} - frcell{2, itime} ./...
            [max([frcell{1, itime}, frcell{2, itime}]')]';
    end
    
    fh = figure;
    gainmat = cell2nanmat(gainfactor, 1);
    plot([1 : size(gainmat, 2)], mean(gainmat, 1, 'omitnan'),...
        'kd', 'markerfacecolor', 'k')
    boxplot(gainmat, 'PlotStyle', 'traditional', 'Whisker', 3);
%     set(gca, 'YScale', 'log')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
for istate = 1 : size(frcell, 1)
    subplot(1, length(stateIdx), istate)
    hold on
    frmat = cell2nanmat(frcell(istate, :));
    plot([1 : size(frmat, 2)], mean(frmat, 1, 'omitnan'),...
        'kd', 'markerfacecolor', 'k')
    boxplot(frmat, 'PlotStyle', 'traditional', 'Whisker', 1.5);
    bh = findobj(gca, 'Tag', 'Box');
    bh = flipud(bh);
    for ibox = 1 : length(bh)
        patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
            cfg_colors{stateIdx(istate)}, 'FaceAlpha', 0.5)
    end
    xticklabels(xLabel)
    xtickangle(45)
    if istate == 1
        ylabel([dataType, ' firing rate [Hz]'])
    end
    title(cfg_names(stateIdx(istate)))
    ylim([0 ceil(maxY)])
end

if saveFig
    figpath = fullfile(figpath, 'graphics', 'sleepState');
    mkdir(figpath)
    switch dataType
        case 'mu'
            figname = fullfile(figpath, ['mu_states_times']);
        case 'su'
            figname = fullfile(figpath, [unitClass, '_states_times']);
    end
    export_fig(figname, '-tif', '-transparent', '-r300')
end

% EOF 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arrange to prism --------------------------------------------------------
clear prismData
prismData = [];
for istate = 1 : size(frcell, 1)
    x = cell2nanmat(frcell(istate, :));
    prismData = [prismData, x];
    
%     x = x(:);
%     prismData(istate, :) = x;
end
