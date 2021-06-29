
% compares distribution of state-dependent firing rate during different
% times of the same recording or during different sessions. can compare
% bins of mu activity or su firing rate. uses box plots with diamonds
% marking the mean. relies on vars loaded from fr_sessions, session params,
% as_configFile, selecUnits.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general
saveFig = true;
grp = [1 : 4];
stateidx = [1 : 6];
dataType = 'su';        % can be 'mu' (bins of sr) or 'su' (fr)
unitClass = 'pyr';      % plot 'int', 'pyr', or 'all'
suFlag = 1;             % plot only su or all units
frBoundries = [0 Inf];  % include only units with mean fr in these boundries

% initialize
maxY = 0;
[cfg_colors, cfg_names, ~] = as_loadConfig([]);

% selection of sessions. if sessionidx longer than one will ingore tstamps
% above
sessionidx = [1, 2, 3, 4];

% selection of timestamps from the same recording for comparison
if length(sessionidx) == 1
    assignVars(varArray, sessionidx(1))
    csamps = cumsum(datInfo.nsamps) / fs;
    tstamps = [ceil(csamps(1)), floor(csamps(2));...
        ceil(csamps(2)), floor(csamps(3));...
        ceil(csamps(3)), floor(csamps(4))];
    xLabels = ["Saline"; "K_0.5"; "K_1"];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(sessionidx) > 1
    % compare different sessions
    for istate = 1 : length(stateidx)
        for isession = 1 : length(sessionidx)
            assignVars(varArray, sessionidx(isession))
            switch dataType
                case 'mu'
                    frcell{istate, isession} = sr.states.fr{istate}(grp, :);
                    frcell{istate, isession} = frcell{istate, isession}(:);
                case 'su'
                    units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, unitClass);
                    nunits = sum(units);
                    frcell{istate, isession} = mean(fr.states.fr{istate}(units, :), 2, 'omitnan');
            end
            maxY = max([prctile(frcell{istate, isession}, 95), maxY]);
        end
    end
    xLabels = sessionDate(sessionidx);
else
    % compare different times from the same session
    for istate = 1 : length(stateidx)
        for itstamps = 1 : size(tstamps, 1)
            tidx = sr.states.tstamps{istate} > tstamps(itstamps, 1) &...
                sr.states.tstamps{istate} < tstamps(itstamps, 2);
            switch dataType
                case 'mu'
                    frcell{istate, itstamps} = sr.states.fr{istate}(grp, tidx);
                    frcell{istate, itstamps} = frcell{istate, itstamps}(:);
                case 'su'
                    units = selectUnits(spikes, cm, fr, suFlag, grp, frBoundries, unitClass);
                    nunits = sum(units);
                    frcell{istate, itstamps} = mean(fr.states.fr{istate}(units, tidx), 2, 'omitnan');
            end
            maxY = max([prctile(frcell{istate, itstamps}, 95), maxY]);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
for istate = 1 : length(stateidx)
    subplot(1, length(stateidx), istate)
    hold on
    frmat = cell2nanmat(frcell(stateidx(istate), :));
    plot([1 : size(frmat, 2)], mean(frmat, 1, 'omitnan'),...
        'kd', 'markerfacecolor', 'k')
    boxplot(frmat, 'PlotStyle', 'traditional', 'Whisker', 1.5);
    bh = findobj(gca, 'Tag', 'Box');
    bh = flipud(bh);
    for ibox = 1 : length(bh)
        patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
            cfg_colors{stateidx(istate)}, 'FaceAlpha', 0.5)
    end
    xticklabels(xLabels)
    xtickangle(45)
    if istate == 1
        ylabel([dataType, ' firing rate [Hz]'])
    end
    title(cfg_names(stateidx(istate)))
    ylim([0 ceil(maxY)])
end

if saveFig
    figpath = fullfile(basepath, 'graphics', 'sleepState');
    mkdir(figpath)
    figname = fullfile(figpath, ['mu_states_times']);
    export_fig(figname, '-tif', '-transparent', '-r300')
end
