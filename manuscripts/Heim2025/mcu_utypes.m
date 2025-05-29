% mcu_cellCalss


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RE-ANALYZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all files in study
basepaths = mcu_sessions('all');
nPaths = length(basepaths);

% vars
vars = {'spikes'};

% load state vars
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

for iPath = 1 : nPaths

    % files
    basepath = basepaths{iPath};
    [~, basename] = fileparts(basepath);
    cd(basepath)

    % get mfr per bout
    % fr = calc_fr(v(iPath).spikes.times, 'basepath', basepath,...
    %     'graphics', false, 'binsize', 60, 'saveVar', 'fr', 'forceA', true,...
    %     'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf],...
    %     'btimes', []);

    % calculate waveform metrices
    swv = spkwv_metrics('basepath', basepath, 'flgSave', true, 'flgForce', true);

end


% Get only WT basepaths
mNames = mcu_sessions('wt');
clear mPaths
for iMouse = 1 : length(mNames)
    mPaths(iMouse, :) = string(mcu_sessions(mNames{iMouse}))';
end
basepaths = mPaths(:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reclassify units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fetTbl = utypes_features('basepaths', basepaths, 'flgPlot', true);

for iAlt = 1 : 3
    
    altClassify = iAlt;
    unitType = utypes_classify('basepaths', basepaths, 'altClassify', altClassify,...
        'flgSave', true);
    nPyr(iAlt) = sum(unitType == 1);

    fh = figure;
    swvFld = 'tpSlope';
    stFld = 'lidor';
    hAx = subplot(1, 2, 1); hold on
    plot_utypes('basepaths', basepaths, 'flgRaw', false,...
        'plotType', 'scatter3', 'swvFld', swvFld, 'stFld', stFld,...
        'unitIdx', unitType, 'hAx', hAx)

    hAx = subplot(1, 2, 2); hold on
    plot_utypes('basepaths', basepaths, 'flgRaw', false,...
        'plotType', 'wv', 'swvFld', swvFld, 'stFld', stFld,...
        'unitIdx', unitType, 'hAx', hAx)

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST FR RESTULS FROM VARIOUS CLASSIFICATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Baseline

grps = {'wt_bsl'; 'mcu_bsl'};

% Preload data
vars = {'fr', 'units'};
clear vCell
for iGrp = 1 : length(grps)
    basepaths = string(mcu_sessions(grps{iGrp})');
    vCell{iGrp}{1} = basepaths2vars('basepaths', basepaths, 'vars', vars);
end

% Initialize empty table before the loop
tblVars = ["Classification",...
    "CTL:pPYR (#)", "MCU:pPYR (#)",...
    "CTL:pINT (#)", "MCU:pINT (#)"...
    "CTL:pPYR (MFR)", "MCU:pPYR (MFR)",...
    "CTL:pINT (MFR)", "MCU:pINT (MFR)"]';
tblStats = table('Size', [0, length(tblVars)], 'VariableTypes', ...
    {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', tblVars);

% Go over all classification options
for iAlt = 1 : 3
    
    altClassify = iAlt;   
    uFld = ['clean', num2str(altClassify)];

    % Switch default units classification scheme
    clear grppaths
    for iGrp = 1 : length(grps)
        grppaths{iGrp} = string(mcu_sessions(grps{iGrp})');
        for iPath = 1 : length(grppaths{iGrp})
            vCell{iGrp}{1}(iPath).units.clean = vCell{iGrp}{1}(iPath).units.(uFld);
        end
    end

    % lme
    frml = 'FR ~ Group * UnitType + (1|Mouse)';
    [lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
        'flgEmg', false, 'varFld', '', 'vCell', vCell);

    [lmeStats, lmeCfg] = lme_analyse(lmeData, lmeCfg, 'contrasts', 'all');
    hFig = lme_plot(lmeData, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'square');

    % Get all statistics in one groupsummary call
    stats = groupsummary(lmeData, ["Group", "UnitType"], "mean", "FR",...
        'IncludeMissingGroups', false, 'IncludeEmptyGroups', false);

    % Create row for current classification
    tblStats = [tblStats;...
        table(altClassify,...
        stats.GroupCount(1), stats.GroupCount(3),...                                        % pPYR counts
        stats.GroupCount(2), stats.GroupCount(4),...                                        % pINT counts
        stats.mean_FR(1), stats.mean_FR(3),...                                              % pPYR FR
        stats.mean_FR(2), stats.mean_FR(4),...                                              % pINT FR
        'VariableNames', tblVars)];

end



% -------------------------------------------------------------------------
% Baclofen

                
% Organize sessions
grps = {'wt', 'mcu'}; 
nGrp = length(grps);
idxRm = [2, 6];    
clear grppaths mPaths
for iGrp = 1 : nGrp
    mNames = mcu_sessions(grps{iGrp});
    for iMouse = 1 : length(mNames)
        tmpPaths = mcu_sessions(mNames{iMouse});
        tmpPaths(idxRm) = [];
        mPaths(iMouse, :) = string(tmpPaths)';
    end
    grppaths{iGrp} = mPaths;
end

% Preload data
vars = {'fr', 'units'};
vCell = cell(1, nGrp);
for iGrp = 1:nGrp
    [~, nDays] = size(grppaths{iGrp});
    vCell{iGrp} = cell(1, nDays);

    for iday = 1:nDays
        basepaths = grppaths{iGrp}(:, iday);
        vCell{iGrp}{iday} = basepaths2vars('basepaths', basepaths, 'vars', vars);
    end
end

% Initialize empty table before the loop
tblVars = [...
    "CTL:pPYR Recovery (%BSL)", "MCU:pPYR Recovery (%BSL)",...
    "CTL:pINT Recovery (%BSL)", "MCU:pINT Recovery (%BSL)"]';
tblStats2 = table('Size', [0, length(tblVars)], 'VariableTypes', ...
    {'double', 'double', 'double', 'double'}, ...
    'VariableNames', tblVars);

% Go over all classification options
for iAlt = 1 : 3
    
    altClassify = iAlt;   
    uFld = ['clean', num2str(altClassify)];

    % Switch default units classification scheme
    for iGrp = 1 : nGrp
        for iDay = 1 : 5
            nMice = length(vCell{iGrp}{iDay});
            for iMouse = 1 : nMice
                vCell{iGrp}{iDay}(iMouse).units.clean = vCell{iGrp}{iDay}(iMouse).units.(uFld);
            end
        end
    end

    % lme
    frml = 'FR ~ Group * Day + (1|Mouse)';
    [lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
        'flgEmg', false, 'varFld', '', 'vCell', vCell);

    % normalize
    lmeData = lme_normalize('lmeData', lmeData, 'normVar', 'Day',...
        'groupVars', {'Group', 'UnitType'});
    
    % select unit type
    plotTbl = lmeData(lmeData.UnitType == 'pINT', :);

    [lmeStats, lmeCfg] = lme_analyse(plotTbl, lmeCfg, 'contrasts', 'all');
    hFig = lme_plot(plotTbl, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'square');
    
    % Get all statistics in one groupsummary call
    stats = groupsummary(lmeData, ["Group", "UnitType", "Day"], "mean", "FR",...
        'IncludeMissingGroups', false, 'IncludeEmptyGroups', false);

    % Create row for current classification
    tblStats2 = [tblStats2;...
        table(...
        stats.mean_FR(4), stats.mean_FR(14),...                                      
        stats.mean_FR(9), stats.mean_FR(19),...                                                                       
        'VariableNames', tblVars)];

end

tblStats = [tblStats, tblStats2];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Formula
frml = 'SWV ~ Group * UnitType + (1|Mouse)';

% organize for lme
varFld = 'tp';
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all files in study
basepaths = mcu_sessions('all');
nPaths = length(basepaths);

altClassify = 3;
swvFld = 'tp';
stFld = 'lidor';
hAx = plot_utypes('basepaths', basepaths, 'flgRaw', false,...
    'plotType', 'scatter', 'swvFld', swvFld, 'stFld', stFld,...
    'unitIdx', altClassify);
hFig = gcf;
xlabel(hAx, 'Trough-to-Peak (ms)')
ylabel(hAx, 'Burstiness (a.u.)')
axis tight
ylim([-1 1])
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square', 'axHeight', 300)

% Save
fname = ['UnitTypes~Scatter_Alt', num2str(altClassify)];
lme_save('fh', hFig, 'fname', fname, 'frmt', {'fig', 'svg'},...
    'lmeData', table(), 'lmeStats', table(), 'lmeCfg', []) 


hAx = plot_utypes('basepaths', basepaths, 'flgRaw', false,...
    'plotType', 'wv', 'swvFld', swvFld, 'stFld', stFld,...
    'unitIdx', altClassify);
hFig = gcf;
xlabel(hAx, 'Time (ms)')
ylabel(hAx, 'Amplitude (a.u.)')
yticks(hAx, [])
ylim([-0.75, 0.1])
xlim([-0.5, 0.5])
hAx.Legend.Location = 'southwest';
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'tall', 'axHeight', 300)
% hAx.Legend.FontSize = 16;
hAx.Legend.ItemTokenSize = [14 10];

% Save
fname = ['UnitTypes~Wv_Alt', num2str(altClassify)];
lme_save('fh', hFig, 'fname', fname, 'frmt', {'fig', 'svg'},...
    'lmeData', table(), 'lmeStats', table(), 'lmeCfg', []) 
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create copy of files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get all files in study
basepaths = mcu_sessions('all');  
nPaths = length(basepaths);

fNames = {...
    '*.acceleration.mat',...
    '*.datInfo*',...
    '*.cell_metrics*',...
    '*.fr.mat',...
    '*.frEmg.mat',...
    '*.nrs',...
    '*.ripp.*',...
    '*.session.mat*',...
    '*.spikes.*',...
    '*.spktimes.mat',...
    '*.sr.mat',...
    '*.st_*',...
    '*.swv_metrics*',...
    '*.swv_raw*',...
    '*.units.mat',...
    '*.xml',...
    '*.sleep_states*',...
    '*.sleep_labelsMan*',...
    % '*.lfp',...
    % '*.sleep_*',...
    };

for iPath = 1 : nPaths
    basepath = basepaths{iPath};
    cp_basepath('fNames', fNames, 'basepath', basepath,...
        'overwrite', false, 'verbose', true)
end
