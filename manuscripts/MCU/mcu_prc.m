clear params
params.winLim = [0 60 * 60];        % Analysis window [s]
params.binSize = 0.001;             % 1ms bins
params.gkHw = 0.012;                % 12ms sigma
params.winStpr = 1.0;               % 1s window
params.nShuffles = 100;            % Number of shuffles
params.spkLim = 1000;

basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
basepaths = [mcu_basepaths('mea_bac')];
nFiles = length(basepaths);
vars = {'spikes', 'st_metrics', 'fr', 'units'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

tic
for iFile = 1 : nFiles

    basepath = basepaths{iFile};
    cd(basepath);
    tic
    [prc] = prc_calc(v(iFile).spikes.times, params, 'flgSave', false);
    toc
    prc_plot(prc, 'basepath', basepath, 'flgSave', false);

end
toc





tmppaths = basepaths(1);
nFiles = length(tmppaths);
vars = {'prc', 'st_metrics', 'fr', 'units'};
v = basepaths2vars('basepaths', tmppaths, 'vars', vars);

prc = catfields([v(:).prc], 1);
units = catfields([v(:).units], 2);
st = catfields([v(:).st], 2);
fr = catfields([v(:).fr], 1);

% Burstiness
figure;
xLbl = 'Burstiness';
yLbl = 'PRC';
cnt = 1;
for iUnit = 1 : 2

    subplot(1, 2, cnt);
    unitIdx = units.clean(iUnit, :);
    x = st.royer(unitIdx);
    % x = fr.mfr(unitIdx);
    y = prc.prc0_norm(unitIdx);
    plot_scatterCorr(x, y, 'xLbl', xLbl, 'yLbl', yLbl, 'flgOtl', true);
    cnt = cnt + 1;
end
sgtitle('Burstiness')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Organize files
grps = {'wt_bsl'; 'mcu_bsl'};

clear grppaths
for igrp = 1 : length(grps)
    grppaths{igrp} = string(mcu_basepaths(grps{igrp})');
end
grppaths{1}(end) = [];
grppaths{1} = [grppaths{1}; "D:\Data\lh123\lh123_221219_094508";
    "D:\Data\lh126\lh126_230111_091208"];


% -------------------------------------------------------------------------
% FR per unit, WT vs MCU for RS vs FS
frml = 'PRC ~ Group * UnitType + (1|Mouse)';
varFld = 'prc0_norm';

% get data
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);

% run lme
contrasts = 'all';
contrasts = [1 : 6];
[lmeStats, lmeCfg] = lme_analyse(lmeData, lmeCfg, 'contrasts', contrasts);

% plot
hFig = lme_plot(lmeData, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'square');

% Update labels
hAx = gca;
ylabel(hAx, 'Firing Rate (Hz)')
xlabel(hAx, '')
title(hAx, '')
hAx.Legend.Location = 'northwest';

% add significance lines
barIdx = {[1, 2], [3, 4]};
barLbl = {'NS', '****'};
plot_sigLines(hAx, barIdx, barLbl, 'lineGap', 0.15)
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square', 'axHeight', 300)

% save
fname = lme_frml2char(frml, 'rmRnd', true,...
    'sfx', ['_altClassify', num2str(altClassify)]);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Grab data
frml = 'Burst ~ Group * UnitType + (1|Mouse)';
varFld = 'mizuseki';
varFld = 'royer';
[brstData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);
frml = 'FR ~ Group * UnitType + (1|Mouse)';
[frData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false);
frml = 'PRC ~ Group * UnitType + (1|Mouse)';
varFld = 'prc0_norm';
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);

% Find shared UnitIDs
sharedIDs = intersect(brstData.UnitID, lmeData.UnitID);
[~, sharedIdx] = ismember(sharedIDs, brstData.UnitID);

% Select only shared units and combine table
brstData = brstData(sharedIdx, :);
frData = frData(sharedIdx, :);
lmeData.Burst = brstData.Burst;
lmeData.FR = frData.FR;


% Plot
grpStr = {'Control', 'MCU-KO'};
unitStr = {'pPYR', 'pINT'};

% Burstiness
figure;
xLbl = 'Burstiness';
yLbl = 'PRC';
cnt = 1;
for iGrp = 1 : 2
    grpIdx = lmeData.Group == grpStr{iGrp};
    for iUnit = 1 : 2

        unitIdx = lmeData.UnitType == unitStr{iUnit};
        plotIdx = unitIdx & grpIdx;
        plotIdx = unitIdx;
        x = lmeData.Burst(plotIdx);
        y = lmeData.PRC(plotIdx);

        subplot(2,2,cnt);
        plot_scatterCorr(x, y, 'xLbl', xLbl, 'yLbl', yLbl, 'flgOtl', false);
        cnt = cnt + 1;
    end
end
sgtitle('Burstiness')

% Firing rate
figure;
xLbl = 'FR';
yLbl = 'PRC';
cnt = 1;
for iGrp = 1 : 2
    grpIdx = lmeData.Group == grpStr{iGrp};
    for iUnit = 1 : 2

        unitIdx = lmeData.UnitType == unitStr{iUnit};
        plotIdx = unitIdx & grpIdx;
        x = lmeData.FR(plotIdx);
        y = lmeData.PRC(plotIdx);

        subplot(2,2,cnt);
        plot_scatterCorr(x, y, 'xLbl', xLbl, 'yLbl', yLbl, 'flgOtl', true);
        cnt = cnt + 1;
    end
end
sgtitle('FR')






