% mcu_fr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate burstiness per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go over each mouse and analyze all experiment days
grps = [unique(get_mname(mcu_basepaths('wt'))), unique(get_mname(mcu_basepaths('mcu')))];
vars = {'spikes'};

% iterate
for igrp = 1 : length(grps)

    % get group basepaths
    queryStr = grps{igrp};
    basepaths = mcu_basepaths(queryStr);
    nfiles = length(basepaths);

    % load state vars
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);

    for ifile = 1 : nfiles

        % files
        basepath = basepaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)

        % brst (mea)
        brst = spktimes_meaBrst(v(ifile).spikes.times, 'binsize', [],...
            'isiThr', 0.02, 'minSpks', 2, 'bins', [0 Inf],...
            'saveVar', true, 'flg_force', true, 'flg_all', false);

        % spike timing metrics
        st = spktimes_metrics('spktimes', v(ifile).spikes.times, 'sunits', [],...
            'bins', [0 Inf], 'flg_force', true, 'saveVar', true, 'flg_all', false);

    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grps = {'wt_bsl'; 'mcu_bsl'};

clear grppaths
for igrp = 1 : length(grps)
    grppaths{igrp} = string(mcu_basepaths(grps{igrp})');
end


% Burst ~ Group * UnitType + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'Burst ~ Group * UnitType + (1|Mouse)';
varFld = 'royer';
varFld = 'rate';

% get data
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);

% run lme
contrasts = 'all';
contrasts = [1 : 6]; % Example of specific contrasts
[lmeStats, lmeCfg] = lme_analyse(lmeData, lmeCfg, 'contrasts', contrasts);

% plot
hFig = lme_plot(lmeData, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'square');

% add significance lines
barIdx = {};
barLbl = {};
txtY = varFld;
if strcmp(varFld , 'royer')
    txtY = 'Burstiness (a.u.)';
    barIdx = {[3, 4], [1, 2]};
    barLbl = {'NS', '***'};
elseif strcmp(varFld, 'rate')
    txtY = 'Burstiness (a.u.)';
    barIdx = {[3, 4], [1, 2]};
    barLbl = {'NS', '***'};
end

% Update labels
hAx = gca;
ylabel(hAx, txtY)
xlabel(hAx, '')
title(hAx, '')
hAx.Legend.Location = 'northeast';

% add significance lines
plot_sigLines(hAx, barIdx, barLbl, 'lineGap', 0.08)
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square', 'axHeight', 300)

% save
altClassify = 3;
fname = lme_frml2char(frml, 'rmRnd', true, 'sfx', ['_', varFld, '_altClassify', num2str(altClassify)]);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BACLOFEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for each group
grps = {'wt', 'mcu'};
idxRm = [2, 6];                         % remove bac on and off
grappaths = cell(1,length(grps));       % Initialize
for igrpBac = 1 : length(grps)

    mNames = unique(get_mname(mcu_basepaths(grps{igrpBac})));
    mPaths = strings(length(mNames), length(mcu_basepaths(mNames{1})) - length(idxRm)); % preallocate
    for imouse = 1 : length(mNames)
        basepaths = mcu_basepaths(mNames{imouse});
        basepaths(idxRm) = [];
        mPaths(imouse, :) = string(basepaths)';
    end
    grppaths{igrpBac} = mPaths;
end

% Burstiness of wt vs mcu across days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'Burst ~ Group * Day + (1|Mouse)';

varFlds = {'royer', 'bspks', 'rate'};
txtUnit = {'pPYR', 'pINT'};

for iFld = 1 : length(varFlds)
    varFld = varFlds{iFld};

    % get data
    [lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
        'flgEmg', false, 'varFld', varFld);

    % select unit
    for iUnit = 1 : length(txtUnit)
        unitType = categorical(txtUnit(iUnit));
        plotTbl = lmeData(lmeData.UnitType == unitType, :);

        % normalize
        plotTbl = lme_normalize('lmeData', plotTbl, 'normVar', 'Day',...
            'groupVars', {'Group', 'State', 'UnitType'});

        % run lme
        contrasts = 'all';
        contrasts = [1 : 11, 16 : 19]; % Example of specific contrasts
        [lmeStats, lmeCfg] = lme_analyse(plotTbl, lmeCfg, 'contrasts', contrasts);

        % add significance lines
        barIdx = {};
        barLbl = {};
        txtY = varFld;
        if strcmp(varFld , 'royer')
            txtY = 'Burstiness (% BSL)';
            if strcmp(char(unitType), 'pPYR')
                barIdx = {[2, 4], [2, 8], [1, 3], [1, 7]};
                barLbl = {'*', '*', '***', 'NS'};
                lineOffset = 0.5;
            else
                barIdx = {[1, 3], [1, 7], [2, 4]};
                barLbl = {'**', 'NS', 'NS'};
                lineOffset = 0.4;
            end
        elseif strcmp(varFld, 'bspks')
            txtY = 'Burst Spikes (% BSL)';
            if strcmp(char(unitType), 'pPYR')
                barIdx = {[1, 3], [2, 8]};
                barLbl = {'*', '****'};
            else
                barIdx = {[1, 7], [2, 8]};
                barLbl = {'*', '*'};
            end
        end

        % plot
        hFig = lme_plot(plotTbl, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'wide');
        plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300)

        % Update labels
        hAx = gca;
        ylabel(hAx, txtY)
        xlabel(hAx, '')
        title(hAx, txtUnit{iUnit})

        plot_sigLines(hAx, barIdx, barLbl,...
            'lineOffset', lineOffset, 'lineGap', 0.02, 'txtOffset', -0.01)
        plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300)
        hAx.Legend.Location = 'northeast';

        % save
        fname = lme_frml2char(frml, 'rmRnd', true, 'sfx', ['_', varFld, '_', char(unitType)]);
        lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
            'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)


    end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BURST vs. FR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'Burst ~ Group * Day + (1|Mouse)';

txtUnit = {'pPYR', 'pINT'};
iUnit = 1;

varFld = 'royer';

% get data
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);

% select unit
unitType = categorical(txtUnit(iUnit));
plotTbl = lmeData(lmeData.UnitType == unitType, :);
