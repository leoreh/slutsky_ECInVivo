% mcu_fr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate fr per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uses bout times from psd (after cleaning bouts)

flgEmg = false;

% go over each mouse and analyze all experiment days
grps = [mcu_sessions('wt'), mcu_sessions('mcu')];

% go over baseline for wt vs. mcu. during experiment, only psdEmg should be
grps = {'mcu_bsl'; 'wt_bsl'};
grps = {'wt_bsl_ripp'};

% vars
if flgEmg
    vars = {'spikes'; 'psdEmg'; 'sleep_statesEmg'};
    fnameFr = 'frEmg';
else
    vars = {'spikes'; 'psd'; 'sleep_states'};
    fnameFr = 'fr';
end

% iterate
for igrp = 1 : length(grps)

    % get group basepaths
    queryStr = grps{igrp};
    grppaths = mcu_sessions(queryStr);
    nfiles = length(grppaths);
    
    % load state vars
    v = basepaths2vars('basepaths', grppaths, 'vars', vars);

    for ifile = 1 : nfiles

        % files
        basepath = grppaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)
        
        % bout time
        % Ensure v(ifile).psd and subsequent fields exist before accessing btimes
        if isfield(v(ifile), 'psd') && isfield(v(ifile).psd, 'bouts') && isfield(v(ifile).psd.bouts, 'times')
            btimes = v(ifile).psd.bouts.times;
        else
            btimes = [];
        end

        % get mfr per bout
        fr = calc_fr(v(ifile).spikes.times, 'basepath', basepath,...
            'graphics', false, 'binsize', 60, 'saveVar', fnameFr, 'forceA', true,...
            'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf],...
            'btimes', btimes);             
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASELINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Organize files
grps = {'wt_bsl'; 'mcu_bsl'};

clear grppaths
for igrp = 1 : length(grps)
    grppaths{igrp} = string(mcu_sessions(grps{igrp})');
end

% -------------------------------------------------------------------------
% FR per unit, WT vs MCU for RS vs FS 
frml = 'FR ~ Group * UnitType + (1|Mouse)';
varFld = '';

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
altClassify = 2;
fname = lme_frml2char(frml, 'rmRnd', true,...
    'sfx', ['_altClassify', num2str(altClassify)]);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BACLOFEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% load data for each group
grps = {'wt', 'mcu'}; 
idxRm = [2, 6];                         % remove bac on and off
grppaths = cell(1,length(grps));       % Initialize
for iGrp = 1 : length(grps)

    mNames = mcu_sessions(grps{iGrp});  
    mPaths = strings(length(mNames), length(mcu_sessions(mNames{1})) - length(idxRm)); % preallocate
    for imouse = 1 : length(mNames)
        tmpPaths = mcu_sessions(mNames{imouse});
        tmpPaths(idxRm) = [];
        mPaths(imouse, :) = string(tmpPaths)';
    end
    grppaths{iGrp} = mPaths;
end

% -------------------------------------------------------------------------
% FR ~ Group * Day + (1|Mouse)
frml = 'FR ~ Group * Day + (1|Mouse)'; 
varFld = '';

altClassify = 3;

% get data
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld);

% select unit
unitType = 'pINT';
iUnit = categorical({unitType}); 
plotTbl = lmeData(lmeData.UnitType == iUnit, :);

% select state
% stateType = 'High EMG';
% iState = categorical({stateType}); 
% plotTbl = plotTbl(plotTbl.State == iState, :);

% normalize
plotTbl = lme_normalize('lmeData', plotTbl, 'normVar', 'Day',...
    'groupVars', {'Group', 'UnitType', 'State'});

% run lme
contrasts = 'all'; 
contrasts =[1 : 19]; 

frmlCompare1 = 'FR ~ Group * Day + (1|Mouse)';
frmlCompare2 = 'FR ~ Group * Day + (Day|Mouse)';
lmeCfg.frml = {frmlCompare1, frmlCompare2}; 
% lmeCfg.frml = 'FR ~ Group * Day + (Day|Mouse)'; 
[lmeStats, lmeCfg] = lme_analyse(plotTbl, lmeCfg, 'contrasts', contrasts);

% add significance lines
if strcmp(unitType, 'pPYR')
    if altClassify == 2
        barIdx = {[2, 4], [2, 8], [1, 3], [1, 7]};
        barLbl = {'****', '****' '****', 'NS'};
    else
        barIdx = {[3, 4], [7, 8], [1.5, 3.5], [2, 8], [1, 7]};
        barLbl = {'NS', '****', '****', '****', '*'};
    end
    lineOffset = 0.2;
    lineGap = 0;
else
    if altClassify == 2
        barIdx = {[2, 4], [2, 8], [1, 7]};
        barLbl = {'NS', 'NS', '*'};
    else
        barIdx = {[2, 4], [2, 8], [1, 3], [1, 7]};
        barLbl = {'***', '*', 'NS', '***'};
    end
    lineOffset = 0.3;
    lineGap = -0.01;
end

% plot
hFig = lme_plot(plotTbl, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'wide');

% Update labels
hAx = gca; % Renamed
ylabel(hAx, 'Firing Rate (% BSL)')
xlabel(hAx, '')
title(hAx, unitType)
hAx.Legend.Location = 'northeast';

% Assert Size
plot_sigLines(hAx, barIdx, barLbl, 'lineGap', lineGap, 'lineOffset', lineOffset)
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300);

% save
fname = lme_frml2char(lmeCfg.frmlMdl, 'rmRnd', true, 'resNew', '',...
    'sfx', [' _', unitType, '_Norm', '_altClassify', num2str(altClassify)]);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', plotTbl, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)










% -------------------------------------------------------------------------
% FR ~ Group * UnitType * Day + (1|Mouse)
frml = 'FR ~ Group * UnitType * Day + (1|Mouse)'; 
varFld = '';

altClassify = 3;

% get data
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', true, 'varFld', varFld);

% normalize
plotTbl = lme_normalize('lmeData', lmeData, 'normVar', 'Day',...
    'groupVars', {'Group', 'UnitType', 'State'});

% run lme
contrasts = 'all'; 
contrasts =[1 : 11, 16 : 19]; 

frmlCompare1 = 'FR ~ Group * UnitType * Day + (1|Mouse)'; 
frmlCompare2 = 'FR ~ Group * UnitType * Day + (Day|Mouse)'; 
lmeCfg.frml = {frmlCompare1, frmlCompare2}; 
% lmeCfg.frml = 'FR ~ Group * Day + (Day|Mouse)'; 
[lmeStats, lmeCfg] = lme_analyse(plotTbl, lmeCfg, 'contrasts', contrasts);

% add significance lines
if strcmp(unitType, 'pPYR')
    if altClassify == 2
        barIdx = {[2, 4], [2, 8], [1, 3], [1, 7]};
        barLbl = {'****', '****' '****', 'NS'};
    else
        barIdx = {[2, 4], [2, 8], [1, 3], [1, 7]};
        barLbl = {'****', '****' '****', 'NS'};
    end
    lineOffset = 0.03;
else
    if altClassify == 2
        barIdx = {[2, 4], [2, 8], [1, 7]};
        barLbl = {'NS', 'NS', '*'};
    else
        barIdx = {[2, 4], [1, 7], [2, 8]};
        barLbl = {'NS', '*', 'NS'};
    end
    lineOffset = 0.44;
end

% plot
hFig = lme_plot(plotTbl, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'wide');

% Update labels
hAx = gca; % Renamed
ylabel(hAx, 'Firing Rate (% BSL)')
xlabel(hAx, '')
title(hAx, unitType)
hAx.Legend.Location = 'northeast';

% Assert Size
plot_sigLines(hAx, barIdx, barLbl, 'lineGap', 0, 'lineOffset', lineOffset)
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300);

% save
fname = lme_frml2char(lmeCfg.frmlMdl, 'rmRnd', true, 'resNew', '',...
    'sfx', [' _', unitType, '_Norm', '_altClassify', num2str(altClassify)]);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', plotTbl, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FR per unit, WT vs MCU across states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'FR ~ Group * State + (1|Mouse)'; % Renamed to avoid conflict


% FR per unit per bout, WT vs MCU across states 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frmlBout = 'FR ~ Group * State + (1|BoutDur) + (1|Mouse) + (1|UnitID)';
varFieldBout = 'bouts'; % Assuming varFld should be 'bouts' for per-bout analysis

% investigate relation BL and FR
if exist('frTblBout', 'var') && ~isempty(frTblBout)
    % [rMat, pMat, fhCorr] = mcu_FRvBL(frTblBout); % Assuming mcu_FRvBL is a custom function

    % limit bouts so that distribution is the same
    % [idxEq, fhEqBLen] = mcu_eqBLen(frTblBout); % Assuming mcu_eqBLen is a custom function

    % eqTbl = frTblBout(idxEq, :);
    % lmeEq = fitlme(eqTbl, lmeCfgBout.frml);
    % lme_plot(eqTbl, lmeEq)
else
    warning('frTblBout not available for BL correlation analysis.');
end


