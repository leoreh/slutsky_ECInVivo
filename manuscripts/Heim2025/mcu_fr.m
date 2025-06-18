% mcu_fr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go over each mouse and analyze all experiment days
grps = [mcu_sessions('wt'), mcu_sessions('mcu')];

% go over baseline for wt vs. mcu. during experiment, only psdEmg should be
grps = {'mcu_bsl'; 'wt_bsl'};
grps = {'wt_bsl_ripp'};

% vars
flgEmg = false;
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
% Organize Table
grps = {'wt_bsl'; 'mcu_bsl'};
lblGrp = {'Control'; 'MCU-KO'};
lblUnit = {'pPYR', 'pINT'};
vars = {'fr', 'units', 'st_metrics', 'st_brst'};

% Define variable mapping
clear varMap tblCell v tagAll tagFiles
varMap.FR = 'fr.mfr';               % Mean firing rate
varMap.BrRoy = 'st.royer';          % Mean firing rate
varMap.BrPct = 'brst.spkprct';      % Mean firing rate
varMap.UnitType = 'units.clean';    % Unit classification

% Load and organize table
for iGrp = 1 : length(grps)
    basepaths = mcu_sessions(grps{iGrp});
    nFiles = length(basepaths);
    
    v{iGrp} = basepaths2vars('basepaths', basepaths, 'vars', vars);    
        
    tagAll.Group = lblGrp{iGrp};
    tagFiles.Name = get_mname(basepaths);
    tblCell{iGrp} = v2tbl('v', v{iGrp}, 'varMap', varMap, 'tagAll',...
        tagAll, 'tagFiles', tagFiles, 'idxCol', 2);  
end
lmeData = vertcat(tblCell{:});

% Clean up
lmeData.UnitType = categorical(lblUnit(lmeData.UnitType + 1)');
lmeData = rmmissing(lmeData);
lmeData.Group = reordercats(lmeData.Group, lblGrp);
lmeData.UnitType = reordercats(lmeData.UnitType, lblUnit);

% -------------------------------------------------------------------------
% FR per unit, WT vs MCU for RS vs FS 
clear lmeCfg
frml = 'BrRoy ~ Group * UnitType + (1|Name)';
lmeCfg.contrasts = 'all';
lmeCfg.contrasts = [1 : 6]; 
lmeCfg.distribution = 'gamma';

% Fit
[lmeStats, lmeMdl] = lme_analyse(lmeData, frml, lmeCfg);

% plot
idxRow = [2 : 6];
hFig = lme_plot(lmeData, lmeMdl, 'lmeStats', lmeStats,...
    'idxRow', idxRow, 'ptype', 'bar', 'axShape', 'square');

% Update labels
hAx = gca;
ylabel(hAx, 'Firing Rate (Hz)')
xlabel(hAx, '')
title(hAx, '')
hAx.Legend.Location = 'northwest';

% Sig Lines
barIdx = {[1, 2], [3, 4]};
barLbl = {'NS', '****'};
plot_sigLines(hAx, barIdx, barLbl, 'lineGap', 0.15)
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square', 'axHeight', 300)

% Save
altClassify = 2;
fname = lme_frml2char(frml, 'rmRnd', true,...
    'sfx', ['_altClassify', num2str(altClassify)]);
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeMdl', lmeMdl)








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BACLOFEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Load and organize data
grps = {'wt', 'mcu'}; 
vars = {'fr', 'units'};

lblGrp = {'Control'; 'MCU-KO'};
lblDay = {'BSL'; 'BAC1'; 'BAC2'; 'BAC3'; 'WASH'};
idxRm = [2, 6];                     % remove bac on and off

% Define variable mapping for v2tbl
clear varMap tblCell v tagAll tagFiles
varMap.FR = 'fr.mfr';               % Mean firing rate
varMap.UnitType = 'units.clean';    % Unit classification

% Single loop over groups and mice (no day loop)
for iGrp = 1 : length(grps)
    mNames = mcu_sessions(grps{iGrp});
    
    for iMouse = 1 : length(mNames)
        % Get all paths for this mouse and remove specified days
        tmpPaths = mcu_sessions(mNames{iMouse});
        tmpPaths(idxRm) = [];
        
        % Load data for all days for this mouse at once
        v{iGrp, iMouse} = basepaths2vars('basepaths', tmpPaths, 'vars', vars);
        
        % Prepare tag structures for this mouse
        tagAll.Group = lblGrp{iGrp};
        tagAll.Name = mNames{iMouse}; 
        tagFiles.Day = lblDay(1:length(tmpPaths)); 
        
        % Create table for this mouse using new flexible approach
        tblCell{iGrp, iMouse} = v2tbl('v', v{iGrp, iMouse}, 'varMap', varMap, ...
            'tagFiles', tagFiles, 'tagAll', tagAll, 'idxCol', 2);
    end
end

% Combine all tables
lmeData = vertcat(tblCell{:});

% Clean up and organize
lblUnit = {'pPYR', 'pINT'};
lmeData.UnitType = categorical(lblUnit(lmeData.UnitType + 1)');
lmeData = rmmissing(lmeData);
lmeData.Group = reordercats(lmeData.Group, lblGrp);
lmeData.UnitType = reordercats(lmeData.UnitType, lblUnit);
lmeData.Day = reordercats(lmeData.Day, lblDay);

% -------------------------------------------------------------------------
% FR ~ Group * Day + (1|Mouse) 

% select unit
unitType = 'pINT';
iUnit = categorical({unitType}); 
plotTbl = lmeData(lmeData.UnitType == iUnit, :);

% normalize
plotTbl = tbl_transform(plotTbl, 'varNorm', 'Day',...
    'varsGrp', {'Group', 'UnitType'}, 'flgNorm', true);

% run lme
clear lmeCfg
frml = 'FR ~ Group * Day + (1|Name)';
lmeCfg.contrasts = 'all'; 
lmeCfg.distribution = 'Gamma'; 
[lmeStats, lmeMdl] = lme_analyse(plotTbl, frml, lmeCfg);

% plot
idxRow = [12 : 19];
hFig = lme_plot(plotTbl, lmeMdl, 'lmeStats', lmeStats,...
    'idxRow', idxRow, 'ptype', 'bar', 'axShape', 'wide');

% Update labels
hAx = gca;
ylabel(hAx, 'Firing Rate (% BSL)')
xlabel(hAx, '')
title(hAx, unitType)
hAx.Legend.Location = 'northeast';

% Assert Size
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300);

% save
altClassify = 3;
fname = lme_frml2char(frml, 'rmRnd', true, 'resNew', '',...
    'sfx', [' _', unitType, '_Norm', '_altClassify', num2str(altClassify)]);
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', plotTbl, 'lmeStats', lmeStats, 'lmeMdl', lmeMdl)










