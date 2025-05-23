% mcu_fr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate burstiness per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go over each mouse and analyze all experiment days
grps = [mcu_sessions('wt'), mcu_sessions('mcu')];
vars = {'spikes'};

% iterate
for igrp = 1 : length(grps)

    % get group basepaths
    queryStr = grps{igrp};
    basepaths = mcu_sessions(queryStr);
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
    grppaths{igrp} = string(mcu_sessions(grps{igrp})');
end


% Burst ~ Group + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'Burst ~ Group * UnitType + (1|Mouse)';

% organize for lme
varField = 'lvr';
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', varField, 'vCell', {});

% run lme
contrasts = 'all';
contrasts = [1:5, 6,9]; 
[lmeStats, lmeCfg] = lme_analyse(lmeData, lmeCfg, 'contrasts', contrasts);

% plot
fh = lme_plot(lmeData, lmeCfg.lmeMdl, 'ptype', 'bar'); 
axh = gca;
ylabel(axh, 'Burstiness Index (a.u.)', 'FontSize', 20)
xlabel(axh, '', 'FontSize', 20)
axh.XAxis.FontSize = 20;
title(axh, '')
axh.Legend.Location = 'northeast';

fname = lme_frml2char(lmeCfg.frmlMdl, 'rmRnd', true, 'sfx', '_Royer');

% save
lme_save('fh', fh, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis during baclofen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for each group
grps = {'wt', 'mcu'}; 
idxRm = [2, 6];                         % remove bac on and off
grappaths = cell(1,length(grps));       % Initialize
for igrpBac = 1 : length(grps)

    mNames = mcu_sessions(grps{igrpBac});  
    mPaths = strings(length(mNames), length(mcu_sessions(mNames{1})) - length(idxRm)); % preallocate
    for imouse = 1 : length(mNames)
        basepaths = mcu_sessions(mNames{imouse});
        basepaths(idxRm) = [];
        mPaths(imouse, :) = string(basepaths)';
    end
    grppaths{igrpBac} = mPaths;
end

% Burstiness of wt vs mcu across days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize for lme
frml = 'Burst ~ Group * Day + (1|Mouse)';

% organize for lme
varField = 'spkprct';
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', varField, 'vCell', {});

% select unit
unitType = 'pPV';
iUnit = categorical({unitType}); 
plotTbl = lmeData(lmeData.UnitType == iUnit, :);

% normalize
plotTbl = lme_normalize('lmeData', plotTbl, 'normVar', 'Day',...
    'groupVars', {'Group', 'State'});

% run lme
contrasts = 'all'; 
contrasts = [1 : 11, 16 : 19]; 

frmlCompare1 = 'Burst ~ Group * Day + (1|Mouse)';
frmlCompare2 = 'Burst ~ Group * Day + (Day|Mouse)';
lmeCfg.frml = {frmlCompare1, frmlCompare2}; 

[lmeStats, lmeCfg] = lme_analyse(plotTbl, lmeCfg, 'contrasts', contrasts);

% plot
hndFig = lme_plot(plotTbl, lmeCfg.lmeMdl, 'ptype', 'bar', 'figShape', 'wide');

% Update labels
axh = gca; % Renamed
ylabel(axh, 'Burst Spikes (%)', 'FontSize', 20)
xlabel(axh, '', 'FontSize', 16)
axh.XAxis.FontSize = 20;
title(axh, unitType)
axh.Legend.Location = 'northwest';
% ylim([0 130])
fname = lme_frml2char(lmeCfg.frmlMdl, 'rmRnd', true, 'sfx', [' _', unitType, '_Prcnt']);

% save
lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', plotTbl, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)
