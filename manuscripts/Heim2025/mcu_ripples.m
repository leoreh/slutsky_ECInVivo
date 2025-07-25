
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYZE AND LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % run the analysis on multiple sessions
% % -------------------------------------------------------------------------
% basepaths = [mcu_sessions('wt_bsl_ripp'), mcu_sessions('mcu_bsl')];
% nfiles = length(basepaths);
% 
% v = basepaths2vars('basepaths', basepaths, 'vars', {'session'});
% 
% for ifile = 1 : nfiles
%     basepath = basepaths{ifile};
%     cd(basepath)
%     [~, basename] = fileparts(basepath);
% 
%     rippCh = v(ifile).session.channelTags.Ripple;
% 
%     ripp = ripp_wrapper('basepath', pwd, 'rippCh', rippCh,...
%         'limState', 4, 'flgRefine', true, 'flgGraphics', true,...
%         'flgSaveVar', true, 'flgSaveFig', true);
% end


% -------------------------------------------------------------------------
% Re-run only part of the analysis pipeline (e.g., ripp_spks)
% Assumes data was already loaded
for iGrp = 1 : length(grps)
    basepaths = mcu_sessions(grps{iGrp});
    for ifile = 1 : length(basepaths)

        % File
        basepath = basepaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)

        % Get the preloaded ripple structure
        ripp = v{iGrp}(ifile).ripp;

        % Calculate range of filtered lfp
        mapData = ripp.maps.filt;
        idxWin = round(size(mapData, 2) / 2);
        idxWin = [idxWin - 5 : idxWin + 5];
        ripp.peakRng = range(mapData(:, idxWin), 2);
        
        % Calculate peakEnergy using RMS
        ripp.peakEnergy = sqrt(mean(mapData(:, idxWin).^2, 2)).^2;
        ripp = orderfields(ripp);

        % % Re-run ripp_spks
        % ripp = ripp_spks(ripp, 'basepath', basepath, 'limState', 4,...
        %     'flgGraphics', true, 'flgSaveVar', true);

        % Re-run state analysis
        % ripp = ripp_states(ripp, 'basepath', basepath, 'flgGraphics', false, ...
        %     'flgSaveVar', true, 'flgSaveFig', true);

        % plot detection
        % ripp_plot(ripp, 'basepath', basepath, 'flgSaveFig', false);

        % Plot spikes
        % ripp_plotSpks(ripp, 'basepath', basepath, 'flgSaveFig', false);

        % phase coupling
        % ripp = v(ifile).ripp;
        % lfpTimes = ripp.times;
        % rippCh = v(ifile).session.channelTags.Ripple;
        % spkLfp = spklfp_calc('basepath', basepath, 'lfpTimes', lfpTimes,...
        %     'ch', rippCh, 'fRange', [120 200],...
        %     'flgSave', false, 'flgGraphics', true, 'flgStl', false);
        % ripp.spkLfp = spkLfp;

        % Update the preloaded variable structure and save
        v{iGrp}(ifile).ripp = ripp;
        save(fullfile(basepath, [basename, '.ripp.mat']), 'ripp', '-v7.3')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Load data 
grps = {'wt_bsl_ripp'; 'mcu_bsl'};
vars = {'ripp', 'units', 'session'};
for iGrp = 1 : length(grps)
    basepaths = mcu_sessions(grps{iGrp});
    v{iGrp} = basepaths2vars('basepaths', basepaths, 'vars', vars);
end

% -------------------------------------------------------------------------
% Organize lme data table
[clr, lbl] = mcu_clr();
clear varMap
varMap.peakFreq = 'peakFreq';
varMap.peakAmp = 'peakAmp';
varMap.dur = 'dur';
varMap.peakFilt = 'peakFilt';
varMap.peakPower = 'peakPower';
varMap.peakEnergy = 'peakEnergy';
varMap.peakRng = 'peakRng';

for iGrp = 1 : length(grps)
    basepaths = mcu_sessions(grps{iGrp});

    % Prepare tag structures for this mouse
    tagAll.Group = lbl.grp{iGrp};
    tagFiles.Name = get_mname(basepaths);

    % Create table
    tblCell{iGrp} = v2tbl('v', [v{iGrp}(:).ripp], 'varMap', varMap, ...
        'tagFiles', tagFiles, 'tagAll', tagAll, 'idxCol', 1);
end
% Combine all tables, clean and organize
lmeData = vertcat(tblCell{:});
lmeData = rmmissing(lmeData);
lmeData.Group = reordercats(lmeData.Group, lbl.grp);

% Hard code exclusion of weak ripples, validated by manual inspection
idxBad = lmeData.peakAmp <= 8;
lmeData(idxBad, :) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables
lblY = {'Peak Frequency (Hz)', 'Peak Amplitude (µV)', 'Duration (ms)',...
    'Peak Energy (µV²)', 'Rate (SWR/s)', 'Density (%)'};
varRsp = {'peakFreq', 'peakAmp', 'dur', 'peakEnergy', 'Rate'};
idxVar = 4;

% LME
frml = [varRsp{idxVar}, ' ~ Group + (1|Name)'];
lmeCfg.contrasts = 'all';
lmeCfg.distribution = 'Normal';
[lmeStats, lmeMdl] = lme_analyse(lmeData, frml, lmeCfg);

% Plot
hFig = lme_plot(lmeData, lmeMdl, 'lmeStats', lmeStats,...
    'idxRow', [], 'ptype', 'bar', 'axShape', 'tall');

% Graphics
hAx = gca;
[barIdx, barLbl] = lme_sigLines(lmeStats, lmeMdl, lmeData,...
    'idxRow', [1, 2]);
plot_sigLines(hAx, {[1, 2]}, {'*'}, 'flgNS', true);
ylabel(hAx, lblY{idxVar})
xlabel(hAx, '');    title(hAx, '')
plot_axSize('hFig', hFig, 'szOnly', true, 'axShape',...
    'tall', 'axHeight', 300)

% save
fname = lme_frml2char(frml, 'rmRnd', true);
fPath = fullfile(lbl.savepath, 'Ripp');
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', [], 'lmeStats', lmeStats, 'lmeMdl', [], 'fPath', fPath)

hFig = gcf;
hAx = hFig.Children(2);
hBar = hAx.Children;
hBar(2).CData = clr.grp(2, :);
hBar(4).CData = clr.grp(1, :);
% hBar(2).FaceColor = 'flat';
fname = 'normRates~GroupXUnitType';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE SPIKES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FR increase during ripples
frml = 'RippSpks ~ Group * UnitType + (1|Mouse)';
varFld = 'normRates';

% organize for lme
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld, 'vCell', vCell);

% run lme
contrasts = 'all';
contrasts = [1 : 6]; % Example of specific contrasts
[lmeStats, lmeCfg] = lme_analyse(lmeData, lmeCfg, 'contrasts', contrasts);

% plot
hFig = lme_plot(lmeData, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'square'); 

% Update labels
hAx = gca;
ylabel(hAx, 'FR Modulation (a.u.)')
xlabel(hAx, '')
title(hAx, '')
hAx.Legend.Location = 'northeast';

% add significance lines
barIdx = {[3, 4], [1, 2]};
barLbl = {'*', 'NS'};
plot_sigLines(hAx, barIdx, barLbl)
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square', 'axHeight', 300)

% save
fname = lme_frml2char(frml, 'rmRnd', true, 'resNew', 'Ripp', 'sfx', ['_', varFld]);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIKE RIPPLE COUPLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Formula
frml = 'RippSpkLfp ~ Group * UnitType + (1|Mouse)';

% organize lme table for plotting
varFld = 'mrl';
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', varFld, 'vCell', vCell);

% grab only significantly coupled units
% varFld = 'mrl';
% [sigData, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
%     'flgEmg', false, 'varFld', 'pVal', 'vCell', vCell);
% sigIdx = sigData.RippSpkLfp < 0.05;
% lmeData = lmeData(sigIdx, :);

% run lme
contrasts = 'all';
contrasts = [1 : 5];
[lmeStats, lmeCfg] = lme_analyse(lmeData, lmeCfg, 'contrasts', contrasts);

% plot
hFig = lme_plot(lmeData, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'square'); 
hAx = gca;
ylabel(hAx, 'Mean Resultant Length'); xlabel(hAx, ''); title(hAx, '');
hAx.Legend.Location = 'southwest';
hAx.YTick = [0, 0.1, 0.2, 0.3];
hAx.YLim = [0, 0.3];

% add significance lines
barIdx = {[1.5, 3.5]};
barLbl = {'*'};
plot_sigLines(hAx, barIdx, barLbl)

% Assert Size
plot_axSize('hFig', hFig, 'szOnly', true, 'axShape', 'square');

% Save
fname = lme_frml2char(frml, 'rmRnd', true, 'resNew', 'Ripp', 'sfx', ['_', varFld]);
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RIPPLE PETH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a mean normalized SU PETH for RS and FS units from all mice of
% a group and create a figure the superimposes WT and MCU-KO (separately
% for RS and FS units.

normMet = 'zscore';        % 'max', 'ctrl', 'none', 'zscore', 'modulation'

% Stores matrices of unit-normalized PETHs for each group
normGrp = cell(1, length(grps));
unitGrp = cell(1, length(grps));

% Process data for each group
for iGrp = 1:length(grps)
    nMice = length(vCell{iGrp}{1});
    normMap = cell(nMice, 1);
    unitData = cell(nMice, 1);

    for iMouse = 1:nMice
        mData = vCell{iGrp}{1}(iMouse);
        ripp = mData.ripp;
        units = mData.units;

        % Calculate mean PETH per unit
        rippMap = squeeze(mean(ripp.spks.su.rippMap, 2, 'omitnan'));
        ctrlMap = squeeze(mean(ripp.spks.su.ctrlMap, 2, 'omitnan'));

        % Calculate unit FR params
        ctrlAvg = mean(ctrlMap, 2, 'omitnan'); 

        % Use std of per-trial control rates for z-scoring
        ctrlRates = mean(ripp.spks.su.ctrlMap, 3, 'omitnan'); 
        ctrlSD = std(ctrlRates, [], 2, 'omitnan');
        ctrlSD(ctrlSD == 0) = 1; % Avoid division by zero

        rippMax = max(rippMap, [], 2);
        rippMax(rippMax == 0) = 1;

        % normalize
        switch normMet
            case 'max'
                normData = rippMap ./ rippMax;
            case 'ctrl'
                normData = rippMap ./ ctrlAvg; % ctrlAvg might be zero for some units
            case 'zscore'
                normData = (rippMap - ctrlAvg) ./ ctrlSD;
            case 'modulation'
                normData = (rippMap - ctrlAvg) ./ (rippMap + ctrlAvg);
            case 'none'
                normData = rippMap;
        end

        % Collect PETHs
        normMap{iMouse} = normData;

        % get unit data
        nUnits = size(normData, 1);
        unitIdx = nan(nUnits, 1);
        unitIdx(units.clean(1, :)) = 1;
        unitIdx(units.clean(2, :)) = 2;
        unitData{iMouse} = unitIdx;
    end
    normGrp{iGrp} = cell2padmat(normMap, 1);
    unitGrp{iGrp} = cell2padmat(unitData, 1);
end

% Get time bins 
mapDur = ripp.spks.info.mapDur * 1000;
nBinsMap = ripp.spks.info.nBinsMap;
timeBins = linspace(mapDur(1), mapDur(2), nBinsMap);

% Figure Parameters
clr(1, :) = [0.3 0.3 0.3];          % Control
clr(2, :) = [0.784 0.667 0.392];    % MCU-KO 
fntSize = 16;
txtUnit = {'pPYR', 'pINT'};
txtGrp = {'Control', 'MCU-KO'};

% initialize
[hFig, hAx] = plot_axSize('szOnly', false);

% Plot each group
nGrp = length(grps);
iUnit = 1;
hAx = gca; cla; hold on
set(hAx, 'FontName', 'Arial', 'FontSize', fntSize);

% Plot each group in reverse order so Control appears on top
txtLgd = cell(nGrp,1);
clear hPlt
flgRvrs = false;
if flgRvrs
    grpOrdr = nGrp : -1 : 1;
else
    grpOrdr = 1 : nGrp;
end
for iGrp = grpOrdr 
    % Get data for current unit type and group
    unitIdx = unitGrp{iGrp} == iUnit;
    pethData = normGrp{iGrp}(unitIdx, :);

    hPlt(iGrp) = plot_stdShade('hAx', hAx, 'dataMat', pethData,...
        'alpha', 0.3, 'clr', clr(iGrp, :), 'xVal', timeBins);
    hPlt(iGrp).DisplayName = txtGrp{iGrp}; % Assign DisplayName for legend
    txtLgd{iGrp} = txtGrp{iGrp};
end

% Add zero line
xline(hAx, 0, '--k');

% Update labels
ylabel(hAx, 'Norm. FR (z-score)')
xlabel(hAx, 'Time (ms)')
title(hAx, txtUnit{iUnit})

% Set limits
% ylim(hAx, [-0.5, 0.5])
xlim(hAx, [-75 75])
legend(hPlt, txtLgd, 'Location', 'northwest',...
    'FontName', 'Arial', 'FontSize', fntSize);

% Assert Size
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square');

% Save
fname = ['Ripp~FR_', txtUnit{iUnit}];
lme_save('fh', hFig, 'fname', fname, 'frmt', {'fig', 'svg'},...
    'lmeData', table(), 'lmeStats', table(), 'lmeCfg', []) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ripple LFP traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores matrices of ripple LFP traces for each group
nGrp = length(grps);
lfpGrp = cell(1, nGrp);

% Process data for each group
for iGrp = 1:nGrp
    nMice = length(v{iGrp});
    lfpMap = cell(nMice, 1);

    for iMouse = 1:nMice
        mData = v{iGrp}(iMouse);
        ripp = mData.ripp;
        
        % Select ripp based on peakAmp
        idxRipp = ripp.peakAmp > 8 & ripp.peakAmp < 18;

        % Get mean ripple LFP trace
        lfpMap{iMouse} = mean(ripp.maps.raw(idxRipp, :), 1, 'omitnan');
    end
    lfpGrp{iGrp} = cell2padmat(lfpMap, 1); 
end

% Get time bins 
mapDur = ripp.spks.info.mapDur * 1000;
nBinsMap = ripp.spks.info.nBinsMap;
timeBins = linspace(mapDur(1), mapDur(2), nBinsMap);

% initialize
[clr, lbl] = mcu_clr();
[hFig, hAx] = plot_axSize('szOnly', false);

% Plot each group in reverse order so Control appears on top
txtLgd = cell(nGrp,1);
clear hPlt
for iGrp = nGrp : -1 : 1 % Use a different loop variable
    hPlt(iGrp) = plot_stdShade('hAx', hAx, 'dataMat', lfpGrp{iGrp},...
        'alpha', 0.3, 'clr', clr.grp(iGrp, :), 'xVal', timeBins);
    hPlt(iGrp).DisplayName = lbl.grp{iGrp}; % Assign DisplayName for legend
end

% Graphics
xline(hAx, 0, '--k');
ylabel(hAx, 'LFP (µV)')
xlabel(hAx, 'Time (ms)')
xlim(hAx, [-50 50])
ylim(hAx, [-45, 45])
legend(hPlt, lbl.grp, 'Location', 'southwest');
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square');

% Save
fname = 'Ripp~LFP Trace';
lme_save('fh', hFig, 'fname', fname, 'frmt', {'fig', 'svg'},...
    'lmeData', table(), 'lmeStats', table(), 'lmeCfg', []) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT firing rate in ripples vs random
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Re-use frml from RippSpks section or define if not available
frml = 'RippSpks ~ Group * UnitType + (1|Mouse)';

% organize lme table for plotting
[lmeData, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', 'rippRates', 'vCell', vCell);
rippFr = lmeData.RippSpks;
[lmeData, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', 'ctrlRates', 'vCell', vCell);
lmeData = addvars(lmeData, rippFr, 'After', 'UnitID'); 
lmeData.Properties.VariableNames{'RippSpks'} = 'randFr';
lmeData = movevars(lmeData, 'rippFr', 'Before', 'randFr');

% Figure Parameters
clr(1, :) = [0.3 0.3 0.3];          % Control
clr(2, :) = [0.784 0.667 0.392];    % MCU-KO 
fntSize = 16;
txtUnit = {'pPYR', 'pINT'};
txtGrp = {'Control', 'MCU-KO'};

% initialize
[hFig, hAx] = plot_axSize('szOnly', false);

% Plot each group
nGrp = length(vCell);
iUnit = 1;
clear hSct
for iGrp = 1 : nGrp % Use different loop var
    idxTbl = lmeData.UnitType == categorical(txtUnit(iUnit)) &...
        lmeData.Group == categorical(txtGrp(iGrp));
    grpTbl = lmeData(idxTbl, :);
    hSct(iGrp) = scatter(grpTbl.randFr, grpTbl.rippFr, 40, ...
        clr(iGrp, :), 'filled', ...
        'MarkerFaceAlpha', 0.5);
    hSct(iGrp).AlphaData = ones(sum(idxTbl), 1) * 0.9;  % Set a single alpha value for all points
    hSct(iGrp).MarkerFaceAlpha = 'flat';
end

set(hAx, 'XScale', 'log', 'YScale', 'log')
eqLim = [0.005, 100];
xlim(eqLim)
ylim(eqLim)
plot(hAx, eqLim, eqLim, '--k', 'LineWidth', 2, 'HandleVisibility', 'off')
tickVals = get(hAx, 'YTick');
xticks(tickVals)

% Update labels
ylabel(hAx, 'FR in SWR (Hz)')
xlabel(hAx, 'FR in Random (Hz)')
title(hAx, txtUnit{iUnit})
legend(hSct, txtGrp, 'Location', 'southeast');

% Assert Size
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square');

% Save
fname = ['Ripp~FRvsRand_', txtUnit{iUnit}];

lme_save('fh', hFig, 'fname', fname, 'frmt', {'fig', 'svg'},...
    'lmeData', table(), 'lmeStats', table(), 'lmeCfg', []) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIKE RIPPLE COUPLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Formula
frml = 'RippSpkLfp ~ Group * UnitType + (1|Mouse)';

% organize lme table for plotting
[lmeData, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', 'mrl', 'vCell', vCell);
MRL = lmeData.RippSpkLfp;
[lmeData, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', 'pVal', 'vCell', vCell);
pVal = lmeData.RippSpkLfp;
[lmeData, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varFld', 'theta', 'vCell', vCell);
lmeData.Properties.VariableNames{'RippSpkLfp'} = 'Theta';
lmeData = addvars(lmeData, MRL, 'After', 'Theta'); 
lmeData = addvars(lmeData, pVal, 'After', 'MRL'); 

% Figure Parameters
clr(1, :) = [0.3 0.3 0.3];          % Control
clr(2, :) = [0.784 0.667 0.392];    % MCU-KO 
fntSize = 16; FntName = 'Arial';
txtUnit = {'pPYR', 'pINT'};
txtGrp = {'Control', 'MCU-KO'};

% initialize
hFig = figure;

% Plot each group
nGrp = length(vCell);
iUnit = 2;
for iGrp = 1 : nGrp
    % Get specific data from table
    idxUnit = lmeData.UnitType == categorical(txtUnit(iUnit));
    idxGrp = lmeData.Group == categorical(txtGrp(iGrp));
    idxSgn = lmeData.pVal < 0.05;
    idxTbl = idxUnit & idxGrp & idxSgn;
    grpTbl = lmeData(idxTbl, :);

    % Plot units, colored by type if population info is available.
    hPlt = polarscatter(grpTbl.Theta, grpTbl.MRL, 50, ...
                       clr(iGrp, :), 'filled', ...
                       'MarkerFaceAlpha', 0.3);
    hold on;
end
rlim([0 0.6])
rticks(0 : 0.3 : 1)
thetaticks(0:90:270)
hAx = gca;
hAx.ThetaAxisUnits = 'degrees';
hAx.GridAlpha = 0.2;
title(txtUnit{iUnit});
legend(txtGrp, 'Location', 'northwest', 'Interpreter', 'none');
set(hAx, 'FontName', 'Arial', 'FontSize', fntSize);
set(hFig, 'Color', 'w');
hTtl = get(hAx, 'Title');
set(hTtl, 'FontName', FntName, 'FontSize', fntSize + 4, 'FontWeight', 'bold');

% Assert Size
plot_axSize('hFig', hFig, 'szOnly', true, 'axShape', 'square', 'axHeight', 300);

% Save
fname = ['Ripp~SpkPolar_', txtUnit{iUnit}];
lme_save('fh', hFig, 'fname', fname, 'frmt', {'fig', 'svg'},...
    'lmeData', table(), 'lmeStats', table(), 'lmeCfg', []) 




% Plot rate map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSgn = nan(2, 2);
prctSgn = nan(2, 2);

% get map data
for iGrp = 1 : 2
    for iUnit = 1 : 2
        nMice = length(grppaths{iGrp});
        mapData = cell(nMice, 1);
        for iMouse = 1 : nMice
            ripp = vCell{iGrp}{1}(iMouse).ripp;
            mapData{iMouse} = ripp.spkLfp.rateMap.rate;
        end
        mapData = cell2padmat(mapData, 3);
        rateMap = ripp.spkLfp.rateMap;

        % get unit indices
        frml = 'RippSpkLfp ~ Group * UnitType + (1|Mouse)';
        [lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
            'flgEmg', false, 'varFld', 'pVal', 'vCell', vCell);
        idxGrp = lmeData.Group == categorical(txtGrp(iGrp));
        grpTbl = lmeData(idxGrp, :);
        idxUnit = grpTbl.UnitType == categorical(txtUnit(iUnit));
        idxSgn = grpTbl.RippSpkLfp < 0.05;
        idxMap = idxUnit & idxSgn;
        
        % store number of significant units
        nSgn(iGrp, iUnit) = sum(idxUnit & idxSgn);
        prctSgn(iGrp, iUnit) = sum(idxUnit & idxSgn) / sum(idxUnit) * 100;

        % Plot Mean Power-Phase Rate Map (averaged across cells).
        % This 2D heatmap shows the average firing rate of neurons as a function of
        % LFP phase (x-axis) and LFP power (y-axis). The phase axis is duplicated
        % (0 to 4*pi) to visualize cyclic nature. A cosine wave is overlaid as a phase reference.

        % Figure Parameters
        clr(1, :) = [0.3 0.3 0.3];          % Control
        clr(2, :) = [0.784 0.667 0.392];    % MCU-KO
        fntSize = 16;
        txtUnit = {'pPYR', 'pINT'};
        txtGrp = {'Control', 'MCU-KO'};

        % initialize
        [hFig, hAx] = plot_axSize('szOnly', false);

        mapAvg = mean(mapData(:, :, idxMap), 3, 'omitnan'); % Average rate map across units.
        imagesc(hAx, rateMap.phaseBins, rateMap.powBins, mapAvg);
        hold on
        % Overlay cosine wave for phase reference.
        plot(hAx, linspace(0, 2*pi, 100), ...
            cos(linspace(0, 2*pi, 100)) * (range(rateMap.powBins)/4) + mean(rateMap.powBins), ...
            'k--', 'LineWidth', 0.5);
        axis xy
        hCb = colorbar;
        hCb.Label.String = 'Firing Rate (Hz)';
        colormap(hAx, "pink")
        ylim(hAx, [min(rateMap.powBins), max(rateMap.powBins)])
        xlim([0 2 * pi])
        xticks(0:pi/2:2*pi)
        hAx.XTickLabel = {'0', '90', '180', '270', '360'};
        xlabel('Phase (°)')
        ylabel('LFP Power (z-score)');
        title([txtGrp{iGrp}, ': ', txtUnit{iUnit}]);

        % Assert Size
        plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square', 'axHeight', 300);

        % Save
        fname = ['Ripp~SpkPhaseMap_', txtGrp{iGrp}, '_', txtUnit{iUnit}];
        lme_save('fh', hFig, 'fname', fname, 'frmt', {'fig', 'svg'},...
            'lmeData', table(), 'lmeStats', table(), 'lmeCfg', [])

    end
end










