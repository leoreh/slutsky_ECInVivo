
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



% load and organize data from all sessions
% -------------------------------------------------------------------------
% Define groups and get basepaths
grps = {'wt_bsl_ripp'; 'mcu_bsl'};
clear grppaths
for iGrp = 1 : length(grps)
    grppaths{iGrp} = string(mcu_sessions(grps{iGrp})');
end
basepaths = vertcat(grppaths{:});

% Preload existing ripple data
v = basepaths2vars('basepaths', basepaths, 'vars', {'ripp', 'units', 'session'});

% Organize data struct for lmeOrg
idxStart = 1;
vCell = cell(1, length(grps)); % Initialize vCell
for iGrp = 1 : length(grppaths)
    nPaths = length(grppaths{iGrp});
    idxEnd = idxStart + nPaths - 1;

    % Assign the relevant subset of the preloaded 'v' to vCell
    vCell{iGrp}{1} = v(idxStart : idxEnd);

    % Update start index for the next group
    idxStart = idxEnd + 1;
end


% Re-run only part of the analysis pipeline (e.g., ripp_spks)
% -------------------------------------------------------------------------
% Assumes data was already loaded
for ifile = 1 : length(basepaths)

    % File
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)

    % oldfile = fullfile(basepath, [basename, '.mat']);
    % newfile = fullfile(basepath, [basename, '.ripp.mat']);
    % movefile(oldfile, newfile); % Moving this line as it might cause issues if run multiple times
                                % And check if the old file exists before moving
    % if exist(oldfile, 'file')
    %     movefile(oldfile, newfile);
    % end

    % % Get the preloaded ripple structure
    % ripp = v(ifile).ripp;
    % ripp = rmfield(ripp, 'spks');
    % 
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

    % Update the preloaded variable structure if needed later
    % v(ifile).ripp = ripp; % This line should be active if ripp is modified above

    % phase coupling
    ripp = v(ifile).ripp;
    lfpTimes = ripp.times;
    rippCh = v(ifile).session.channelTags.Ripple;
    spkLfp = spklfp_calc('basepath', basepath, 'lfpTimes', lfpTimes,...
        'ch', rippCh, 'fRange', [120 200],...
        'flgSave', false, 'flgGraphics', true, 'flgStl', false);
    
    ripp.spkLfp = spkLfp;
    save(fullfile(basepath, [basename, '.ripp.mat']), 'ripp', '-v7.3')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Formula
frml = 'Ripp ~ Group + (1|Mouse)';

% Variables to examine
varFields = {'peakFreq', 'dur', 'peakAmp', 'rate', 'density'};
txtY = {'Peak Frequency (Hz)', 'Duration (ms)', 'Peak Amplitude (uV)', 'Rate (SWR/s)', 'Density (%)'};
nFlds = length(varFields);

for iFld = 1 : nFlds
    varField = varFields{iFld};
    
    % organize for lme
    [lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
        'flgEmg', false, 'varField', varField, 'vCell', vCell);

    % run lme
    contrasts = 'all';
    [lmeStats, lmeCfg] = lme_analyse(lmeData, lmeCfg, 'contrasts', contrasts);

    % plot
    hFig = lme_plot(lmeData, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'tall');
    hAx = gca;
    ylabel(hAx, txtY{iFld})
    xlabel(hAx, '')
    title(hAx, '')

    fname = lme_frml2char(frml, 'rmRnd', true, 'sfx', ['_', varFields{iFld}]);

    % save
    lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
        'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeCfg', lmeCfg)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE SPIKES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FR increase during ripples
frml = 'RippSpks ~ Group * UnitType + (1|Mouse)';
varField = 'normRates';

% organize for lme
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', varField, 'vCell', vCell);

% run lme
contrasts = 'all';
% contrasts = [1 : 5, 8]; % Example of specific contrasts
[lmeStats, lmeCfg] = lme_analyse(lmeData, lmeCfg, 'contrasts', contrasts);

% plot
hFig = lme_plot(lmeData, lmeCfg.lmeMdl, 'ptype', 'bar', 'figShape', 'square'); 

% Update labels
hAx = gca;
ylabel(hAx, 'FR Modulation', 'FontSize', 20)
xlabel(hAx, '', 'FontSize', 16)
hAx.XAxis.FontSize = 20;
title(hAx, '')
hAx.Legend.Location = 'northeast';

ylbl = 'Ripp FR';
fname = lme_frml2char(frml, 'rmRnd', true, 'resNew', ylbl); % rm_rnd to rmRnd

% save
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
        mouseData = vCell{iGrp}{1}(iMouse);
        ripp = mouseData.ripp;
        units = mouseData.units;

        % Calculate mean PETH per unit
        rippMap = squeeze(mean(ripp.spks.su.rippMap, 2, 'omitnan'));
        ctrlMap = squeeze(mean(ripp.spks.su.ctrlMap, 2, 'omitnan'));

        % Calculate unit FR params
        ctrlAvg = mean(ctrlMap, 2, 'omitnan'); 

        % Use std of per-trial control rates for z-scoring
        ctrlSD = std(mouseData.ripp.spks.su.ctrlRates, [], 2, 'omitnan');
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
clr = zeros(2,3); % Initialize clr
clr(1, :) = [0.3 0.3 0.3];          % Control
clr(2, :) = [0.784 0.667 0.392];    % MCU-KO 
fntSize = 16;
txtUnit = {'pPYR', 'pPV'};
txtGrp = {'Control', 'MCU-KO'};

% initialize
hFig = figure;
set(hFig, 'Color', 'w');
fhUnits = get(hFig, 'Units');
set(hFig, 'Units', 'pixels');
fhPos = get(hFig, 'Position');
fhPos(3) = fhPos(4);
set(hFig, 'Position', fhPos);
set(hFig, 'Units', fhUnits);

% Plot each group
nGrp = length(grps);
iUnit = 1;
hAx = gca; cla; hold on
set(hAx, 'FontName', 'Arial', 'FontSize', fntSize);
hndPlt = [];
for iGrp = 1 : nGrp 
    % Get data for current unit type and group
    unitIdx = unitGrp{iGrp} == iUnit;
    pethData = normGrp{iGrp}(unitIdx, :);

    % Plot with std shade
    hndPlt(end + 1) = plot_stdShade('axh', hAx, 'dataMat', pethData,...
        'alpha', 0.3, 'clr', clr(iGrp, :), 'xVal', timeBins);
end

% Add zero line
xline(hAx, 0, '--k');

% Update labels
ylabel(hAx, 'Norm. FR', 'FontSize', 20)
xlabel(hAx, 'Time (ms)', 'FontSize', 20)
title(hAx, txtUnit{iUnit}, 'FontSize', fntSize + 4, 'FontName', 'Arial')

% Set limits
% ylim(axh, [-0.5, 0.5])
xlim(hAx, [-50 50])
legend(hndPlt, txtGrp{iGrp}, 'Location', 'northeast',...
    'FontName', 'Arial', 'FontSize', fntSize);

% Save
fname = 'Ripp PETH';
% lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
%     'lmeData', table(), 'lmeStats', table(), 'lmeCfg', [])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ripple LFP traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure Parameters
% clr already defined
fntSize = 16;
% txtGrp already defined

% Stores matrices of ripple LFP traces for each group
lfpGrp = cell(1, length(grps));

% Process data for each group
for iGrp = 1:length(grps)
    nMice = length(vCell{iGrp}{1});
    lfpMap = cell(nMice, 1);

    for iMouse = 1:nMice
        mouseData = vCell{iGrp}{1}(iMouse);
        ripp = mouseData.ripp;
        
        % Get mean ripple LFP trace
        if isfield(ripp, 'maps') && isfield(ripp.maps, 'raw') && ~isempty(ripp.maps.raw)
            lfpMap{iMouse} = mean(ripp.maps.raw, 1, 'omitnan');
        else
            lfpMap{iMouse} = []; % Handle case where data might be missing
        end
    end
    lfpGrp{iGrp} = cell2padmat(lfpMap, 1); 
end

% Get time bins 
mapDur = ripp.spks.info.mapDur * 1000;
nBinsMap = ripp.spks.info.nBinsMap;
timeBins = linspace(mapDur(1), mapDur(2), nBinsMap);

% initialize
hFig = figure;
set(hFig, 'Color', 'w');
fhUnits = get(hFig, 'Units');
set(hFig, 'Units', 'pixels');
fhPos = get(hFig, 'Position');
fhPos(3) = fhPos(4);
set(hFig, 'Position', fhPos);
set(hFig, 'Units', fhUnits);
hndTil = tiledlayout(1, 1);
hndTil.TileSpacing = 'tight';
hndTil.Padding = 'tight';

% Plot
hAx = nexttile(hndTil, 1, [1, 1]); cla; hold on
set(hAx, 'FontName', 'Arial', 'FontSize', fntSize);

ph = gobjects(nGrp,1); % Preallocate plot handle array
legendTxt = cell(nGrp,1);
validPlots = 0;
% Plot each group in reverse order so Control appears on top
for iGrp = nGrp : -1 : 1 % Use a different loop variable
    if ~isempty(lfpGrp{iGrp})
        hndPlt(iGrp) = plot_stdShade('axh', hAx, 'dataMat', lfpGrp{iGrp},...
            'alpha', 0.3, 'clr', clr(iGrp, :), 'xVal', timeBins);
        % plotHandles(iGrpPlot).DisplayName = txtGrp{iGrpPlot}; % Assign DisplayName for legend
        validPlots = validPlots + 1;
        ph(validPlots) = hndPlt(iGrp);
        legendTxt{validPlots} = txtGrp{iGrp};
    end
end

% Add zero line
xline(hAx, 0, '--k');

% Update labels
ylabel(hAx, 'LFP (µV)', 'FontSize', 20)
xlabel(hAx, 'Time (ms)', 'FontSize', 20)
title(hAx, '', 'FontSize', fntSize + 4, 'FontName', 'Arial')

% Set limits
xlim(hAx, mapDur)
ylim(hAx, [-150, 200])

% Add legend with custom order for valid plots
if validPlots > 0
    legend(ph(1:validPlots), legendTxt(1:validPlots), 'Location', 'northeast',...
        'FontName', 'Arial', 'FontSize', fntSize);
end

% Save
fname = 'Ripp LFP';
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg'},...
    'lmeData', table(), 'lmeStats', table(), 'lmeCfg', []) % Removed mat/xlsx for this plot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT firing rate in ripples vs random
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure Parameters
% clr already defined
fntSize = 16;
% txtUnit already defined
% txtGrp already defined

% Re-use frml from RippSpks section or define if not available
if ~exist('frml','var') || ~strcmp(frml, 'RippSpks ~ Group * UnitType + (1|Mouse)')
    frml = 'RippSpks ~ Group * UnitType + (1|Mouse)';
else
    frml = frml;
end

% organize lme table for plotting
[lmeData, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', 'rippRates', 'vCell', vCell);
rippFr = lmeData.RippSpks;

[lmeTblCtrl, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', 'ctrlRates', 'vCell', vCell);
lmeTblCtrl = addvars(lmeTblCtrl, rippFr, 'After', 'UnitID'); % Ensure consistent column order for merging/plotting
% Make sure UnitID, Mouse, Group, UnitType are present for proper identification, then rename
if ismember('RippSpks', lmeTblCtrl.Properties.VariableNames)
    lmeTblCtrl.Properties.VariableNames{'RippSpks'} = 'randFr';
end
lmeTblPlot = lmeTblCtrl;
% lmeTblPlot = movevars(lmeTblPlot, "rippFr", 'Before', "randFr"); % This was causing issues if randFr was not yet the name
% Ensure variable names are correct before moving
if ismember('rippFr', lmeTblPlot.Properties.VariableNames) && ismember('randFr', lmeTblPlot.Properties.VariableNames)
    lmeTblPlot = movevars(lmeTblPlot, 'rippFr', 'Before', 'randFr');
else
    warning('Could not move rippFr column, check variable names in lmeTblPlot');
end

% initialize
hFig = figure;
set(hFig, 'Color', 'w');
fhUnits = get(hFig, 'Units');
set(hFig, 'Units', 'pixels');
fhPos = get(hFig, 'Position');
fhPos(3) = fhPos(4) * 2;
set(hFig, 'Position', fhPos);
set(hFig, 'Units', fhUnits);
hndTil = tiledlayout(1, 2);
hndTil.TileSpacing = 'tight';
hndTil.Padding = 'tight';

% plot
nGrp = length(vCell);
legendEntriesScatter = {};
scatterHandles = [];
for iUnit = 1 : 2
    hAx = nexttile(hndTil, iUnit, [1, 1]); cla; hold on
    axis(hAx, 'square');
    set(hAx, 'FontName', 'Arial', 'FontSize', fntSize);
    
    % Plot each group
    tempHandles = [];
    tempEntries = {};
    for iGrp = 1 : nGrp % Use different loop var
        idxTbl = lmeTblPlot.UnitType == categorical(txtUnit(iUnit)) &...
            lmeTblPlot.Group == categorical(txtGrp(iGrp));
        if any(idxTbl)
            grpTbl = lmeTblPlot(idxTbl, :);
            scatterMarkerSize = 400; 
            hndSc = polarscatter(grpTbl.randFr, grpTbl.rippFr, scatterMarkerSize, ...
                               clr(iGrp, :), 'filled', ...
                               'MarkerFaceAlpha', 0.5);
            hndSc.AlphaData = ones(sum(idxTbl), 1) * 0.9;  % Set a single alpha value for all points
            hndSc.MarkerFaceAlpha = 'flat';
            tempHandles(end+1) = hndSc;
            tempEntries{end+1} = txtGrp{iGrp};
        end
    end
    if iUnit == 1 % Store legend info from first plot only
        scatterHandles = tempHandles;
        legendEntriesScatter = tempEntries;
    end
    
    set(hAx, 'XScale', 'log', 'YScale', 'log')
    currentXLim = xlim(hAx);
    currentYLim = ylim(hAx);
    allMin = min([currentXLim(1), currentYLim(1), 0.01]); % Ensure 0.01 is considered
    allMax = max([currentXLim(2), currentYLim(2), 100]);   % Ensure 100 is considered (or 200 based on original)
    eqLim = [allMin, allMax];
    eqLim = [max(0.01, eqLim(1)), min(200, eqLim(2))]; % Clamp to original 0.01, 200 range if possible

    xlim(eqLim)
    ylim(eqLim)
    plot(hAx, eqLim, eqLim, '--k', 'LineWidth', 2)
    
    % Update labels
    ylabel(hAx, 'FR in SWR (Hz)', 'FontSize', 20)
    xlabel(hAx, 'FR in Random (Hz)', 'FontSize', 20)
    title(hAx, txtUnit{iUnit}, 'FontSize', fntSize + 4, 'FontName', 'Arial')
end
if ~isempty(scatterHandles)
    legend(scatterHandles, legendEntriesScatter, 'Location', 'southeast',...
        'FontName', 'Arial', 'FontSize', fntSize);
end

% save
fname = 'Ripp FR vs Rand FR';
lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg'},...
    'lmeData', lmeTblPlot, 'lmeStats', table(), 'lmeCfg', []) % Save plot table

spklfp_plot(ripp.spkLfp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIKE RIPPLE COUPLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure Parameters
clr(1, :) = [0.3 0.3 0.3];          % Control
clr(2, :) = [0.784 0.667 0.392];    % MCU-KO 
fntSize = 16;
txtUnit = {'pPYR', 'pPV'};
txtGrp = {'Control', 'MCU-KO'};

% Formula
frml = 'RippSpkLfp ~ Group * UnitType + (1|Mouse)';

% organize lme table for plotting
[lmeData, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', 'mrl', 'vCell', vCell);
MRL = lmeData.RippSpkLfp;
[lmeData, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', 'pVal', 'vCell', vCell);
pVal = lmeData.RippSpkLfp;
[lmeData, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', 'theta', 'vCell', vCell);
lmeData.Properties.VariableNames{'RippSpkLfp'} = 'Theta';
lmeData = addvars(lmeData, MRL, 'After', 'Theta'); 
lmeData = addvars(lmeData, pVal, 'After', 'MRL'); 

% Initialize figure
hFig = figure;
set(hFig, 'Color', 'w');
fhUnits = get(hFig, 'Units');
set(hFig, 'Units', 'pixels');
fhPos = get(hFig, 'Position');
fhPos(3) = fhPos(4);
set(hFig, 'Position', fhPos);
set(hFig, 'Units', fhUnits);

% Plot each group
nGrp = length(vCell);
iUnit = 1;
for iGrp = 1 : nGrp
    % Get specific data from table
    idxUnit = lmeData.UnitType == categorical(txtUnit(iUnit));
    idxGrp = lmeData.Group == categorical(txtGrp(iGrp));
    idxSgn = lmeData.pVal < 0.05;
    idxTbl = idxUnit & idxGrp & idxSgn;
    grpTbl = lmeData(idxTbl, :);

    % Plot units, colored by type if population info is available.
    hndPlt = polarscatter(grpTbl.Theta, grpTbl.MRL, 50, ...
                       clr(iGrp, :), 'filled', ...
                       'MarkerFaceAlpha', 0.3);
    hold on;
end
rlim([0 0.6])
rticks(0 : 0.3 : 1)
thetaticks(0:90:270)
hndAx = gca;
hndAx.ThetaAxisUnits = 'degrees';
hndAx.GridAlpha = 0.2;
title(txtUnit{iUnit});
legend(txtGrp, 'Location', 'northwest', 'Interpreter', 'none');
set(hndAx, 'FontName', 'Arial', 'FontSize', fntSize);

% save
% fname = 'Ripp FR vs Rand FR';
% lme_save('fh', hFig, 'fname', fname, 'frmt', {'svg'},...
%     'lmeData', lmeTblPlot, 'lmeStats', table(), 'lmeCfg', []) % Save plot table


% LME analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Formula
frml = 'RippSpkLfp ~ Group * UnitType + (1|Mouse)';

% organize lme table for plotting
[lmeData, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', 'mrl', 'vCell', vCell);

% run lme
contrasts = 'all';
[lmeStats, lmeCfg] = lme_analyse(lmeData, lmeCfg, 'contrasts', contrasts);

% plot
hFig = lme_plot(lmeData, lmeCfg.lmeMdl, 'ptype', 'bar', 'axShape', 'square'); 
hAx = gca;
ylabel(hAx, 'Mean Resultant Length'); xlabel(hAx, ''); title(hAx, '');
hAx.Legend.Location = 'southeast';
plot_axSize('hFig', hFig, 'szOnly', false)

% add significance lines
barIdx = {[1, 2], [1.5, 3.5], [3, 4]};
barLbl = {'NS', '*', '***'};
plot_sigLines(hAx, barIdx, barLbl)

ylbl = 'Ripp Theta';
fname = lme_frml2char(frml, 'rmRnd', true, 'resNew', ylbl); 

% Plot rate map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get map data
iGrp = 1;
iUnit = 1;
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
    'flgEmg', false, 'varField', 'pVal', 'vCell', vCell);
idxGrp = lmeData.Group == categorical(txtGrp(iGrp));
grpTbl = lmeData(idxGrp, :);
idxUnit = grpTbl.UnitType == categorical(txtUnit(iUnit));
idxSgn = grpTbl.RippSpkLfp < 0.05;
idxMap = idxUnit & idxSgn;

% Plot Mean Power-Phase Rate Map (averaged across cells).
% This 2D heatmap shows the average firing rate of neurons as a function of
% LFP phase (x-axis) and LFP power (y-axis). The phase axis is duplicated
% (0 to 4*pi) to visualize cyclic nature. A cosine wave is overlaid as a phase reference.

% Initialize figure
hFig = figure;
set(hFig, 'Color', 'w');
fhUnits = get(hFig, 'Units');
set(hFig, 'Units', 'pixels');
fhPos = get(hFig, 'Position');
fhPos(3) = fhPos(4);
set(hFig, 'Position', fhPos);
set(hFig, 'Units', fhUnits);
hndAx = gca;

mapAvg = mean(mapData(:, :, idxMap), 3, 'omitnan'); % Average rate map across units.
imagesc(hndAx, rateMap.phaseBins, rateMap.powBins, mapAvg); 
hold on
% Overlay cosine wave for phase reference.
plot(hndAx, linspace(0, 2*pi, 100), ...
     cos(linspace(0, 2*pi, 100)) * (range(rateMap.powBins)/4) + mean(rateMap.powBins), ...
     'k--', 'LineWidth', 0.5);
axis xy
colorbar
xlim([0 2 * pi])
xticks(0:pi/2:2*pi)
hndAx.XTickLabel = {'0', '90', '180', '270', '360'};
xlabel('Phase [o]')
ylabel('Norm. LFP Power (Z-score)');
title(txtUnit{iUnit});
