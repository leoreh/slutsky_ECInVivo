
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the analysis on multiple sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths = [mcu_sessions('wt_bsl_ripp'), mcu_sessions('mcu_bsl')];
nfiles = length(basepaths);

v = basepaths2vars('basepaths', basepaths, 'vars', {'session'});

for ifile = 1 : nfiles
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    rippCh = v(ifile).session.channelTags.Ripple;

    ripp = ripp_wrapper('basepath', pwd, 'rippCh', rippCh,...
        'limState', 4, 'flgRefine', true, 'flgGraphics', true,...
        'flgSaveVar', false, 'flgSaveFig', true);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and organize data from all sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define groups and get basepaths
grps = {'wt_bsl_ripp'; 'mcu_bsl'};
clear grppaths
for iGrp = 1 : length(grps)
    grppaths{iGrp} = string(mcu_sessions(grps{iGrp})');
end
basepaths = vertcat(grppaths{:});

% Preload existing ripple data
v = basepaths2vars('basepaths', basepaths, 'vars', {'ripp', 'units'});

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-run only part of the analysis pipeline (e.g., ripp_spks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumes data was already loaded

for ifile = 1 : length(basepaths)

    % File
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    
    oldfile = fullfile(basepath, [basename, '.mat']);
    newfile = fullfile(basepath, [basename, '.ripp.mat']);
    % movefile(oldfile, newfile); % Moving this line as it might cause issues if run multiple times
                                % And check if the old file exists before moving
    if exist(oldfile, 'file')
        movefile(oldfile, newfile);
    end


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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FR increase during ripples
frml = 'RippSpks ~ Group * UnitType + (1|Mouse)';
varField = 'normRates';

% organize for lme
[lmeTbl, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', varField, 'vCell', vCell);

% run lme
contrasts = 'all';
% contrasts = [1 : 5, 8]; % Example of specific contrasts
[lmeResults, lmeCfg] = lme_analyse(lmeTbl, lmeCfg, 'contrasts', contrasts);

% plot
hndFig = lme_plot(lmeTbl, lmeCfg.lmeMdl, 'ptype', 'bar', 'figShape', 'square'); % Changed lmeCfg.mdl to lmeCfg.lmeMdl

% Update labels
axh = gca;
ylabel(axh, 'FR Modulation', 'FontSize', 20)
xlabel(axh, '', 'FontSize', 16)
axh.XAxis.FontSize = 20;
title(axh, '')
axh.Legend.Location = 'northeast';

ylbl = 'Ripp FR';
fname = lme_frml2char(frml, 'rmRnd', true, 'resNew', ylbl); % rm_rnd to rmRnd

% save
lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeTbl', lmeTbl, 'lmeResults', lmeResults, 'lmeCfg', lmeCfg)



% -------------------------------------------------------------------------


% analyze ripple parameter
frml = 'Ripp ~ Group + (1|Mouse)';

% organize for lme
varField = 'rate';

% organize for lme
[lmeTbl, lmeCfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flgEmg', false, 'varField', varField, 'vCell', vCell);

% run lme
contrasts = 'all';
[lmeResults, lmeCfg] = lme_analyse(lmeTbl, lmeCfg, 'contrasts', contrasts);

% plot
hndFig = lme_plot(lmeTbl, lmeCfg.lmeMdl, 'ptype', 'bar', 'figShape', 'tall'); % Changed lmeCfg.mdl to lmeCfg.lmeMdl
axh = gca;
ylabel(axh, 'Rate (SWR/s)', 'FontSize', 20)
xlabel(axh, '', 'FontSize', 20)
title(axh, '')
axh.XAxis.FontSize = 20;
axh.XTickLabelRotation = 0;

ylbl = 'Ripp Rate';
fname = lme_frml2char(frml, 'rmRnd', true, 'resNew', ylbl); % rm_rnd to rmRnd

% save
lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeTbl', lmeTbl, 'lmeResults', lmeResults, 'lmeCfg', lmeCfg)









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concate and plot ripple PETH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a mean normalized SU PETH for RS and FS units from all mice of
% a group and create a figure the superimposes WT and MCU-KO (separately
% for RS and FS units.

normMet = 'zscore';        % 'max', 'ctrl', 'none', 'zscore', 'modulation'

% Figure Parameters
clr = zeros(2,3); % Initialize clr
clr(1, :) = [0.3 0.3 0.3];          % Control
clr(2, :) = [0.784 0.667 0.392];    % MCU-KO 
fntSize = 16;
txtUnit = {'pPYR', 'pPV'};
txtGrp = {'Control', 'MCU-KO'};

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
        
        % Calculate mean PETH per unit
        rippMap = squeeze(mean(ripp.spks.su.rippMap, 2, 'omitnan'));
        ctrlMap = squeeze(mean(ripp.spks.su.ctrlMap, 2, 'omitnan'));
        
        % Calculate unit FR params
        ctrlAvg = mean(ctrlMap, 2, 'omitnan'); % NUnit x 1, average rate over all bins in control PETH for each unit
        
        % Use std of per-trial control rates for z-scoring, as in lme_load/org_rippSpks
        if isfield(mouseData.ripp.spks.su, 'ctrlRates') && ~isempty(mouseData.ripp.spks.su.ctrlRates)
            ctrlSD = std(mouseData.ripp.spks.su.ctrlRates, [], 2, 'omitnan'); % NUnit x 1
        else % Fallback if ctrlRates is not available or empty
            ctrlSD = std(ctrlMap, 0, 2, 'omitnan'); % Std over bins for each unit as a fallback
        end
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
        units = mouseData.units;
        unitIdx = nan(nUnits, 1);
        if isfield(units, 'clean') && ~isempty(units.clean) % Check if units.clean exists and is not empty
            if size(units.clean,1) >=1 && ~isempty(units.clean(1, :))
                 unitIdx(units.clean(1, units.clean(1,:) <= nUnits)) = 1; % Added boundary check
            end
            if size(units.clean,1) >=2 && ~isempty(units.clean(2, :))
                unitIdx(units.clean(2, units.clean(2,:) <= nUnits)) = 2; % Added boundary check
            end
        end
        unitData{iMouse} = unitIdx;

    end
    normGrp{iGrp} = cell2padmat(normMap, 1); 
    unitGrp{iGrp} = cell2padmat(unitData, 1);
end

% Get time bins 
mapDur = ripp.spks.info.mapDur * 1000;
nBinsMap = ripp.spks.info.nBinsMap;
timeBins = linspace(mapDur(1), mapDur(2), nBinsMap);

% initialize
hndFig = figure;
set(hndFig, 'Color', 'w');
fhUnits = get(hndFig, 'Units');
set(hndFig, 'Units', 'pixels');
fhPos = get(hndFig, 'Position');
fhPos(3) = fhPos(4) * 2;
set(hndFig, 'Position', fhPos);
set(hndFig, 'Units', fhUnits);
hndTil = tiledlayout(1, 2);
hndTil.TileSpacing = 'tight';
hndTil.Padding = 'tight';

nGrp = length(grps); % Define nGrp
% Plot for each unit type
for iUnit = 1 : 2
    axh = nexttile(hndTil, iUnit, [1, 1]); cla; hold on
    set(axh, 'FontName', 'Arial', 'FontSize', fntSize);
    
    % Plot each group
    legendEntries = {}; % For legend
    plotHandles = [];
    for iGrpPlot = 1 : nGrp % Use a different loop variable to avoid conflict with outer iGrp if any
        % Get data for current unit type and group
        if ~isempty(unitGrp{iGrpPlot}) % Check if unitGrp data exists for this group
            unitIdx = unitGrp{iGrpPlot} == iUnit;
            if any(unitIdx) % Check if there are any units of this type in this group
                pethData = normGrp{iGrpPlot}(unitIdx, :);
                
                % Plot with std shade
                ph = plot_stdShade('axh', axh, 'dataMat', pethData,...
                    'alpha', 0.3, 'clr', clr(iGrpPlot, :), 'xVal', timeBins);
                legendEntries{end+1} = txtGrp{iGrpPlot};
                plotHandles(end+1) = ph; 
            end
        end
    end
    
    % Add zero line
    xline(axh, 0, '--k');
    
    % Update labels
    ylabel(axh, 'Norm. FR', 'FontSize', 20)
    xlabel(axh, 'Time (ms)', 'FontSize', 20)
    title(axh, txtUnit{iUnit}, 'FontSize', fntSize + 4, 'FontName', 'Arial')
    
    % Set limits
    % ylim(axh, [-0.5, 0.5])
    xlim(axh, mapDur)
    if iUnit == 1 && ~isempty(plotHandles) % Add legend only to the first plot with actual data
        legend(plotHandles, legendEntries, 'Location', 'northeast',...
            'FontName', 'Arial', 'FontSize', fntSize);
    end
end

% Save
fname = 'Ripp PETH';
% lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
%     'lmeTbl', table(), 'lmeResults', table(), 'lmeCfg', [])


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
hndFig = figure;
set(hndFig, 'Color', 'w');
fhUnits = get(hndFig, 'Units');
set(hndFig, 'Units', 'pixels');
fhPos = get(hndFig, 'Position');
fhPos(3) = fhPos(4);
set(hndFig, 'Position', fhPos);
set(hndFig, 'Units', fhUnits);
hndTil = tiledlayout(1, 1);
hndTil.TileSpacing = 'tight';
hndTil.Padding = 'tight';

% Plot
axh = nexttile(hndTil, 1, [1, 1]); cla; hold on
set(axh, 'FontName', 'Arial', 'FontSize', fntSize);

ph = gobjects(nGrp,1); % Preallocate plot handle array
legendTxt = cell(nGrp,1);
validPlots = 0;
% Plot each group in reverse order so Control appears on top
for iGrpPlot = nGrp : -1 : 1 % Use a different loop variable
    if ~isempty(lfpGrp{iGrpPlot})
        plotHandles(iGrpPlot) = plot_stdShade('axh', axh, 'dataMat', lfpGrp{iGrpPlot},...
            'alpha', 0.3, 'clr', clr(iGrpPlot, :), 'xVal', timeBins);
        % plotHandles(iGrpPlot).DisplayName = txtGrp{iGrpPlot}; % Assign DisplayName for legend
        validPlots = validPlots + 1;
        ph(validPlots) = plotHandles(iGrpPlot);
        legendTxt{validPlots} = txtGrp{iGrpPlot};
    end
end

% Add zero line
xline(axh, 0, '--k');

% Update labels
ylabel(axh, 'LFP (µV)', 'FontSize', 20)
xlabel(axh, 'Time (ms)', 'FontSize', 20)
title(axh, '', 'FontSize', fntSize + 4, 'FontName', 'Arial')

% Set limits
xlim(axh, mapDur)
ylim(axh, [-150, 200])

% Add legend with custom order for valid plots
if validPlots > 0
    legend(ph(1:validPlots), legendTxt(1:validPlots), 'Location', 'northeast',...
        'FontName', 'Arial', 'FontSize', fntSize);
end

% Save
fname = 'Ripp LFP';
lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg'},...
    'lmeTbl', table(), 'lmeResults', table(), 'lmeCfg', []) % Removed mat/xlsx for this plot


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
    frmlRippSpks = 'RippSpks ~ Group * UnitType + (1|Mouse)';
else
    frmlRippSpks = frml;
end

% organize lme table for plotting
[lmeTblRipp, ~] = lme_org('grppaths', grppaths, 'frml', frmlRippSpks,...
    'flgEmg', false, 'varField', 'rippRates', 'vCell', vCell);
rippFr = lmeTblRipp.RippSpks;

[lmeTblCtrl, ~] = lme_org('grppaths', grppaths, 'frml', frmlRippSpks,...
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
hndFig = figure;
set(hndFig, 'Color', 'w');
fhUnits = get(hndFig, 'Units');
set(hndFig, 'Units', 'pixels');
fhPos = get(hndFig, 'Position');
fhPos(3) = fhPos(4) * 2;
set(hndFig, 'Position', fhPos);
set(hndFig, 'Units', fhUnits);
hndTil = tiledlayout(1, 2);
hndTil.TileSpacing = 'tight';
hndTil.Padding = 'tight';

% plot
nGrp = length(vCell);
legendEntriesScatter = {};
scatterHandles = [];
for iUnit = 1 : 2
    axh = nexttile(hndTil, iUnit, [1, 1]); cla; hold on
    axis(axh, 'square');
    set(axh, 'FontName', 'Arial', 'FontSize', fntSize);
    
    % Plot each group
    tempHandles = [];
    tempEntries = {};
    for iGrpPlot = 1 : nGrp % Use different loop var
        idxTbl = lmeTblPlot.UnitType == categorical(txtUnit(iUnit)) &...
            lmeTblPlot.Group == categorical(txtGrp(iGrpPlot));
        if any(idxTbl)
            grpTbl = lmeTblPlot(idxTbl, :);
            hndSc = scatter(axh, grpTbl.randFr, grpTbl.rippFr,...
                30, clr(iGrpPlot, :), 'filled');
            hndSc.AlphaData = ones(sum(idxTbl), 1) * 0.9;  % Set a single alpha value for all points
            hndSc.MarkerFaceAlpha = 'flat';
            tempHandles(end+1) = hndSc;
            tempEntries{end+1} = txtGrp{iGrpPlot};
        end
    end
    if iUnit == 1 % Store legend info from first plot only
        scatterHandles = tempHandles;
        legendEntriesScatter = tempEntries;
    end
    
    set(axh, 'XScale', 'log', 'YScale', 'log')
    currentXLim = xlim(axh);
    currentYLim = ylim(axh);
    allMin = min([currentXLim(1), currentYLim(1), 0.01]); % Ensure 0.01 is considered
    allMax = max([currentXLim(2), currentYLim(2), 100]);   % Ensure 100 is considered (or 200 based on original)
    eqLim = [allMin, allMax];
    eqLim = [max(0.01, eqLim(1)), min(200, eqLim(2))]; % Clamp to original 0.01, 200 range if possible

    xlim(eqLim)
    ylim(eqLim)
    plot(axh, eqLim, eqLim, '--k', 'LineWidth', 2)
    
    % Update labels
    ylabel(axh, 'FR in SWR (Hz)', 'FontSize', 20)
    xlabel(axh, 'FR in Random (Hz)', 'FontSize', 20)
    title(axh, txtUnit{iUnit}, 'FontSize', fntSize + 4, 'FontName', 'Arial')
end
if ~isempty(scatterHandles)
    legend(scatterHandles, legendEntriesScatter, 'Location', 'southeast',...
        'FontName', 'Arial', 'FontSize', fntSize);
end

% save
fname = 'Ripp FR vs Rand FR';
lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg'},...
    'lmeTbl', lmeTblPlot, 'lmeResults', table(), 'lmeCfg', []) % Save plot table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bit2uv = 0.195;

for ifile = 1:length(v) % Iterate from 1, check for basename inside
        
    basepath = char(basepaths(ifile)); % Ensure basepath is char for fileparts
    [~, basename] = fileparts(basepath);
    rippFile = fullfile(basepath, [basename, '.ripp.mat']);

    % Skip if this is not the problematic file or if correction already applied
    if contains(basename, 'lh107') % Assuming lh107 was an example of already correct
        if ~isfield(v(ifile).ripp.info, 'bit2uv') || v(ifile).ripp.info.bit2uv ~= 1
             v(ifile).ripp.info.bit2uv = 1; % Mark as 'corrected' or special case
             ripp = v(ifile).ripp;
             save(rippFile, 'ripp', '-v7.3');
        end
        continue
    end

    % Check if correction has already been applied by checking for the bit2uv field and its value
    if isfield(v(ifile).ripp, 'info') && isfield(v(ifile).ripp.info, 'bit2uv') && v(ifile).ripp.info.bit2uv == bit2uv
        % fprintf('Correction already applied to %s\n', basename);
        continue;
    end

    % Apply correction factor
    v(ifile).ripp.peakFilt = v(ifile).ripp.peakFilt * bit2uv;
    v(ifile).ripp.peakAmp  = v(ifile).ripp.peakAmp  * bit2uv;

    if isfield(v(ifile).ripp, 'maps') % Check if maps field exists
        if isfield(v(ifile).ripp.maps, 'raw')
            v(ifile).ripp.maps.raw  = v(ifile).ripp.maps.raw  * bit2uv;
        end
        if isfield(v(ifile).ripp.maps, 'ripp')
             v(ifile).ripp.maps.ripp = v(ifile).ripp.maps.ripp * bit2uv;
        end
        if isfield(v(ifile).ripp.maps, 'amp')
            v(ifile).ripp.maps.amp  = v(ifile).ripp.maps.amp  * bit2uv;
        end
    end

    % Recalculate correlations if the corr field and necessary data exist
    if isfield(v(ifile).ripp, 'corr')
        v(ifile).ripp.corr.AmpFreq = corr(v(ifile).ripp.peakAmp, v(ifile).ripp.peakFreq, 'rows', 'complete');
        v(ifile).ripp.corr.DurAmp  = corr(v(ifile).ripp.dur, v(ifile).ripp.peakAmp, 'rows', 'complete');
    end

    % Add/update informational fields regarding the correction
    if ~isfield(v(ifile).ripp, 'info') || ~isstruct(v(ifile).ripp.info)
        v(ifile).ripp.info = struct(); % Initialize info if it doesn't exist or not a struct
    end
    v(ifile).ripp.info.bit2uv = bit2uv;

    % Save the updated ripp structure
    ripp = v(ifile).ripp;
    save(rippFile, 'ripp', '-v7.3');
    fprintf('Applied bit2uv correction and saved: %s\n', rippFile);
end








