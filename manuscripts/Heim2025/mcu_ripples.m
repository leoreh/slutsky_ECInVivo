
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
    movefile(oldfile, newfile);


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
    v(ifile).ripp = ripp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FR increase during ripples
frml = 'RippSpks ~ Group * UnitType + (1|Mouse)';
var_field = 'normRates';

% organize for lme
[lme_tbl, lme_cfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flg_emg', false, 'var_field', var_field, 'vCell', vCell);

% run lme
contrasts = 'all';
contrasts = [1 : 5, 8];
[lme_results, lme_cfg] = lme_analyse(lme_tbl, lme_cfg, 'contrasts', contrasts);

% plot
hndFig = lme_plot(lme_tbl, lme_cfg.mdl, 'ptype', 'bar', 'figShape', 'square');

% Update labels
axh = gca;
ylabel(axh, 'FR Modulation', 'FontSize', 20)
xlabel(axh, '', 'FontSize', 16)
axh.XAxis.FontSize = 20;
title(axh, '')
axh.Legend.Location = 'northeast';

ylbl = 'Ripp FR';
fname = lme_frml2char(frml, 'rm_rnd', true, 'resNew', ylbl);

% save
lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lme_tbl', lme_tbl, 'lme_results', lme_results, 'lme_cfg', lme_cfg)



% -------------------------------------------------------------------------


% analyze ripple parameter
frml = 'Ripp ~ Group + (1|Mouse)';

% organize for lme
var_field = 'rate';

% organize for lme
[lme_tbl, lme_cfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flg_emg', false, 'var_field', var_field, 'vCell', vCell);

% run lme
contrasts = 'all';
[lme_results, lme_cfg] = lme_analyse(lme_tbl, lme_cfg, 'contrasts', contrasts);

% plot
hndFig = lme_plot(lme_tbl, lme_cfg.mdl, 'ptype', 'bar', 'figShape', 'tall');
axh = gca;
ylabel(axh, 'Rate (SWR/s)', 'FontSize', 20)
xlabel(axh, '', 'FontSize', 20)
title(axh, '')
axh.XAxis.FontSize = 20;
axh.XTickLabelRotation = 0;

ylbl = 'Ripp Rate';
fname = lme_frml2char(frml, 'rm_rnd', true, 'resNew', ylbl);

% save
lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lme_tbl', lme_tbl, 'lme_results', lme_results, 'lme_cfg', lme_cfg)









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concate and plot ripple PETH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a mean normalized SU PETH for RS and FS units from all mice of
% a group and create a figure the superimposes WT and MCU-KO (separately
% for RS and FS units.

normMet = 'zscore';        % 'max', 'ctrl', 'none', 'zscore', 'modulation'

% Figure Parameters
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
        idxRS = mouseData.units.clean(1, :);
        idxFS = mouseData.units.clean(2, :);

        % Calculate unit FR params
        ctrlAvg = mean(ctrlMap, 2);

        ctrlSe = std(frUnit, [], 2);
        ctrlSe(ctrlSe == 0) = 1;
        rippMax = max(rippMap, [], 2);
        rippMax(rippMax == 0) = 1;

        % normalize
        switch normMet
            case 'max'
                normData = rippMap ./ rippMax;
            case 'ctrl'
                normData = rippMap ./ ctrlAvg;
            case 'zscore'
                normData = (rippMap - ctrlAvg) ./ ctrlSe;
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

% Plot for each unit type
for iUnit = 1 : 2
    axh = nexttile(hndTil, iUnit, [1, 1]); cla; hold on
    set(axh, 'FontName', 'Arial', 'FontSize', fntSize);
    
    % Plot each group
    for iGrp = 1 : nGrp
        % Get data for current unit type and group
        unitIdx = unitGrp{iGrp} == iUnit;
        pethData = normGrp{iGrp}(unitIdx, :);
        
        % Plot with std shade
        plot_stdShade('axh', axh, 'dataMat', pethData,...
            'alpha', 0.3, 'clr', clr(iGrp, :), 'xVal', timeBins);
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
end
% Add legend
legend(txtGrp, 'Location', 'northeast',...
    'FontName', 'Arial', 'FontSize', fntSize);

% Save
fname = 'Ripp PETH';
% lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
%     'lme_tbl', table(), 'lme_results', table(), 'lme_cfg', [])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ripple LFP traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure Parameters
clr(1, :) = [0.3 0.3 0.3];          % Control
clr(2, :) = [0.784 0.667 0.392];    % MCU-KO 
fntSize = 16;
txtGrp = {'Control', 'MCU-KO'};

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
        lfpMap{iMouse} = mean(ripp.maps.raw, 1, 'omitnan');
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

% Plot each group in reverse order so Control appears on top
for iGrp = nGrp : -1 : 1
    % Plot with std shade
    ph(iGrp) = plot_stdShade('axh', axh, 'dataMat', lfpGrp{iGrp},...
        'alpha', 0.3, 'clr', clr(iGrp, :), 'xVal', timeBins);
    ph(iGrp).DisplayName = txtGrp{igrp};
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

% Add legend with custom order
legend(ph, txtGrp(2 : -1 : 1), 'Location', 'northeast',...
    'FontName', 'Arial', 'FontSize', fntSize);

% Save
fname = 'Ripp LFP';
lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lme_tbl', table(), 'lme_results', table(), 'lme_cfg', [])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT firing rate in ripples vs random
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure Parameters
clr(1, :) = [0.3 0.3 0.3];          
clr(2, :) = [0.784 0.667 0.392];    
fntSize = 16;
txtUnit = {'pPYR', 'pPV'};
txtGrp = {'Control', 'MCU-KO'};

% organize lme table for plotting
[lme_tbl, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flg_emg', false, 'var_field', 'rippRates', 'vCell', vCell);
rippFr = lme_tbl.RippSpks;

[lme_tbl, ~] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flg_emg', false, 'var_field', 'ctrlRates', 'vCell', vCell);
lme_tbl = addvars(lme_tbl, rippFr);
lme_tbl.Properties.VariableNames{1} = 'randFr';
lme_tbl = movevars(lme_tbl, "rippFr", 'Before', "randFr");

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
for iUnit = 1 : 2
    axh = nexttile(hndTil, iUnit, [1, 1]); cla; hold on
    axis(axh, 'square');
    set(axh, 'FontName', 'Arial', 'FontSize', fntSize);
    
    % Plot each group
    for iGrp = 1 : nGrp
        idxTbl = lme_tbl.UnitType == categorical(txtUnit(iUnit)) &...
            lme_tbl.Group == categorical(txtGrp(iGrp));
        grp_tbl = lme_tbl(idxTbl, :);
        hndSc = scatter(axh, grp_tbl.randFr, grp_tbl.rippFr,...
            30, clr(iGrp, :), 'filled');
        hndSc.AlphaData = ones(sum(idxTbl), 1) * 0.9;  % Set a single alpha value for all points
        hndSc.MarkerFaceAlpha = 'flat';
    end
    
    set(axh, 'XScale', 'log', 'YScale', 'log')
    eqLim = [min([xlim, ylim]), max([xlim, ylim])];
    eqLim = [0.01, 200];
    xlim(eqLim)
    ylim(eqLim)
    plot(axh, eqLim, eqLim, '--k', 'LineWidth', 2)
    
    % Update labels
    ylabel(axh, 'FR in SWR (Hz)', 'FontSize', 20)
    xlabel(axh, 'FR in Random (Hz)', 'FontSize', 20)
    title(axh, txtUnit{iUnit}, 'FontSize', fntSize + 4, 'FontName', 'Arial')
end
legend(txtGrp, 'Location', 'southeast',...
    'FontName', 'Arial', 'FontSize', fntSize);

% save
fname = 'Ripp FR vs Rand FR';
lme_save('fh', hndFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lme_tbl', lme_tbl, 'lme_results', table(), 'lme_cfg', [])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bit2uv = 0.195;

for ifile = 2:length(v)
        
    basepath = char(basepaths(ifile)); % Ensure basepath is char for fileparts
    [~, basename] = fileparts(basepath);
    rippFile = fullfile(basepath, [basename, '.ripp.mat']);

    if contains(basename, 'lh107')
        v(ifile).ripp.info.bit2uv = 1;
        ripp = v(ifile).ripp;
        save(rippFile, 'ripp', '-v7.3');
        continue
    end

    % Apply correction factor
    v(ifile).ripp.peakFilt = v(ifile).ripp.peakFilt * bit2uv;
    v(ifile).ripp.peakAmp  = v(ifile).ripp.peakAmp  * bit2uv;

    v(ifile).ripp.maps.raw  = v(ifile).ripp.maps.raw  * bit2uv;
    v(ifile).ripp.maps.ripp = v(ifile).ripp.maps.ripp * bit2uv;
    v(ifile).ripp.maps.amp  = v(ifile).ripp.maps.amp  * bit2uv;

    % Recalculate correlations if the corr field and necessary data exist
    % corr function with 'rows', 'complete' handles NaNs and should be okay with empty inputs (returning NaN)
    v(ifile).ripp.corr.AmpFreq = corr(v(ifile).ripp.peakAmp, v(ifile).ripp.peakFreq, 'rows', 'complete');
    v(ifile).ripp.corr.DurAmp  = corr(v(ifile).ripp.dur, v(ifile).ripp.peakAmp, 'rows', 'complete');

    % Add/update informational fields regarding the correction
    v(ifile).ripp.info.bit2uv = bit2uv;

    % Save the updated ripp structure
    ripp = v(ifile).ripp;
    save(rippFile, 'ripp', '-v7.3');
end








