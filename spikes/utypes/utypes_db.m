function utypes_db
% UTYPES_DB Systematically investigates effects of cell classification schemes
%
% This script loads a defined database of sessions, ensures metrics are
% calculated (using field comparison against a reference), and compares
% two different unit classification methods (2 and 3).
%
% Comparisons include:
%   - Firing rates of pINT and pPYR populations across methods.
%   - Unit identity stability (swaps between pPYR and pINT).
%   - Visualizations using plot_utypes.
%
% HISTORY:
%   Dec 2024 - Systematic classification testing

%% ========================================================================
%  DATABASE
%  ========================================================================

clear

% WT
paths{1} = {...
    'D:\Data\lh123\lh123_221219_094508',...
    'D:\Data\lh126\lh126_230111_091208',...
    'D:\Data\lh119\lh119_221114_081305',...
    'E:\Data\lh86\lh86_210227_070000',...       % OP saline
    'E:\Data\lh87\lh87_210523_100607',...       % OP saline
    'E:\Data\lh93\lh93_210819_221608',...       % IT saline
    'E:\Data\lh95\lh95_210824_083300',...       % IT saline
    'E:\Data\lh111\lh111_220823_094417',...
    'E:\Data\lh112\lh112_220831_100435',...     % OP CNO (no effect)
    'E:\Data\lh119\lh119_221114_081305',...     % OP ketamine (no effect)
    'E:\Data\lh120\lh120_221113_090526',...     % OP ketamine (no effect)
    'E:\Data\lh121\lh121_221210_090041',...     % OP ketamine (maybe effect)
    'E:\Data\lh129\lh129_230123_095540',...     % w/ 20 min LT
    'E:\Data\lh130\lh130_230322_084541',...     % OP MK801 (no effect)
    'E:\Data\lh106\lh106_220512_102302',...     % IP ketmamine
    'E:\Data\lh99\lh99_211218_090630',...       % IT Saline (8 tet OE)
    };

% WT (MCU control)
paths{2} = {...
    'D:\Data\lh96\lh96_220120_090157',...
    'D:\Data\lh100\lh100_220413_111004',...
    'D:\Data\lh107\lh107_220518_091200',...
    'D:\Data\lh122\lh122_221223_092656',...
    'D:\Data\lh142\lh142_231005_091832',...
    };

% MCU-KO
paths{3} = {...
    'D:\Data\lh132\lh132_230413_094013',...
    'D:\Data\lh133\lh133_230413_094013',...
    'D:\Data\lh134\lh134_230504_091744',...
    'D:\Data\lh136\lh136_230519_090043',...
    'D:\Data\lh140\lh140_230619_090023',...
    'D:\Data\lh137\lh137_230516_091852',...
    };

% WT (Refaela)
paths{4} = {...
    'E:\Data\Colleagues\RA\raWT5Gq\raWT5Gq_040221_0655',...
    'E:\Data\Colleagues\RA\raWT2Gq\raWT2Gq_200820_bslDay1',...
    'E:\Data\Colleagues\RA\raTgGq\raTgGq_210820',...
    'E:\Data\Colleagues\RA\raTg4\raTg4_210730_155700',...
    };

% MCU-KO (Refaela)
paths{5} = {...
    'E:\Data\Colleagues\RA\raMCU1\raMCU1_080621_0930',...
    'E:\Data\Colleagues\RA\raMCU2\raMCU2_080621_0930',...
    'E:\Data\Colleagues\RA\raMCU3\raMCU3_20211203_084720',...
    'E:\Data\Colleagues\RA\raMCU4\raMCU4_211220_0834',...
    'E:\Data\Colleagues\RA\raMCU5\raMCU5_220322_1906',...
    };

basepaths = [paths{:}];
nPaths = length(basepaths);


%% ========================================================================
%  ANAYLZE (unit features)
%  ========================================================================

flgSave = true;
flgAnalyze = false;

if flgAnalyze

    % load vars
    vars = {'swv_metrics', 'st_metrics', 'fr', 'spikes'};
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);

    % reference path for checking most-recent fields
    refIdx = 17;

    for iPath = 1 : nPaths

        % files
        basepath = basepaths{iPath};
        [~, basename] = fileparts(basepath);
        cd(basepath)

        % firing rate
        if isempty(v(iPath).fr) || ~isempty(setdiff(fieldnames(v(refIdx).fr), fieldnames(v(iPath).fr)))
            fr = calc_fr(v(iPath).spikes.times, 'basepath', basepath,...
                'graphics', false, 'binsize', 60, 'saveVar', flgSave,...
                'smet', 'GK', 'winBL', [0, Inf], 'winCalc', [0, Inf], 'forceA', true);
        end

        % waveform metrices
        if isempty(v(iPath).swv) || ~isempty(setdiff(fieldnames(v(refIdx).swv), fieldnames(v(iPath).swv)))
            swv = spkwv_metrics('basepath', basepath, 'flgSave', flgSave, 'flgForce', true);
        end

        % Spike timing metrics
        if isempty(v(iPath).st) || ~isempty(setdiff(fieldnames(v(refIdx).st), fieldnames(v(iPath).st)))
            st = spktimes_metrics('spktimes', v(iPath).spikes.times, 'sunits', [],...
                'bins', {[0, Inf]}, 'flgForce', true, 'flgSave', flgSave, 'flgAll', false);
        end

    end
end

%% ========================================================================
%  CLASSIFY
%  ========================================================================

% Re-load data
vars = {'swv_metrics', 'st_metrics', 'fr', 'units'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Create table of features for classification
fetTbl = utypes_features('basepaths', basepaths, 'flgPlot', false, 'v', v);

% Classify
unitType = utypes_classify('basepaths', basepaths, 'altClassify', 2,...
    'flgSave', false, 'fetTbl', fetTbl);
fetTbl.unitType = unitType;

%% ========================================================================
%  INSPECT
%  ========================================================================

% Separation by features
hFig = figure;
hAx = gca;
xFet = 'tp';
yFet = 'lidor';
zFet = 'asym';
plot_utypes('basepaths', basepaths, 'fetTbl', fetTbl, 'flgOther', true,...
    'plotType', 'scatter', 'xFet', xFet, 'yFet', yFet, 'zFet', zFet,...
    'unitType', unitType, 'hAx', hAx);


%% ========================================================================
%  VISUALIZATION (MFR vs File)
%  ========================================================================

% Colors
[clr, lbl] = mcu_clr();
clr = clr.unitType; % [RS; FS; Other]

frLim = [-2.5 2];

figure('Name', 'FR vs File', 'Position', [100 100 1600 800]);

% Prepare data: Map 0 (Other) to 3 for indexing
uType = fetTbl.unitType;
uType(uType == 0) = 3;
yVal = fetTbl.mfr;

% Define Groups based on 'paths' structure
grpIdx = cell(length(paths), 1);
lastIdx = 0;
for iUnit = 1:length(paths)
    nFiles = length(paths{iUnit});
    grpIdx{iUnit} = lastIdx + (1:nFiles);
    lastIdx = lastIdx + nFiles;
end
grpNames = {'WT', 'WT (MCU)', 'MCU-KO', 'WT (Ref)', 'MCU-KO (Ref)'};
nGroups = length(paths);

hTile = tiledlayout(3, nGroups, 'TileSpacing', 'tight', 'Padding', 'tight');

% Iterate Groups (Columns)
axSwarm = gobjects(1, nGroups);
axBox = gobjects(1, nGroups);
axHist = gobjects(1, nGroups);

for iGrp = 1:nGroups

    % Identify files in this group
    currPathIdx = grpIdx{iGrp};
    currBase = basepaths(currPathIdx);
    currNames = get_mname(currBase);

    % Filter Table for this Group
    inGroup = ismember(fetTbl.Name, currNames);

    subY = yVal(inGroup);

    % Remove unused categories so swarmchart only shows files in this group
    subX = removecats(fetTbl.Name(inGroup));

    subType = uType(inGroup);

    % --- Top Plot: Swarm (FR vs File) ---
    axSwarm(iGrp) = nexttile(iGrp);
    hold on

    % Plot loop (RS, FS, Other)
    for iUnit = 1:3
        unitIdx = subType == iUnit;
        if ~any(unitIdx)
            continue;
        end

        swarmchart(subX(unitIdx), subY(unitIdx), 10, clr(iUnit, :), 'filled', ...
            'MarkerFaceAlpha', 0.6, 'XJitterWidth', 0.6);
    end

    if iGrp == 1
        ylabel('FR (log_{10} Hz)', 'Interpreter', 'tex');
    else
        ylabel('');
        set(gca, 'YTickLabel', []);
    end
    ylim(frLim);

    title(grpNames{iGrp}, 'Interpreter', 'none');
    xtickangle(45);
    grid on

    % --- Middle Plot: Input for Box Plot ---
    axBox(iGrp) = nexttile(iGrp + nGroups);

    % Prepare data (Rows for plot_boxMean)
    rsData = subY(subType == 1)';
    fsData = subY(subType == 2)';

    if isempty(rsData), rsData = nan; end
    if isempty(fsData), fsData = nan; end

    plot_boxMean({rsData; fsData}, 1, clr([1,2], :), 0.5, 'box', axBox(iGrp), {'RS', 'FS'});

    if iGrp == 1
        ylabel('FR (log_{10} Hz)', 'Interpreter', 'tex');
    else
        ylabel('');
        set(gca, 'YTickLabel', []);
    end
    ylim(frLim);
    grid on

    % --- Bottom Plot: Distribution (FR) ---
    axHist(iGrp) = nexttile(iGrp + 2*nGroups);
    hold on

    % Histogram for RS and FS
    hHist = gobjects(1, 2);
    nCounts = zeros(1, 2);

    for iUnit = 1:2
        unitIdx = subType == iUnit;
        nCounts(iUnit) = sum(unitIdx);

        if ~any(unitIdx)
            continue;
        end

        currData = subY(unitIdx);

        % Histogram (FR on X-axis)
        hHist(iUnit) = histogram(currData, 'BinWidth', 0.2, 'FaceColor', clr(iUnit,:), ...
            'FaceAlpha', 0.4, 'EdgeColor', 'none', 'Normalization', 'probability');

        % Mean Line
        xline(nanmean(currData), '--', 'Color', clr(iUnit,:), 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
    end

    % Legend for this distribution
    legTxt = {sprintf('RS (n=%d)', nCounts(1)), sprintf('FS (n=%d)', nCounts(2))};
    % Filter for existing handles
    hasH = isgraphics(hHist);
    if any(hasH)
        legend(hHist(hasH), legTxt(hasH), 'Location', 'northeast', 'FontSize', 8);
    end

    xlabel('FR (log_{10} Hz)', 'Interpreter', 'tex');
    if iGrp > 1
        ylabel('');
        set(gca, 'YTickLabel', []);
    else
        ylabel('Prob.');
    end

    xlim(frLim);
    grid on
end

% Formatting
linkaxes([axSwarm, axBox], 'y');
linkaxes(axHist, 'xy');
