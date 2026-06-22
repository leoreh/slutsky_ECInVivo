
%  MCU_SPIKES - Analyze in vivo MCU spike data

%% ========================================================================
%  LOAD_DATA
%  ========================================================================

basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
cfg = mcu_cfg();

% Load table
presets = {'burst'};
tbl = mcu_tblVivo('basepaths', basepaths, 'flgClean', true, ...
    'presets', presets);

% Assert no zero values
tblLme = tbl_trans(tbl, 'flg0', true, 'verbose', true);

% Limit to RS units
uIdx = tblLme.unitType == 'RS';
tblLme = tblLme(uIdx, :);

% logit pBurst
tblTrans = tbl_trans(tblLme, 'varsInc', {'pBurst'}, 'logBase', 'logit');
tblLme.pBurst_trans = tblTrans.pBurst;

% Indices
idxGrp = tblLme.genotype == 'Control';
idxGrp = tblLme.genotype == 'MCU-KO';

% Plots
tblGUI_scatHist(tblLme(:, :), 'xVar', 'fr', 'yVar', 'pBurst', 'grpVar', 'genotype')
tblGUI_bar(tblLme);


%% ========================================================================
%  BACLOFEN (y ~ Group * Day + (Day|Name))
%  ========================================================================

% Remove WASH
tblLme(tblLme.day == 'WASH', :) = [];
tblLme.day = removecats(tblLme.day, {'WASH'});

tblLme(tblLme.day == 'BAC1', :) = [];
tblLme.day = removecats(tblLme.day, {'BAC1'});

tblLme(tblLme.day == 'BAC2', :) = [];
tblLme.day = removecats(tblLme.day, {'BAC2'});


% Select Params
varRsp = 'pBurst';

% run lme
frml = [varRsp, ' ~ genotype * day + (day|sbjID)'];
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'logit-normal');

% Plot
hFig = tblGUI_bar(tblLme, 'yVar', varRsp, 'xVar', 'day', 'grpVar', 'genotype');


% Prism
iGrp = 1;
[prismMat] = tbl2prism(tblLme(tblLme.genotype == cfg.lbl.grp{iGrp}, :), ...
    'yVar', varRsp, 'grpVar', 'day');



%% ========================================================================
%  BASELINE (y ~ Group * unitType + (1|Name))
%  ========================================================================


% Select data
tblLme = tbl(tbl.day == 'BSL', :);
tblLme = tblLme(tblLme.unitType == 'RS', :);

% Select Params
varRsp = 'fr';

% Fit
frml = [varRsp, ' ~ genotype + (1|sbjID)'];
[lmeMdl, lmeStats, lmeInfo, tblMdl] = lme_analyse(tblLme, frml, 'dist', 'gamma');


% Plot
hFig = tblGUI_bar(tblLme, 'yVar', varRsp, 'xVar', 'genotype');

% Prism
[prismMat] = tbl2prism(tblLme, 'yVar', varRsp, 'grpVar', 'genotype');

log10(mean(prismMat, 'omitnan'))



%% ========================================================================
%  FR VS TIME
%  ========================================================================

% Files
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];

% Grab FR vs Time data
[tAxis, tbl] = mcu_frTbl(basepaths, 'flgPlot', false);

% Define Baseline Window (Indices where Time < 0)
winNorm = [24 * 60, find(tAxis >= 0, 1) - 1];

% Calculate Floor Value (1 event per max time)
% Prevents clipping of valid low-rate dynamics when using Geometric Mean
floorVal = 1 / (diff(winNorm) * 60);

% Filter Units
idxUnits = tbl.unitType == 'RS';
tblPlot = tbl(idxUnits, :);

% Denoise
tblPlot.FRt = fr_denoise(tblPlot.FRt, 'flgPlot', false, 'frameLen', 300);

% Normalize
% We group by 'sbjID' (Mouse) so that all units from the same mouse are
% normalized by the pooled baseline statistics of that mouse.
tblPlot = tbl_tNorm(tblPlot, 'varsInc', 'FRt', 'varsGrp', 'sbjID', ...
    'winNorm', winNorm, 'Method', 'percentage', ...
    'flgGeom', false, 'floorVal', floorVal);

% Plot (Log Scale)
hFig = tblGUI_xy(tAxis, tblPlot, 'yVar', 'FRt', 'tileVar', 'genotype', 'grpVar', 'sbjID');

% Grab to prism
idxUnits = tblPlot.unitType == 'RS' & tblPlot.genotype == 'Control';
frMat = tblPlot.FRt_bins(idxUnits, :)';

prismMat = [mean(frMat, 2, 'omitnan'), ...
    std(frMat, [], 2, 'omitnan'), ...
    sum(~isnan(frMat), 2, 'omitnan')];

miceID = unique(tblPlot.sbjID);
prismMat = [];
for iSbj = 1 : length(miceID)
    idxUnits = tblPlot.unitType == 'RS';
    idxUnits = idxUnits & tblPlot.sbjID == miceID(iSbj);
    frMat = tblPlot.FRt(idxUnits, :)';
    prismMat = [prismMat, mean(frMat, 2, 'omitnan')];
end

% Table per day
statType = 'mean';
tblDay = groupsummary(tblPlot, {'genotype', 'day'}, statType, ...
    vartype("numeric"));

% -------------------------------------------------------------------------
% Binned Bar Plot (mean ± SEM per 6-hr bin, grouped by genotype)
% -------------------------------------------------------------------------

% Bin parameters
binSizeH  = 6;                                    % Bin width [hr]
tBinEdges = -33 : binSizeH : 96;                  % Bin edges [hr]
tBinCents = tBinEdges(1:end-1) + binSizeH / 2;   % Bin centers [hr]
nBins     = length(tBinCents);

% Assign each tAxis point to a bin; points outside [-24, 72] yield NaN
binIdx = discretize(tAxis, tBinEdges);

% Average FRt within each temporal bin per unit → new tblPlot column
FRt_bins = nan(height(tblPlot), nBins);
for iBin = 1 : nBins
    idxT = binIdx == iBin;
    FRt_bins(:, iBin) = mean(tblPlot.FRt(:, idxT), 2, 'omitnan');
end
tblPlot.FRt_bins = FRt_bins;

hFig = tblGUI_xy(tBinCents, tblPlot, 'yVar', 'FRt', 'tileVar', 'genotype');


miceID = unique(tblPlot.sbjID);
prismMat = [];
for iSbj = 1 : length(miceID)
    idxUnits = tblPlot.unitType == 'RS';
    idxUnits = idxUnits & tblPlot.sbjID == miceID(iSbj);
    frMat = tblPlot.FRt_bins(idxUnits, :)';
    prismMat = [prismMat, mean(frMat, 2, 'omitnan')];
end


%% ========================================================================
%  REPRESENTATIVE RASTER
%  ========================================================================

% Files
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
basepaths = natsort(basepaths);
idxFiles = [1, 5, 50, 54];
% idxFiles = [15, 19, 50, 54];

% Load
vars = {'spikes', 'units', 'burst'};
v = basepaths2vars('basepaths', basepaths(idxFiles), 'vars', vars);

% Config
cfg = mcu_cfg();
winPlot = [0, 60];
lnW = 1;
lnH = 1;

% Exact inner axis size (cm) and surrounding white-space
axW   = 6;  axH   = 3.5;
lMarg = 1.8;  bMarg = 1.5;   % left / bottom margins (room for labels)
hGap  = 2.2;  vGap  = 2.0;   % horizontal / vertical gap between tiles
rMarg = 0.5;  tMarg = 1.2;   % right / top margins

figW = lMarg + axW + hGap + axW + rMarg;
figH = bMarg + axH + vGap + axH + tMarg;

hFig = figure;
set(hFig, 'Units', 'centimeters');
hFig.Position(2) = 5;
hFig.Position(3:4) = [figW, figH] * 1.4;

% Bottom-left corner of each tile [x y] in cm, left-to-right / top-to-bottom
axOrigins = [lMarg, bMarg + axH + vGap;     % (1,1) top-left
    lMarg + axW + hGap, bMarg + axH + vGap; % (1,2) top-right
    lMarg, bMarg;                           % (2,1) bottom-left
    lMarg + axW + hGap, bMarg];             % (2,2) bottom-right

% Plot
for iFile = 1 : length(idxFiles)

    hAx = axes('Units', 'centimeters', ...
                'Position', [axOrigins(iFile,:), axW, axH]); %#ok<LAXES>

    % Prep spktimes
    uIdx = v(iFile).units.type == 'RS';
    spktimes = cellfun(@(x) x', ...
        v(iFile).spikes.times, 'uni', false)';
    spktimes = spktimes(uIdx);
    btimes = v(iFile).burst.spktimes(uIdx);

    plot_raster(spktimes, 'PlotType', 'vertline', ...
        'lineHeight', lnH, ...
        'lineWidth', lnW, ...
        'hAx', hAx, ...
        'clr', [0 0 0], ...
        'xLim', winPlot, ...
        'spkDur', 0.0005);

    plot_raster(btimes, 'PlotType', 'vertline', ...
        'lineHeight', lnH, ...
        'lineWidth', lnW, ...
        'hAx', hAx, ...
        'clr', [1 0 0], ...
        'xLim', winPlot, ...
        'spkDur', 0.0005);

    set(gca, 'YDir', 'normal');
    title(hAx, basepaths(idxFiles(iFile)), 'Interpreter', 'none')

    if iFile == 1
        xlim([18.5 19.3])
    elseif iFile == 2
        xlim([6.1 6.9])
    elseif iFile == 3
        xlim([15.2 16])
    elseif iFile == 4
        xlim([14 14.8])
    end

    % Typography: Arial, tick labels 10 pt, axis labels 12 pt
    set(hAx, 'FontName', 'Arial', 'FontSize', 10);
    xlabel(hAx, 'Time (s)', 'FontName', 'Arial', 'FontSize', 12);
    ylabel(hAx, 'Unit No.', 'FontName', 'Arial', 'FontSize', 12);
end



%% ========================================================================
%  COLLAPSE UNIT TABLE BY DAY
%  ========================================================================

statType = 'mean';

% Variable names
varsTbl = tblLme.Properties.VariableNames;
isNum = cellfun(@(x) isnumeric(tblLme.(x)) && ~iscategorical(tblLme.(x)), varsTbl);
varsNum = varsTbl(isNum);

% Table per day
tblDay = groupsummary(tblLme, {'genotype', 'sbjID', 'day'}, statType, ...
    vartype("numeric"));

% Replace var names
tblDay(:, "GroupCount") = [];
varsTbl = tblDay.Properties.VariableNames;
isNum = cellfun(@(x) isnumeric(tblDay.(x)) && ~iscategorical(tblDay.(x)), varsTbl);
tblDay.Properties.VariableNames(isNum) = varsNum;


%% ========================================================================
%  RECOVERY METRICS (dBrst, dSngl)
%  ========================================================================

% Select variables to unstack
varsUnstack = {'frBurst', 'frSingle'};

% Unstack table (Wide format)
% Creates columns: frBurst_BSL, frBurst_BAC3, etc.
tblWide = unstack(tblDay(:, [{'genotype', 'sbjID', 'day'}, varsUnstack]), ...
    varsUnstack, 'day');

% Calculate Log-Ratios
% dBrst = log(frBurst_BAC3 / frBurst_BSL)
% dSngl = log(frSingle_BAC3 / frSingle_BSL)
tblWide.dBrst = log(tblWide.frBurst_BAC3 ./ tblWide.frBurst_BSL);
tblWide.dSngl = log(tblWide.frSingle_BAC3 ./ tblWide.frSingle_BSL);


% Plot
tblGUI_bar(tblWide, 'xVar', 'genotype', 'yVar', 'dBrst');
tblGUI_scatHist(tblWide, 'xVar', 'dSngl', 'yVar', 'dBrst', 'grpVar', 'genotype');



%% ========================================================================
%  CAG:MCU-KO - METRICS (st & bursts)
%  ========================================================================
% Three MCU-KO recordings (CAG cohort, 'ra' path key; TDT, fs ~24.4 kHz).
% Spike timing metrics, bursts and burst stats, computed with the same params
% as the wt / mcu pipeline (see mcu_wrapper). Unit classification is handled
% separately in the next section.

basepaths = mcu_basepaths('ra');
nFiles = length(basepaths);

% vars
vars = {'spikes'};

% load state vars
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Burst detection params (identical to mcu_wrapper)
isiStart = 0.006;
minSpks  = 2;
isiEnd   = isiStart * 2;
minIBI   = isiEnd;
minDur   = 0;

for iFile = 1 : nFiles

    % file
    basepath = basepaths{iFile};
    [~, basename] = fileparts(basepath);
    cd(basepath)

    % spktimes
    spktimes = v(iFile).spikes.times;

    % Spike timing metrics
    st = spktimes_metrics('spktimes', spktimes, 'sunits', [], ...
        'bins', {[0, Inf]}, 'flgForce', true, 'flgSave', true, ...
        'flgAll', false);

    % Burst detection
    burst = burst_detect(spktimes, ...
        'minSpks', minSpks, ...
        'isiStart', isiStart, ...
        'isiEnd', isiEnd * 2, ...
        'minDur', minDur, ...
        'minIBI', minIBI, ...
        'flgForce', true, 'flgSave', true, 'flgPlot', false);

    % Burst statistics
    stats = burst_stats(burst, spktimes, 'winCalc', [], 'flgSave', true);

end


%% ========================================================================
%  CAG:MCU-KO - CLASSIFICATION
%  ========================================================================
% Assign RS / FS unit types and write units.mat. Two methods, toggled by flgCE:
%   flgCE = false : my pipeline. utypes_classify (GMM on waveform features)
%                   opens utypes_gui. Inspect the scatter / waveforms, reassign
%                   points, then click "Push Units" to save units.mat.
%   flgCE = true  : fall back to CellExplorer putativeCellType (Pyramidal -> RS,
%                   Narrow Interneuron -> FS, anything else -> Other).

flgCE = false;

basepaths = mcu_basepaths('ra');
nFiles = length(basepaths);

if ~flgCE

    % My pipeline: GMM classification, then manual curation in the GUI
    fetSelect = {'Asym', 'Hpk', 'TP'};
    tblUnit = utypes_classify('basepaths', basepaths, 'fetSelect', fetSelect, ...
        'rsPrior', 0.97, 'regVal', 0.01, 'flgPlot', true);

else

    % Fallback: CellExplorer putativeCellType
    v = basepaths2vars('basepaths', basepaths, 'vars', {'cell_metrics'});

    for iFile = 1 : nFiles

        % file
        basepath = basepaths{iFile};
        [~, basename] = fileparts(basepath);

        % Unit types from CellExplorer
        ctype = v(iFile).cell_metrics.putativeCellType(:);
        nUnits = numel(ctype);
        typeNum = zeros(nUnits, 1);                                  % Other
        typeNum(contains(ctype, 'Pyramidal'))          = 1;         % RS
        typeNum(contains(ctype, 'Narrow Interneuron')) = 2;         % FS

        units = struct();
        units.clean = false(2, nUnits);
        units.clean(1, :) = (typeNum == 1)';
        units.clean(2, :) = (typeNum == 2)';
        units.type = categorical(typeNum, [0, 1, 2], {'Other', 'RS', 'FS'});
        units.date = datetime('now');
        save(fullfile(basepath, [basename, '.units.mat']), 'units')

    end
end


%% ========================================================================
%  CAG:MCU-KO - BURSTINESS COMPARISON
%  ========================================================================
% Compare baseline burstiness of the three CAG:MCU-KO mice against the wt and
% mcu baseline groups. RS units only; CAG:MCU-KO is kept as its own group.

% CAG:MCU-KO group
basepaths = mcu_basepaths('ra');
tblCag = mcu_tblVivo('basepaths', basepaths, 'presets', {'burst'}, ...
    'flgClean', false);
tblCag.genotype = addcats(tblCag.genotype, {'CAG:MCU-KO'});
tblCag.genotype(:) = 'CAG:MCU-KO';

% Reference baseline groups (wt and mcu)
basepathsRef = [mcu_basepaths('wt_bsl'), mcu_basepaths('mcu_bsl')];
tblRef = mcu_tblVivo('basepaths', basepathsRef, 'presets', {'burst'}, ...
    'flgClean', false);
tblRef.genotype = addcats(tblRef.genotype, {'CAG:MCU-KO'});

% Combine
tblBrst = [tblRef; tblCag];
tblBrst.genotype = reordercats(tblBrst.genotype, ...
    {'Control', 'MCU-KO', 'CAG:MCU-KO'});

% Compare (switch yVar in the GUI for frBurst, br, bSize, etc.)
tblGUI_bar(tblBrst, 'yVar', 'pBurst', 'xVar', 'genotype', 'grpVar', 'unitType');
