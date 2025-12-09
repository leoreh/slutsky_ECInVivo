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
    'E:\Data\Colleagues\RA\hDLX_Gq_WT5\040221_0655_24hr',...
    'E:\Data\Colleagues\RA\hDLX_Gq_WT2\200820_bslDay1',...
    'E:\Data\Colleagues\RA\hDLX_Gq_Tg\210820_bslDay2Raw2',...
    'E:\Data\Colleagues\RA\tg4_210730_155700',...
    };

% MCU-KO (Refaela)
paths{5} = {...
    'E:\Data\Colleagues\RA\MCU1_080621_0930',...
    'E:\Data\Colleagues\RA\MCU2_080621_0930',...
    'E:\Data\Colleagues\RA\MCU3_20211203_084720',...
    'E:\Data\Colleagues\RA\MCU4_211220_0834',...
    'E:\Data\Colleagues\RA\MCU5_220322_1906',...
    };

basepaths = [paths{:}];
nPaths = length(basepaths);


%% ========================================================================
%  ANAYLZE (unit features)
%  ========================================================================

flgSave = true;

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


%% ========================================================================
%  CLASSIFY
%  ========================================================================

% Re-load data
vars = {'swv_metrics', 'st_metrics', 'fr', 'units'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Create table of features for classification
fetTbl = utypes_features('basepaths', basepaths, 'flgPlot', true);

% Classify
unitType = utypes_classify('basepaths', basepaths, 'altClassify', 3,...
    'flgSave', false, 'fetTbl', fetTbl);


%% ========================================================================
%  INSPECT
%  ========================================================================

% Separation by features
fh = figure;
xFet = 'tp';
yFet = 'lidor';
zFet = 'asym';
hAx = subplot(1, 2, 1); hold on
plot_utypes('basepaths', basepaths, 'fetTbl', fetTbl,...
    'plotType', 'scatter3', 'xFet', xFet, 'yFet', yFet, 'zFet', zFet,...
    'unitIdx', unitType, 'hAx', hAx)

hAx = subplot(1, 2, 2); hold on
plot_utypes('basepaths', basepaths, 'flgRaw', false,...
    'plotType', 'wv', 'unitIdx', unitType, 'hAx', hAx)


%% ========================================================================
%  VISUALIZATION
%  ========================================================================

