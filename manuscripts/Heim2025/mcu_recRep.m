% TO DO LIST
% # xlim of fr graph to -24 : +24
% # horizontal bar on fr graph automatic update

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mname = 'lh96';
% mname = 'lh133';
mname = 'lh132';
% mname = 'lh107';

forceL = false;
forceL = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mouse params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch mname
    case 'lh96'
        fname{1} = 'F:\Data\lh96\lh96_220121_090213';
        fname{2} = 'F:\Data\lh96\lh96_220125_090041';

        rawCh{1} = [1 : 4];
        rawCh{2} = [9 : 11, 13];

        xStart{1}(1) = 6.9 * 60 * 60;
        xStart{1}(2) = xStart{1}(1);
        xStart{1}(3) = xStart{1}(1);
        xStart{1}(4) = 25200;
        xStart{1}(5) = 36012;
        xStart{1}(6) = xStart{1}(4) + 5;
        xStart{1}(7) = xStart{1}(5) + 1.7;

        xStart{2}(1) = 4.4 * 60 * 60;
        xStart{2}(2) = xStart{2}(1);
        xStart{2}(3) = xStart{2}(1);
        xStart{2}(4) = 16215;
        xStart{2}(5) = 26500;
        xStart{2}(6) = xStart{2}(4) + 2.4;
        xStart{2}(7) = xStart{2}(5) + 0.8;

    case 'lh133'
        fname{1} = 'F:\Data\lh133\lh133_230414_090645';
        fname{2} = 'F:\Data\lh133\lh133_230418_094147';

        rawCh{1} = [5 : 8];
        rawCh{2} = [5 : 8];

        xStart{1}(1) = 4.2 * 60 * 60;
        xStart{1}(2) = xStart{1}(1);
        xStart{1}(3) = xStart{1}(1);
        xStart{1}(4) = 16200;
        xStart{1}(5) = 45000;
        xStart{1}(6) = xStart{1}(4) + 5.2;
        xStart{1}(7) = xStart{1}(5) + 2.2;

        xStart{2}(1) = 2.2 * 60 * 60;
        xStart{2}(2) = xStart{2}(1);
        xStart{2}(3) = xStart{2}(1);
        xStart{2}(4) = 8800;
        xStart{2}(5) = 36700;
        xStart{2}(6) = xStart{2}(4) + 1.8;
        xStart{2}(7) = xStart{2}(5) + 4.4;

    case 'lh107'
        fname{1} = 'F:\Data\lh107\lh107_220519_091300';
        fname{2} = 'F:\Data\lh107\lh107_220523_102100';

        rawCh{1} = [1 : 4];
        rawCh{2} = [13 : 16];

        xStart{1}(1) = 4 * 60 * 60;
        xStart{1}(2) = xStart{1}(1);
        xStart{1}(3) = xStart{1}(1);
        xStart{1}(4) = 25200;
        xStart{1}(5) = 36012;
        xStart{1}(6) = xStart{1}(4) + 5;
        xStart{1}(7) = xStart{1}(5) + 1.7;

        xStart{2}(1) = 3.5 * 60 * 60;
        xStart{2}(2) = xStart{2}(1);
        xStart{2}(3) = xStart{2}(1);
        xStart{2}(4) = 16215;
        xStart{2}(5) = 26500;
        xStart{2}(6) = xStart{2}(4) + 2.4;
        xStart{2}(7) = xStart{2}(5) + 0.8;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general plotting params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure panel and tile arrangement
panels2plot =["hypnogram", "spec", "emg", "raster", "raster", "raw", "raw"];
npanels = length(panels2plot);

% xDur determines the window size for each plot and xStart determined the
% timepoint where the plot starts. the actual duration that is plotted is
% 3-times xDur to allow for manual scrolling through the data. All window
% params should be given in seconds. Panels that share the same duration
% and start point will be x-linked
xDur(1) = 10 * 60 * 60;
xDur(2) = xDur(1);
xDur(3) = xDur(1);
xDur(4) = 10;
xDur(5) = xDur(4);
xDur(6) = 0.5;
xDur(7) = xDur(6);

% set relationships between panels for the lines that represent zoom in
% views. for example, xRelations{3} = 4, will add a line to panel 3 that
% represents the xlim of panel 4
xRelations = cell(npanels, 1);
xRelations{3} = [4, 5];
xRelations{4} = 6;
xRelations{5} = 7;

% general params
saveFig = false;

% downsample eeg
fs_eeg = 250;

% generate figure
fh = figure;
tlayout = [6, 4];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot fr of entire experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axh = nexttile(th, 1, [1, 4]);

% get fr data
[expData, info] = sessions_catVarTime('mname', mname,...
    'dataPreset', 'fr', 'graphics', false, 'dataAlt', 1,...
    'basepaths', {}, 'xTicksBinsize', 24, 'markRecTrans', false,...
    'axh', axh, 'saveFig', false);
expData = mean(expData, 2, 'omitnan');


xData = info.x_data;
xt = info.x_ticks;


% replace nan values
if strcmp(mname, 'lh133')
    [tmp, info_tmp] = sessions_catVarTime('mname', 'lh107',...
        'dataPreset', 'fr', 'graphics', false, 'dataAlt', 1,...
        'basepaths', {}, 'xTicksBinsize', 24, 'markRecTrans', false,...
        'axh', axh, 'saveFig', false);
    tmpx = info_tmp.x_data;
    nanIdx = isnan(expData');
    gapStarts = find(diff([false, nanIdx, false]) == 1);
    gapEnds = find(diff([false, nanIdx, false]) == -1) - 1;
    gapSizes = gapEnds - gapStarts + 1;
    [maxGapSize, maxGapIdx] = max(gapSizes);
    gapIdx = gapStarts(maxGapIdx):gapEnds(maxGapIdx);
    fillIdx = find(tmpx == xData(gapIdx(1))) + 2200;
    expData(gapIdx) = tmp(fillIdx : fillIdx + length(gapIdx) - 1);
end

% plot
plotData = movmean(expData, 13, 1);
plot(axh, xData, plotData, 'k', 'LineWidth', 0.5)
legend('off')
xlabel('')
ylabel('MFR (Hz)')
yLimit = [0 3.5];
ylim(yLimit)
hold on

% load session data from entire experiment
varsFile = ["datInfo"; "session"];
varsName = ["datInfo"; "session"];
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"]);

% get time point of op on / off. this can be very specific for each mouse
recStart(1) = xt(2);
recStart(2) = xt(6);
t_on = recStart(1) + v(2).session.general.timepnt;
t_off = recStart(2) + v(6).session.general.timepnt;

% set xticks and labels
dt = 60 * 60 * 6;
newTicks = unique([t_on : -dt : 1, t_on : dt : xData(end)]);
xticks(newTicks)
switch mname
    case 'lh96'
        xLabels = string([-30 : 6 : 6 * length(newTicks) - 30]);
        xLimit = [newTicks(2), newTicks(26)];
    
    case 'lh133'
        xLabels = string([-24 : 6 : 6 * length(newTicks) - 30]);
        xLimit = [newTicks(1), newTicks(25)];

    case 'lh107'
        xLabels = string([-24 : 6 : 6 * length(newTicks) - 30]);
        xLimit = [newTicks(1), newTicks(25)];
end
xticklabels(xLabels)
xlim(xLimit)

% add dashed line at time points of interest (op on / off)
plot([t_on, t_on], yLimit, '--r', 'LineWidth', 2)
plot([t_off, t_off], yLimit, '--r', 'LineWidth', 2)

% add horizontal lines on fr experiment graph 
for ifile = 1 : 2
    xwin = [xStart{ifile}(1), xStart{ifile}(1) + xDur(1)] + recStart(ifile);
    plot(axh, xwin, [0, 0], 'b', 'LineWidth', 10);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot recording - WT BAC ON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
if forceL
    for ifile = 1 : 2
        [emg{ifile}, spec{ifile}, stateEpochs{ifile}] = recRep_prep('basepath', fname{ifile},...
            'panels2plot', panels2plot, 'fs_eeg', fs_eeg);
    end
end

% plot
for ifile = 1 : 2
       
    % get time point of interest
    basepath = fname{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    load([basename, '.session.mat'])
    tLine = session.general.timepnt;

    % organization of panels across tiles
    if ifile == 1
        panelTiles{1} = [5, 6];
        panelTiles{2} = [9, 10];
        panelTiles{3} = [13, 14];
        panelTiles{4} = [17];
        panelTiles{5} = [18];
        panelTiles{6} = [21];
        panelTiles{7} = [22];

        xOffset = 0;
    else
        panelTiles{1} = [7, 8];
        panelTiles{2} = [11, 12];
        panelTiles{3} = [15, 16];
        panelTiles{4} = [19];
        panelTiles{5} = [20];
        panelTiles{6} = [23];
        panelTiles{7} = [24];

        xOffset = t_off - t_on;
    end

    % plot
    recRep_plot('basepath', basepath,...
        'xDur', xDur, 'xStart', xStart{ifile}, 'panels2plot', panels2plot,...
        'spec', spec{ifile}, 'emg', emg{ifile}, 'stateEpochs', stateEpochs{ifile},...
        'saveFig', saveFig, 'rawCh', rawCh{ifile},...
        'xRelations', xRelations, 'tlayout', tlayout, 'panelTiles', panelTiles,...
        'tLine', tLine, 'th', th, 'adjustT', true, 'xOffset', xOffset);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cleanup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove redundant labels
axh = findobj(th, 'Type', 'axes');
ylabel(axh(3) , '')     % raster 4
ylabel(axh(4) , '')     % raster 3
ylabel(axh(5) , '')     % emg 2
ylabel(axh(6) , '')     % spec 2
ylabel(axh(1) , '')     % raw 4
ylabel(axh(2) , '')     % raw 3
ylabel(axh(10) , '')    % raster 2
ylabel(axh(8) , '')     % raw 2
% ylabel(axh(9) , '')   % raw 1
% ylabel(axh(11) , '')    % raster 1

lgdh = findobj(th, 'Type', 'legend');
delete(lgdh(1))

title(th, '')

allAxes = findall(gcf, 'type', 'axes');
set(allAxes, 'FontSize', 14);
















