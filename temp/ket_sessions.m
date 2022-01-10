% ket_sessions

forceL = true;
normFlag = false;
frBoundries = [0.1, Inf; 0.1, Inf];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% local mk801
basepaths = [{'F:\Data\Processed\lh96\lh96_211205_072000'}];

% local acsf
basepaths = [{'K:\Data\lh99\lh99_211218_090630'},...
    {'F:\Data\Processed\lh96\lh96_211201_070100'},...
    {'K:\Data\lh95\lh95_210824_083300'},...
    {'G:\Data\lh93\lh93_210811_102035'}];

% baclofen
basepaths = [{'F:\Data\Processed\lh96\lh96_211207_071500'},...
    {'K:\Data\lh99\lh99_211220_091903'},...
    {'G:\Data\lh98\lh98_211220_104619'}];
basepaths = basepaths(1 : 2);

% local ket
basepaths = [{'I:\lh96\lh96_211126_072000'},...
    {'F:\Data\Processed\lh96\lh96_211202_070500'},...
    {'K:\Data\lh95\lh95_210825_080400'},...
    {'K:\Data\lh99\lh99_211219_085802'}];

% load vars from each session
varsFile = ["fr"; "sr"; "spikes"; "st_metrics"; "swv_metrics";...
    "cell_metrics"; "sleep_states"; "datInfo"; "session"];
varsName = ["fr"; "sr"; "spikes"; "st"; "swv"; "cm"; "ss";...
    "datInfo"; "session"];
if ~exist('v', 'var') || forceL
    v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName);
end
nsessions = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

units = [];
injIdx = [];
tLen = [];
for isession = 1 : nsessions
    
    % session params
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    
    fs = v(isession).session.extracellular.sr;
    fsLfp = v(isession).session.extracellular.srLfp;
    spkgrp = v(isession).session.extracellular.spikeGroups.channels;
    nchans = v(isession).session.extracellular.nChannels;
    
    if contains(basename, 'lh99')
        grp = [1, 3 : 4, 7];
    else
        grp = [];
    end
    
    % plot fr vs. time
    %     plot_FRtime_session('basepath', pwd, 'grp', grp,...
    %         'frBoundries', [0.01 Inf; 0.01 Inf], 'muFlag', false, 'saveFig', false,...
    %         'dataType', 'strd')
    % %
    %     plot_FRtime_session('basepath', pwd, 'grp', grp,...
    %         'frBoundries', [0.01 Inf; 0.01 Inf], 'muFlag', false, 'saveFig', false,...
    %         'dataType', 'norm')
    
    % timebins
    fileinfo = dir([basename, '.dat']);
    recLen = floor(fileinfo.bytes / 2 / nchans / fs);
    csec = floor(cumsum(v(isession).datInfo.nsamps / fs));
    [~, pntIdx] = min(abs(csec - 5.5 * 60 * 60));
    timepoints = csec(pntIdx);
    % timepoints = v(isession).datInfo.nsec;
    chunks = n2nchunks('n', recLen, 'nchunks', 8, 'timepoints', timepoints);
    
    [~, injIdx(isession)] = min(abs(v(isession).fr.tstamps - timepoints));
    tLen(isession) = length(v(isession).fr.tstamps);
    
    clear tmp_units
    tmp_units(1, :) = selectUnits(v(isession).spikes, v(isession).cm,...
        v(isession).fr, 1, grp, frBoundries, 'pyr');
    tmp_units(2, :) = selectUnits(v(isession).spikes, v(isession).cm,...
        v(isession).fr, 1, grp, frBoundries, 'int');
    units = [units, tmp_units];
    nunits(isession) = length(tmp_units);
    
end

units = logical(units);

maxIdx = max(injIdx);
frMat = nan(length(units), maxIdx + max(tLen));

cnt = 1;
dataType = 'strd';
for isession = 1 : nsessions
    startIdx = maxIdx - injIdx(isession) + 1;
    frIdx = startIdx : tLen(isession) + startIdx - 1;
    frMat(cnt : cnt + nunits(isession) - 1, frIdx) = v(isession).fr.(dataType);
    cnt = cnt + nunits(isession);
end

% noramlize
if normFlag
    frNorm = frMat ./ mean(frMat(:, 1 : maxIdx), 2, 'omitnan');
    frNorm(frNorm == Inf) = nan;
    ytxt = 'Norm Firing Rate';
else
    frNorm = frMat;
    ytxt = 'Firing Rate [Hz]';
end

xidx = [1 : length(frNorm)] / 60;
tidx = maxIdx / 60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.3)
set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.5)
set(groot, 'DefaultAxesFontSize', 14)
fh = figure;
% ---------------------------------------------------------------------
% individual cells on a log scale
sb1 = subplot(3, 1, 1);
hold on
% rs
ph = plot(xidx, frNorm(units(1, :), :), 'b', 'LineWidth', 1);
alphaIdx = linspace(1, 0.2, length(ph));
clrIdx = linspace(0.2, 0.6, length(ph));
[~, mfr_order] = sort(mean(frNorm(units(1, :), :), 2, 'omitnan'));
for iunit = 1 : length(ph)
    ph(iunit).Color(1) = clrIdx(mfr_order(iunit));
    ph(iunit).Color(2) = alphaIdx(mfr_order(iunit));
    ph(iunit).Color(4) = alphaIdx(mfr_order(iunit));
end
set(gca, 'YScale', 'log')
plot([tidx tidx], ylim, '--k', 'LineWidth', 2)
axis tight
ylabel(ytxt)
set(gca, 'box', 'off')
legend(sprintf('RS = %d su', sum(units(1, :))))
xticks([0 : 3 : 24])

% fs
sb2 = subplot(3, 1, 2);
hold on
ph = plot(xidx, frNorm(units(2, :), :), 'r', 'LineWidth', 1);
alphaIdx = linspace(1, 0.2, length(ph));
clrIdx = linspace(0.2, 0.6, length(ph));
[~, mfr_order] = sort(mean(frNorm(units(2, :), :), 2, 'omitnan'));
for iunit = 1 : length(ph)
    ph(iunit).Color(3) = clrIdx(mfr_order(iunit));
    ph(iunit).Color(2) = alphaIdx(mfr_order(iunit));
    ph(iunit).Color(4) = alphaIdx(mfr_order(iunit));
end
set(gca, 'YScale', 'log')
plot([tidx tidx], ylim, '--k', 'LineWidth', 2)
axis tight
ylabel(ytxt)
set(gca, 'box', 'off')
legend(sprintf('FS = %d su', sum(units(2, :))));
xticks([0 : 3 : 24])

% ---------------------------------------------------------------------
% mean per cell class on a linear scale
sb3 = subplot(3, 1, 3);
hold on
plot(xidx, mean(frNorm(units(1, :), :), 1, 'omitnan'), 'b', 'LineWidth', 2)
ylabel(['RS ' ytxt])
yyaxis right
plot(xidx, mean(frNorm(units(2, :), :), 1, 'omitnan'), 'r', 'LineWidth', 2)
ylabel(['FS ' ytxt])
yLimit = ylim;
plot([tidx tidx], yLimit, '--k', 'LineWidth', 2)
ax = gca;
set(ax.YAxis(1), 'color', 'b')
set(ax.YAxis(2), 'color', 'r')
xlabel('ZT Time [h]')
set(gca, 'box', 'off')
linkaxes([sb1, sb2, sb3], 'x')
xlim([0 24])
xticks([0 : 3 : 24])

figname = 'fr_time';
figpath = 'D:\Google Drive\PhD\Slutsky';
figname = fullfile(figpath, 'fr_temp');
export_fig(figname, '-tif', '-transparent', '-r300')



