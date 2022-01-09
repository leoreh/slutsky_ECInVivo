% ket_sessions_psd

forceL = true;
normFlag = false;
frBoundries = [0.1, Inf; 0.1, Inf];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% local ket
basepaths = [{'F:\Data\Processed\lh96\lh96_211202_070500'},...
    {'K:\Data\lh95\lh95_210825_080400'},...
    {'K:\Data\lh99\lh99_211219_085802'},...
    {'G:\Data\lh93\lh93_210813_110609'}];

% baclofen
basepaths = [{'F:\Data\Processed\lh96\lh96_211207_071500'},...
    {'K:\Data\lh99\lh99_211220_091903'},...
    {'G:\Data\lh98\lh98_211220_104619'}];
basepaths = basepaths(1 : 2);

% local acsf
basepaths = [{'K:\Data\lh99\lh99_211218_090630'},...
    {'F:\Data\Processed\lh96\lh96_211201_070100'},...
    {'K:\Data\lh95\lh95_210824_083300'},...
    {'G:\Data\lh93\lh93_210811_102035'}];

% ket 10 mg/kg i.p.
basepaths = [{'K:\Data\lh99\lh99_211224_084528'},...
{'G:\Data\lh98\lh98_211224_084528'},...    
{'F:\Data\Processed\lh96\lh96_211206_070400'},...
    {'G:\Data\lh81\lh81_210204_190000'}];
basepaths = basepaths(1 : 3);

% ket 60 mg/kg i.p.
basepaths = [{'G:\Data\lh84\lh84_210504_225608'},...
    {'G:\Data\lh81\lh81_210206_190000'}];



% load vars from each session
varsFile = ["datInfo"; "session"; "psdBins.mat"];
varsName = ["datInfo"; "session"; "psdBins"];
if ~exist('v', 'var') || forceL
    v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName);
end
nsessions = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
faxis = 0.2 : 0.2 : 120;
stateIdx = 1;
freqIdx = faxis > 20 & faxis < 100;
faxis(freqIdx);

dataMat = nan(nsessions * 8, sum(freqIdx));
for isession = 1 : nsessions
    
    % session params
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    
        tbins_txt = {'0-3ZT', '3-6ZT', '6-9ZT', '9-12ZT',...
            '12-15ZT', '15-18ZT', '18-21ZT', '21-24ZT'};
        psdBins = psd_states_timebins('basepath', pwd,...
            'chEeg', [], 'forceA', false, 'graphics', true,...
            'timebins', chunks, 'saveVar', true,...
            'sstates', [1, 4, 5], 'tbins_txt', tbins_txt);
    
end

cnt = 1;
for ibin = 1 : 8
    for isession = 1 : nsessions
        dataVec = squeeze(v(isession).psdBins.psdLfp(ibin, stateIdx, freqIdx));
        dataVec = dataVec / sum(dataVec);
        dataMat(cnt, :) = dataVec;           
        cnt = cnt + 1;
    end
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
plot([tidx tidx], ylim, '--k')
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
plot([tidx tidx], ylim, '--k')
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
plot([tidx tidx], yLimit, '--k')
ax = gca;
set(ax.YAxis(1), 'color', 'b')
set(ax.YAxis(2), 'color', 'r')
xlabel('ZT Time [h]')
set(gca, 'box', 'off')
linkaxes([sb1, sb2, sb3], 'x')
xlim([0 24])
xticks([0 : 3 : 24])
figname = 'fr_time';




