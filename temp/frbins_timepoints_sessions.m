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
clear frBins
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
 
    % timebins
    fileinfo = dir([basename, '.dat']);
    recLen = floor(fileinfo.bytes / 2 / nchans / fs);
    csec = floor(cumsum(v(isession).datInfo.nsamps / fs));
    [~, pntIdx] = min(abs(csec - 5.5 * 60 * 60));
    timepoints = csec(pntIdx);
    % timepoints = v(isession).datInfo.nsec;
    chunks = n2nchunks('n', recLen, 'nchunks', 4, 'timepoints', timepoints);
    
    clear tmp_units
    tmp_units(1, :) = selectUnits(v(isession).spikes, v(isession).cm,...
        v(isession).fr, 1, grp, frBoundries, 'pyr');
    tmp_units(2, :) = selectUnits(v(isession).spikes, v(isession).cm,...
        v(isession).fr, 1, grp, frBoundries, 'int');
    units = [units, tmp_units];
    nunits(isession) = length(tmp_units);
    
    frBins(isession, :) = fr_timebins('basepath', pwd,...
        'forceA', true, 'graphics', true,...
        'timebins', chunks, 'saveVar', true);
    
end
units = logical(units);

sratio = [];
mfrWake = [];
mfrNrem = [];
ridx = [1, 4];
for isession = 1 : nsessions
    clear sratio_temp mfrWake_temp mfrNrem_temp
    for iwin = 1 : size(chunks, 1)
        sratio_temp(iwin, :) = squeeze(frBins(isession, iwin).states.ratio(ridx(1), ridx(2), :));
        mfrWake_temp(iwin, :) = mean(frBins(isession, iwin).states.fr{1}, 2);
        mfrNrem_temp(iwin, :) = mean(frBins(isession, iwin).states.fr{4}, 2);
    end
    sratio = [sratio, sratio_temp];
    mfrWake = [mfrWake, mfrWake_temp];
    mfrNrem = [mfrNrem, mfrNrem_temp];
end

yLimit = [min([sratio], [], 'all'), max([sratio], [], 'all')];
tbins_txt = {'0-3ZT', '3-6ZT', '6-9ZT', '9-12ZT',...
    '12-15ZT', '15-18ZT', '18-21ZT', '21-24ZT'};
        
% graphics
fh = figure;
cfg = as_loadConfig();

subplot(3, 2, 1)
dataMat = mfrWake(:, units(1, :));
plot_boxMean('dataMat', dataMat', 'clr', cfg.colors{1})
ylabel('MFR WAKE')
subtitle('RS units')
xticklabels(tbins_txt)

subplot(3, 2, 2)
dataMat = mfrWake(:, units(2, :));
plot_boxMean('dataMat', dataMat', 'clr', cfg.colors{1})
ylabel('MFR WAKE')
subtitle('FS units')
xticklabels(tbins_txt)

subplot(3, 2, 3)
dataMat = mfrNrem(:, units(1, :));
plot_boxMean('dataMat', dataMat', 'clr', cfg.colors{4})
ylabel('MFR NREM')
xticklabels(tbins_txt)

subplot(3, 2, 4)
dataMat = mfrNrem(:, units(2, :));
plot_boxMean('dataMat', dataMat', 'clr', cfg.colors{4})
ylabel('MFR NREM')
xticklabels(tbins_txt)

subplot(3, 2, 5)
dataMat = sratio(:, units(1, :));
plot_boxMean('dataMat', dataMat', 'clr', 'k')
ylabel({'WAKE - NREM /', 'WAKE + NREM'})
ylim(yLimit)
xticklabels(tbins_txt)

subplot(3, 2, 6)
dataMat = sratio(:, units(2, :));
plot_boxMean('dataMat', dataMat', 'clr', 'k')
ylabel({'WAKE - NREM /', 'WAKE + NREM'})
ylim(yLimit)
xticklabels(tbins_txt)



figname = 'fr_time';
figpath = 'D:\Google Drive\PhD\Slutsky';
figname = fullfile(figpath, 'fr_temp');
export_fig(figname, '-tif', '-transparent', '-r300')



