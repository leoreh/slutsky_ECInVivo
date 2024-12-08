% IED_wrapper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'forceDef', true, 'forceL', true, 'saveVar', true);      
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
[~, basename] = fileparts(basepath);

% load sSig or just the sig info
sSig = load([basename, '.sleep_sig.mat']);
load([basename, '.sleep_sig.mat'], 'fs');



% consider filtering the data
% sig = highpass(sig, 1 , fs);
% sig = lowpass(sig1, 800 , fs);
% emg = lfp.data(:,4);
% emg = bandstop(emg,[48 52],fs); %only when needed
% emg = highpass(sig, 1 , fs); % sometimes higher
% sig = lowpass(sig1, 100 , fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% see Gelinas (Buzsaki), Nat. Med., 2016 for a different pipeline

thr_mv = 0;
thr_Z  = 6;
ied = IED.detect("sig", sSig.eeg, "emg", sSig.emg, "fs", fs, "thr",...
    [thr_Z, thr_mv], "thrDir", "both");

ied.sig = sSig.eeg;
ied.emg = sSig.emg;
ied = IED.curate(ied, "saveVar", true, "basepath", basepath,...
    "basename",basename);

% analyze after curation:
bin_size = ied.fs*60;   % bin for rate, in [samples]
smooth_win = 1 ;        % size of movemean window, in [bins]
clip_marg = 0.05;       % how much time to show around the event
ied = IED.analyze(ied, "binsize", bin_size, "smf", smooth_win,...
    "marg", clip_marg, "saveFig", false, "saveVar", true, "sig", sSig.eeg);



load([basename, '.sleep_states.mat']);
labelsmanfile = [basename, '.sleep_labelsRev.mat'];
AccuSleep_viewer(sSig, ss.labels, labelsmanfile)



% cut EDs from sig
tstamps = [1 : length(sSig.eeg)]' / fs;
durWin = [-150 150] / 1000;
nbinsMap = floor(fs * diff(durWin) / 2) * 2 + 1; % must be odd
centerBin = ceil(nbinsMap / 2);
[r, i] = Sync([tstamps sSig.eeg], ied.pos(ied.accepted) / fs, 'durations', durWin);
edWv = SyncMap(r, i, 'durations', durWin,...
    'nbins', nbinsMap, 'smooth', 0);


fh = figure;
plot(edWv')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over mice 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
bin_size = ied.fs*60;   % bin for rate, in [samples]
smooth_win = 1 ;        % size of movemean window, in [bins]
clip_marg = 0.05;       % how much time to show around the event
fs = 1250;

% files
mname = 'wt_bsl';
basepaths = [mcu_sessions(mname)];
nfiles = length(basepaths);

% inspect
for ifile = 1 : nfiles

    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    load([basename, '.sleep_sig.mat'], 'eeg');
    load([basename, '.ied.mat'], 'ied');

    if sum(ied.accepted) > 0
        ied = IED.analyze(ied, "binsize", bin_size, "smf", smooth_win,...
            "marg", clip_marg, "saveFig", false, "saveVar", true, "sig", eeg);
    end

end


% analyze in relation to states

% files
mname = 'mcu_bsl';
basepaths = [mcu_sessions(mname)];
nfiles = length(basepaths);

% analyze 
stateIdx = cell(nfiles, 1);
for ifile = 1 : nfiles

    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    load([basename, '.ied.mat'], 'ied');
    load([basename, '.sleep_states.mat'], 'ss');

    tstamps = round(ied.pos(ied.accepted) / ied.fs);
    stateIdx{ifile} = ss.labels(tstamps);

end
sMat = cell2padmat(stateIdx, 2);

% plot
sstates = [1 : 5];
sums = zeros(length(sstates), size(sMat, 2));
for i = 1:size(sMat, 2)
    for j = 1:length(sstates)
        sums(j, i) = sum(sMat(:, i) == sstates(j), 'omitnan');
    end
end

fh = figure;
set(fh, 'DefaultAxesFontSize', 16);

bh = bar(sums', 'stacked');
for ib = 1 : length(bh)
    bh(ib).FaceColor = 'flat';
    bh(ib).CData = ss.info.colors{ib};
end
mnames = get_mname(basepaths);
legend(ss.info.names(sstates))
xticks(1 : length(mnames))
xticklabels(mnames)
ylabel('IEDs per State');
ylim([0 12])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze ied 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params 
sstates = [1 : 5];
fs = 1250;

% file
basepath = pwd;
cd(basepath)
[~, basename] = fileparts(basepath);

% load data
sSig = load([basename, '.sleep_sig.mat']);
load([basename, '.sleep_states.mat'], 'ss');
load([basename, '.ied.mat'], 'ied');

% organize data
spec.s = sSig.spec;
spec.freq = sSig.spec_freq;
spec.tstamps = sSig.spec_tstamps;

% vars
tstamps = ied.pos(ied.accepted) / fs;

% cut EDs from sig and emg
sig_tstamps = [1 : length(sSig.eeg)]' / fs;
durWin = [-200 200] / 1000;
nbinsMap = floor(fs * diff(durWin) / 2) * 2 + 1; % must be odd
centerBin = ceil(nbinsMap / 2);
[r, i] = Sync([sig_tstamps sSig.eeg], tstamps, 'durations', durWin);
[edWv, tbins] = SyncMap(r, i, 'durations', durWin, 'nbins', nbinsMap, 'smooth', 0);
[r, i] = Sync([sig_tstamps sSig.emg], tstamps, 'durations', durWin);
edEmg = SyncMap(r, i, 'durations', durWin, 'nbins', nbinsMap, 'smooth', 0);

% open figure
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlay = [4, 3];
th = tiledlayout(tlay(1), tlay(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, basename, 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 20);

% spectrogram
axh1 = nexttile(th, 1, [1, 3]); cla; hold on
plot_spec(spec, 'ch', 1, 'logfreq', false, 'saveFig', false,...
    'axh', axh1, 'xtime', 1)
axis tight
yLimit = ylim;
xticks([]);
xlabel('')

% hypnogram and emg
axh2 = nexttile(th, tlay(2) + 1, [1, 3]); cla; hold on
emg = sSig.emg_rms;
plot([1 : length(emg)], emg, 'k', 'LineWidth', 0.5)
yLimit = ylim;
plot_hypnogram('stateEpochs', ss.stateEpochs,...
    'clr', ss.info.colors(sstates), 'axh', axh2, 'sstates', [1 : length(sstates)],...
    'yshift', 1)
xval = [3600 : 3600 : length(emg)];
xticks([]);
set(axh2, 'YTickMode', 'auto')
set(axh2, 'YColor', 'k')
yticks(axh2, [])
ylabel(axh2, 'EMG')

% ied rate
axh3 = nexttile(th, tlay(2) * 2 + 1, [1, 3]); cla; hold on
[rate_ied, ~, rate_tstamps] = times2rate(tstamps, 'binsize', 60 * 5,...
    'winCalc', [0, Inf], 'c2r', true);
plot(rate_tstamps, rate_ied)
yLimit = ylim;
scatter(tstamps, ones(length(tstamps), 1) * yLimit(2), 4, 'filled', 'r')
xval = [3600 : 3600 : length(emg)];
xticks(xval);
xticklabels(string(xval / 3600))
xlabel('Time [Hr]')
ylabel('IED Rate (Hz)')

linkaxes([axh1, axh2, axh3], 'x')
axis tight

axh = nexttile(th, 10, [1, 1]); cla; hold on
plot(tbins, edWv)
xlabel('Time (s)')
ylabel('Voltage')

axh = nexttile(th, 11, [1, 1]); cla; hold on
plot(tbins, edEmg)
xlabel('Time (s)')
ylabel('EMG')

axh = nexttile(th, 12, [1, 1]); cla; hold on
coff = [prctile(edEmg(:), 20), prctile(edEmg(:), 80)];
PlotColorMap(edEmg, 1, 'bar','on', 'cutoffs', coff, 'x', tbins);
ylabel('IED No.')
xlabel(axh, 'Time (s)')
title('EMG')


