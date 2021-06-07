
% inscopix
% load data
fname = dir('tra*');
load(fname.name)
wake = trace1;
fname = dir('sl*');
load(fname.name)
sleep = trace1(1 : length(wake));

sleep = bz_NormToRange(sleep, [0 1]);
wake = bz_NormToRange(wake, [0 1]);

fs = 10;
tstamps = [1 / fs : 1 / fs : length(wake) / fs];
yLimit = sort([0 1]);

figure
subplot(2, 1, 1)
plot(tstamps, wake)
ylim(yLimit)
box off
yticks([])
xticks([])
title('wt')
subplot(2, 1, 2)
plot(tstamps, sleep)
ylim(yLimit)
box off
title('app')



% lfp

% arrange data
basepath = 'F:\Data\Colleagues\SS';
cd(basepath)
fname = dir('WT*');
load(fname.name)
wt = lfp;
fname = dir('APP*');
load(fname.name)
app = lfp;

% params
fs = lfp.fs;
idx = 5 * 60 * fs : 7 * 60 * fs;
yLimit = sort([-2.5 1]);

figure
subplot(2, 1, 1)
plot(wt.timestamps(idx), wt.data(idx))
ylim(yLimit)
box off
yticks([])
xticks([])
title('wt')
subplot(2, 1, 2)
plot(app.timestamps(idx), app.data(idx))
ylim(yLimit)
box off
title('app')
