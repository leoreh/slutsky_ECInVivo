%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[mousepath, basename] = fileparts(basepath);

% load and assign
varArray = getSessionVars('dirnames', {basename}, 'mousepath', mousepath,...
    'sortDir', false);
assignVars(varArray, 1)

% params
fs = session.extracellular.sr;
fsLfp = session.extracellular.srLfp;

% recording window
csum_blockDur = cumsum(datInfo.nsec);
dark_phase_start = csum_blockDur(2);
inj_time = csum_blockDur(1);
recWin = [0, inj_time;...
    inj_time + 1, dark_phase_start];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distribution of ripple rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rate --------------------------------------------------------------------
fh = figure;
histBins = 50;
lgndstr = {'Baseline', 'aCSF'};

% all ripples
for iwin = 1 : size(recWin, 2)
    idx_win{iwin} = InIntervals(ripp.rate.tstamps, recWin(iwin, :));
end

subplot(2, 2, 1)
plotHist(ripp.rate.rate(idx_win{1}), 'k', histBins)
hold on
plotHist(ripp.rate.rate(idx_win{2}), 'b', histBins)
xlabel('Ripple Rate [Hz]')
subtitle('All Ripples')
legend(lgndstr)

% only in nrem
for iwin = 1 : size(recWin, 2)
    idx_win{iwin} = InIntervals(ripp.states.tstamps{4}, recWin(iwin, :));
end

subplot(2, 2, 2)
plotHist(ripp.states.rate{4}(idx_win{1}), 'k', histBins)
hold on
plotHist(ripp.states.rate{4}(idx_win{2}), 'b', histBins)
xlabel('Ripple Rate [Hz]')
subtitle('NREM Ripples')
legend(lgndstr)

% amplitude ---------------------------------------------------------------
for iwin = 1 : size(recWin, 2)
    idx_rec = InIntervals(ripp.peakPos, recWin(iwin, :));
    idx_nrem = InIntervals(ripp.peakPos, ss.stateEpochs{4});
    idx_win{iwin} = idx_rec & idx_nrem;
end

histBins = 200;
subplot(2, 2, 3)
plotHist(ripp.peakAmp(idx_win{1}), 'k', histBins)
hold on
plotHist(ripp.peakAmp(idx_win{2}), 'b', histBins)
xlabel('Peak Amplitude')
subtitle('NREM Ripples')
legend(lgndstr)

histBins = 200;
subplot(2, 2, 4)
plotHist(ripp.dur(idx_win{1}), 'k', histBins)
hold on
plotHist(ripp.dur(idx_win{2}), 'b', histBins)
xlabel('Ripple Duration [s]')
subtitle('NREM Ripples')
legend(lgndstr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% separate ripples 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get ripples
ripp1 = getRipples('basepath', basepath, 'ch', [12 : 15],...
    'recWin', recWin(1, :), 'graphics', true, 'saveVar', false);
ripp2 = getRipples('basepath', basepath, 'ch', [12 : 15],...
    'recWin', recWin(2, :), 'graphics', true, 'saveVar', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% separate nrem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rippels in NREM
ripp1_idx = InIntervals(ripp1.peakPos, ss.stateEpochs{4});
ripp2_idx = InIntervals(ripp2.peakPos, ss.stateEpochs{4});

% NREM stats
nrem_dur = ss.stateEpochs{4}(:, 2) - ss.stateEpochs{4}(:, 1);
for iwin = 1 : length(recWin)
    nrem_idx(iwin, :) = InIntervals(ss.stateEpochs{4}, recWin(iwin, :));
end

% ripple rate
winCalc = ss.stateEpochs{4};
[ripp1.rate.rate, ripp1.rate.binedges, ripp1.rate.tstamps] =...
    times2rate(ripp1.peakPos, 'binsize', 60, 'winCalc', winCalc,...
    'c2r', true);
[ripp2.rate.rate, ripp2.rate.binedges, ripp2.rate.tstamps] =...
    times2rate(ripp2.peakPos, 'binsize', 60, 'winCalc', winCalc,...
    'c2r', true);

fh = figure;
plot(ripp1.rate.tstamps / 60 / 60, ripp1.rate.rate, 'k')
hold on
plot(ripp2.rate.tstamps / 60 / 60, ripp2.rate.rate, 'b')

mean(ripp1.rate.rate)
mean(ripp2.rate.rate)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure;
subplot(2, 2, 1)
plotHist(ripp1.dur, 'k')
hold on
plotHist(ripp2.dur, 'b')
xlabel('Ripple Duration [s]')

subplot(2, 2, 2)
plotHist(ripp1.peakAmp, 'k')
hold on
plotHist(ripp2.peakAmp, 'b')
xlabel('Peak Amp')

subplot(2, 2, 3)
plotHist(ripp1.peakFreq, 'k')
hold on
plotHist(ripp2.peakFreq, 'b')
xlabel('Peak Frequency [Hz]')

subplot(2, 2, 4)
plotHist(ripp1.maxFreq, 'k')
hold on
plotHist(ripp2.maxFreq, 'b')
xlabel('Max Frequency [Hz]')



% stats as a function of time
fh = figure;
data = movmedian(ripp.peakAmp, 993);
plot(ripp.peakPos /  60 / 60, data)




