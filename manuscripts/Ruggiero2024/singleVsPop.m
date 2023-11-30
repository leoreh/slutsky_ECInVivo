
% -------------------------------------------------------------------------
% estimated variance of delta MFR


% Sample data for MFR before and after perturbation
x = mfr(:, 1);
y = mfr(:, 2);

% Calculate delta MFR
dMFR = y - x;

% Calculate the variance of delta MFR using the formula:
covXY = cov(y, x);
var_delta_MFR = var(y) + var(x) - 2 * covXY(2, 1);

disp(['Expected: ', num2str(var_delta_MFR)]);
disp(['Calculated: ', num2str(var(dMFR))]);

% -------------------------------------------------------------------------
% distribution of delta MFR. comparison of baclofen and baseline

d1 = (mfr(:, 1) - mfr(:, 2)) ./ mfr(:, 1);
d2 = (mfr(:, 1) - mfr(:, 2)) ./ mfr(:, 1);


% distribution of fr changes per unit
fh = figure;
data = d1;
[pdf, x] = ksdensity(data, 'Bandwidth', binWidth); 
% smthCnts = pdf * numel(data) * binWidth; 
plot(x, pdf, 'k')
hold on
xLimit = xlim
xlabel('delta MFR')
ylabel('No. deltas')

hold on
data = d2;
[pdf, x] = ksdensity(data, 'Bandwidth', binWidth); 
% smthCnts = pdf * numel(data) * binWidth; 
plot(x, pdf, 'b')
hold on
xLimit = xlim
xlabel('delta MFR')
ylabel('No. deltas')

legend({'Baseline'; 'Baclofen'})

std(d1)
std(d2(randperm(413, 184)))


% distribution of fr changes per unit
fh = figure;
nbins = 20;
data = d1;
histogram(data, nbins, 'Normalization', 'pdf')

hold on
data = d2;
histogram(data, nbins, 'Normalization', 'pdf')

xlabel('delta MFR')
ylabel('No. deltas')

legend({'Baseline'; 'Baclofen'})


% -------------------------------------------------------------------------
% delta of simulated signal (and shuffled)

% params
trend = [0 0];
f = 100; % sine frequency
N = 1000; % number of data points
windowWidth = 100; % width of window

% create signal
t = linspace(0, 10 * pi, N); % time vector
linearTrend = linspace(trend(1), trend(2), N); % linear trend
sinWave = sin(f * t); % sine wave
sig = linearTrend + sinWave; % combined signal

% Calculate averages at the beginning and end of the signal
avgStart = mean(sig(1 : windowWidth));
avgEnd = mean(sig(end - windowWidth + 1 : end));

% Reshuffle the signal
shfSig = sig(randperm(N));

% Calculate averages of reshuffled signal
avgStartShf = mean(shfSig(1 : windowWidth));
avgEndShf = mean(shfSig(end - windowWidth + 1 : end));

% Plot original and reshuffled signals with averages
close all
figure;
subplot(2,1,1);
plot(t, sig, 'b', 'LineWidth', 2);
hold on;
plot([t(1), t(windowWidth)], [avgStart, avgStart], 'r', 'LineWidth', 2);
plot([t(end-windowWidth+1), t(end)], [avgEnd, avgEnd], 'r', 'LineWidth', 2);
title('Original Signal');
xlabel('Time');
ylabel('Amplitude');
legend('Signal', 'Averages');

subplot(2,1,2);
plot(t, shfSig, 'g', 'LineWidth', 2);
hold on;
plot([t(1), t(windowWidth)], [avgStartShf, avgStartShf], 'm', 'LineWidth', 2);
plot([t(end-windowWidth+1), t(end)], [avgEndShf, avgEndShf], 'm', 'LineWidth', 2);
title('Reshuffled Signal');
xlabel('Time');
ylabel('Amplitude');
legend('Signal', 'Averages');

% Display the results
disp(['Original Signal: Avg at start = ', num2str(avgStart), ', Avg at end = ', num2str(avgEnd)]);
disp(['Reshuffled Signal: Avg at start = ', num2str(avgStartShf), ', Avg at end = ', num2str(avgEndShf)]);











% -------------------------------------------------------------------------
% mfr per unit vs. pMFR. 
% take data from Lee excel, bsl vs. 24 hr bac

% % calc
% d = (mfr(:, 1) - mfr(:, 2));
% nunits = size(mfr, 1);
% nrep = 1000;
% d2 = nan(nunits, nrep);
% for iunit = 1 : nunits
%     for irep = 1 : nrep
%         
%         unitidx = randperm(nunits, iunit);
%         d2(iunit, irep) = mean(mfr(unitidx, 1)) - mean(mfr(unitidx, 2));
% 
% 
%     end
% end
% 
% % fig params
% nbins = 15;
% [histCnts, histCents] = histcounts(mfr(1, :), nbins);
% histCents = histCents(1 : end - 1) + (diff(histCents) / 2)
% binWidth = mean(diff(histCents));
% binWidth = 0.25
% 
% 
% fh = figure;
% 
% % fr distributions before and after compensation
% subplot(2, 2, 1)
% data = log10(mfr(:, 1));
% [pdf, x] = ksdensity(data, 'Bandwidth', binWidth); 
% smthCnts = pdf * numel(data) * binWidth; 
% plot(x, smthCnts, 'k')
% 
% hold on
% data = log10(mfr(:, 2));
% [pdf, x] = ksdensity(data, 'Bandwidth', binWidth); 
% smthCnts = pdf * numel(data) * binWidth; 
% plot(x, smthCnts, 'r')
% xlabel('log firing rate')
% ylabel('No. Units')
% 
% % distribution of fr changes per unit
% subplot(2, 2, 2)
% data = d;
% [pdf, x] = ksdensity(data, 'Bandwidth', binWidth); 
% smthCnts = pdf * numel(data) * binWidth; 
% plot(x, smthCnts, 'b')
% hold on
% xLimit = xlim
% xlabel('delta MFR')
% ylabel('No. deltas')
% 
% % distribution of delta between fr distsribtions given xx units
% subplot(2, 2, 3)
% nunits = [3, 30, 300];
% hold on
% for iunit = nunits
% data = d2(iunit, :);
% [pdf, x] = ksdensity(data, 'Bandwidth', binWidth); 
% smthCnts = pdf * numel(data) * binWidth; 
% plot(x, smthCnts)
% end
% xlim(xLimit)
% xlabel('delta pMFR')
% ylabel('No. deltas')
% nunits_str = arrayfun(@(x) ['n=' num2str(x)], nunits, 'UniformOutput', false);
% legend(nunits_str)
% 
% % std of distribution deltas as function of nunits 
% subplot(2, 2, 4)
% data = std(d2, [], 2);
% plot(data)
% ylabel('STD of delta pMFR')
% xlabel('Population size')
% 
% %%%
% fh = figure;
% % scatter(log10(mfr(:, 1)), d)
% 
% nbins = 50;
% nunits = size(mfr, 1);
% histogram((mfr(:, 2)), nbins)
% hold on
% histogram(((mfr(:, 1) + d(randperm(nunits, nunits)))), nbins)
% set(gca,'xscale','log')
% 
% 
% 
% 
% % -------------------------------------------------------------------------
% % mfr per baclofen experiment
% % take data mea folder
% 
% basepath = 'F:\Data\MEA\baclofen\210527_045600';
% cd(basepath)
% [~, basename] = fileparts(basepath);
% load([basename, '.fr.mat'])
% idx1 = [10 : 70];
% idx2 = [510 : 570];
% xfr{1, 1} = mean(fr.strd(:, idx1), 2);
% xfr{1, 2} = mean(fr.strd(:, idx2), 2);
% 
% 
% basepath = 'F:\Data\MEA\baclofen\210325_082300';
% cd(basepath)
% [~, basename] = fileparts(basepath);
% load([basename, '.fr.mat'])
% % fh = figure;
% % plot(mean(fr.strd))
% idx1 = [10 : 70];
% idx2 = [510 : 570];
% xfr{2, 1} = mean(fr.strd(:, idx1), 2);
% xfr{2, 2} = mean(fr.strd(:, idx2), 2);
% 
% basepath = 'F:\Data\MEA\baclofen\190910_045200';
% cd(basepath)
% [~, basename] = fileparts(basepath);
% load([basename, '.fr.mat'])
% % fh = figure;
% % plot(mean(fr.strd))
% idx1 = [10 : 70];
% idx2 = [480 : 540];
% xfr{3, 1} = mean(fr.strd(:, idx1), 2);
% xfr{3, 2} = mean(fr.strd(:, idx2), 2);
% 
% basepath = 'F:\Data\MEA\baclofen\190813_022600';
% cd(basepath)
% [~, basename] = fileparts(basepath);
% load([basename, '.fr.mat'])
% % fh = figure;
% % plot(mean(fr.strd))
% idx1 = [10 : 70];
% idx2 = [480 : 540];
% xfr{4, 1} = mean(fr.strd(:, idx1), 2);
% xfr{4, 2} = mean(fr.strd(:, idx2), 2);
% 
% 
% mfr = cellfun(@(x) mean(x), xfr, 'uni', true);
% 
% % calc
% d = (mfr(:, 1) - mfr(:, 2)) ./ max(mfr')';
% nunits = size(mfr, 1);
% nrep = 1000;
% d2 = nan(nunits, nrep);
% for iunit = 1 : nunits
%     for irep = 1 : nrep
%         
%         unitidx = randperm(nunits, iunit);
%         d2(iunit, irep) = mean(mfr(unitidx, 1)) - mean(mfr(unitidx, 2));
% 
% 
%     end
% end
% 
% % fig params
% nbins = 15;
% [histCnts, histCents] = histcounts(mfr(1, :), nbins);
% histCents = histCents(1 : end - 1) + (diff(histCents) / 2)
% binWidth = mean(diff(histCents));
% binWidth = 0.25
% 
% 
% fh = figure;
% 
% % fr distributions before and after compensation
% subplot(2, 2, 1)
% data = log10(mfr(:, 1));
% [pdf, x] = ksdensity(data, 'Bandwidth', binWidth); 
% smthCnts = pdf * numel(data) * binWidth; 
% plot(x, smthCnts, 'k')
% 
% hold on
% data = log10(mfr(:, 2));
% [pdf, x] = ksdensity(data, 'Bandwidth', binWidth); 
% smthCnts = pdf * numel(data) * binWidth; 
% plot(x, smthCnts, 'r')
% xlabel('log firing rate')
% ylabel('No. Units')
% 
% % distribution of fr changes per unit
% subplot(2, 2, 2)
% data = d;
% [pdf, x] = ksdensity(data, 'Bandwidth', binWidth); 
% smthCnts = pdf * numel(data) * binWidth; 
% plot(x, smthCnts, 'b')
% hold on
% xLimit = xlim
% xlabel('delta MFR')
% ylabel('No. deltas')
% 
% % distribution of delta between fr distsribtions given xx units
% subplot(2, 2, 3)
% nunits = [1, 2, 4];
% hold on
% for iunit = nunits
%     data = d2(iunit, :);
%     [pdf, x] = ksdensity(data, 'Bandwidth', binWidth);
%     smthCnts = pdf * numel(data) * binWidth;
%     plot(x, smthCnts)
% end
% xlim(xLimit)
% xlabel('delta pMFR')
% ylabel('No. deltas')
% nunits_str = arrayfun(@(x) ['n=' num2str(x)], nunits, 'UniformOutput', false);
% legend(nunits_str)
% 
% % std of distribution deltas as function of nunits 
% subplot(2, 2, 4)
% data = std(d2, [], 2);
% plot(data)
% ylabel('STD of delta pMFR')
% xlabel('Population size')
% 
