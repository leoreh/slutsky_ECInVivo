
nEvents = length(events.peaks);
fs = lfp.fs;
marg = 0.1 * fs;

% find baseline
[~, idx1] = min(abs(events.peaks / fs - info.blockduration(1)));
idx1 = 1 : idx1 - 1;    % remove one event to be sure

% find gabazine
[~, idx2] = min(abs(events.peaks / fs - sum(info.blockduration(1:3))));
[~, idx3] = min(abs(events.peaks / fs - sum(info.blockduration(1:2))));
idx2 = idx3 + 500 : idx2 - 1;

% go over events and calc params
for i = 1 : nEvents
    
    x = events.peaks(i);
    xidx = (x - marg) : (x + marg);
    y(i, :) = double(lfp.data(xidx, 16));
    
    auc(i) = trapz(abs(y(i, :)));
    p2p(i) = abs(max(y(i, :)) - min(y(i, :)));
    % calc psd at f frequencies and normalize to total power  
    f = 1  : 250;
    p(i, :) = pwelch(y(i, :), [], [], f, fs);
    p(i, :) = p(i, :) / sum(p(i, :));

end


% plot population characteristics
nbins = 40;
figure;
subplot(1, 3, 1)
h = histogram(log10(p2p(idx1)), nbins, 'Normalization', 'probability');
h.EdgeColor = 'none';
h.FaceColor = 'k';
hold on
h = histogram(log10(p2p(idx2)), nbins, 'Normalization', 'probability');
h.EdgeColor = 'none';
h.FaceColor = 'c';
legend('Baseline', 'Gabazine')
box off
xlabel('Peak to peak [log(V)]')
ylabel('Probability')
set(gca,'TickLength',[0 0])

subplot(1, 3, 2)
h = histogram(log10(auc(idx1)), nbins, 'Normalization', 'probability');
h.EdgeColor = 'none';
h.FaceColor = 'k';
hold on
h = histogram(log10(auc(idx2)), nbins, 'Normalization', 'probability');
h.EdgeColor = 'none';
h.FaceColor = 'c';
box off
xlabel('AUC [log(V)]')
ylabel('Probability')
set(gca,'TickLength',[0 0])

pp = sum(p(:, 70 : 250), 2);
subplot(1, 3, 3)
h = histogram(log10(pp(idx1)), nbins, 'Normalization', 'probability');
h.EdgeColor = 'none';
h.FaceColor = 'k';
hold on
h = histogram(log10(pp(idx2)), nbins, 'Normalization', 'probability');
h.EdgeColor = 'none';
h.FaceColor = 'c';
box off
xlabel('PSD 50:250Hz [log(V)]')
ylabel('Probability')
set(gca,'TickLength',[0 0])

% plot individual bursts
close all
idx = [idx1(1 : 2) idx2(1 : 2)];
for i = idx    
    figure
    subplot(2, 1, 1)
    plot(xidx / fs, y(i, :))
    axis tight
    ylim([-4000 1000])
    subplot(2, 1, 2)
    pwelch(y(i, :), [], [], f, fs)
    ylim([0 60])
end

% separate bursts and iis
iis = auc > 10^4.9 & p2p > 10^3.1;
bs = find(~iis);
iis = find(iis);

t = 1;
for i = 1 : 50
    figure
    subplot(2, 1, 1)
    plot(xidx / fs, y(bs(t), :))
    axis tight
    ylim([-4000 1000])
    subplot(2, 1, 2)
    plot(xidx / fs, y(iis(t), :))
    axis tight
    ylim([-4000 1000])
    t = t + 1;
end

figure
plot(y(iis(1), :))




% pwelch
[p] = pwelch(y(i, :), [], [], [50 100], fs);

% fft
L = size(y, 2);
n = 2 ^ nextpow2(L);
ft = fft(y(i, :), n, 2);
p2 = abs(ft / L);
p1 = p2(1 : L / 2 + 1);
p1(2 : end - 1) = 2 * p1(2 : end - 1);
f = fs * (0 : (L / 2)) / L;
plot(f, p1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|p1(f)|')








% ch = 16;
% fs = 1250;
% type = 'butter';
% order = 3; 
% winLength = 11;                 % [samples]    
% 
% % filter    
% passband = [1 100];
% sigfilt = filterLFP(double(lfp.data(:, ch)), 'type', type,...
%     'passband', passband, 'order', order, 'graphics', false);
% 
% %%% Time-frequency analysis
% idx = 30*60 * fs : 35*60 * fs;
% % stft
% figure; spectrogram(sigfilt(idx), 1000, [], [0 : 1 : 30], fs, 'yaxis')
% % cwt
% figure; cwt(sigfilt(idx), 'morse', fs, 'FrequencyLimits', [1 500])
% 
% %%% single IIS analysis
% % single IIS
% dur = fs * 0.8; % burst duration in samples
% idx = 2582350 : 2582650;        % iis
% idx = 323000 : 323600;          % burst
% iis = lfp.data(idx, ch);
% figure;
% plot(idx / fs, iis)
% axis tight
% 
% % filter    
% passband = [30 500];
% x = filterLFP(double(iis), 'type', type,...
%     'passband', passband, 'order', order, 'graphics', false);
% 
% % fft
% figure
% l = length(x);
% y = fft(x);
% P2 = abs(y/l);
% P1 = P2(1:l/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(l/2))/l;
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% 
% %%% detect BS
% passband = [2 100];
% x = filterLFP(double(lfp.data(:, 16)), 'order', 3, 'type', 'butter',...
%     'passband', passband, 'graphics', false);
% 
% % x = double(lfp.data(:, 16));
% t = [1 : length(x)] / fs;
% winLength = 30;
% 
% % rectify
% y = x .^ 2;
% % moving average
% y = movmean(y, winLength);
% % standerdize
% y = (y - mean(y(:))) / std(y(:));
% 
% % y = fastrms(x);
% % y = (y - mean(y(:))) / std(y(:));
% 
% figure;
% plot(t, x)
% hold on
% yyaxis right
% p = plot(t, y);
% p.Color(4) = 0.3;
% % axis tight
% xlim([1800 2200])
% 
% idx = 25*60 * fs : 50*60 * fs;
% figure; spectrogram(x(idx), 10000, [], [0 : 10 : 200], fs, 'yaxis')
% figure; cwt(x(idx), 'morse', fs, 'FrequencyLimits', [1 500])
% 
% bs = findLFPevents('lfp', lfp, 'emgThr', 0, 'basepath', basepath, 'ch', [16], 'preset', 'bs');
