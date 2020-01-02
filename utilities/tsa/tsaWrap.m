
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% args
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('freerun[9-9].mat')
t = T9;
orig = Y9';
dt = mean(diff(t));
fs = 1 / dt;
l = t(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect PLI crossings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bandpass = [45 65];
LINE = lineDetect(orig, fs, bandpass);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove interference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tsa = [];
DC = [];
MA = 0;
TW = 1;
[clean tsa wvec xmat] = lineRemove(orig, LINE, tsa, DC, MA, TW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2, 2, 1)
plot(t, orig), hold on, plot(t, clean)
separators(LINE / fs, [], [0 0.7 0 ]);
legend('original', 'clean', 'line crossings')
xlabel('time')
ylabel('amplitude')
title('signal')
axis tight
xlim([25 25.2])
ylim([-0.2 0.2])

n = 2^nextpow2(l);
f = fs * (0 : (n / 2)) / n;
pclean = abs(fft(clean, n) / n);
porig = abs(fft(orig, n) / n);
subplot(2, 2, 2)
plot(f, porig(1:n/2+1))
hold on
plot(f, pclean(1:n/2+1))
legend('original', 'clean')
title('power spectrum')

subplot(2, 2, 3)
plot(tsa)
xlabel('sample')
ylabel('amplitude')
title('average interference')
axis tight

subplot(2, 2, 4)
hist(diff(LINE), min(diff(LINE)) : max(diff(LINE)))
xlabel('length [samples]')
ylabel('counts')
title('cycle duration')