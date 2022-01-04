
[sig, tsaSig, revIdx] = tsa_filter('sig', EMG, 'fs', 1250, 'tw', false,...
    'graphics', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh = figure;

% example of the clean and raw siganl
subplot(2, 2, 1)
tstart = round(length(sig) / 2);     
nRev = 5;       
sigIdx = [tstart : nRev * length(tsaSig) + tstart];          
tstamps = sigIdx / fs;
lidx = find(revIdx > sigIdx(1) & revIdx < sigIdx(end));
plot(tstamps, EMG(sigIdx), 'k')
hold on
plot(tstamps, sig(sigIdx), 'b')
line([revIdx(lidx) revIdx(lidx)] / fs, ylim, 'Color', 'k', 'HandleVisibility','off')
legend('raw', 'filtered')
box off
axis tight
set(gca, 'TickLength', [0 0])
xlabel('Time [s]')

% pwelch 
subplot(2, 2, 2)
win = hann(2 ^ (nextpow2(2 * fs) - 1));
noverlap = floor(0.25 * fs);
faxis = [0.2 : 0.2 : 220];       
sigIdx = [1 : 30 * 60 * fs];    % limit duration
[pow_orig, ~] = pwelch(EMG(sigIdx), win, noverlap, faxis, fs);
[pow_clean, ~] = pwelch(sig(sigIdx), win, noverlap, faxis, fs);
ph = plot(faxis, pow_orig ./ sum(pow_orig, 2), 'k', 'LineWidth', 2);
hold on
ph = plot(faxis, pow_clean ./ sum(pow_clean, 2), 'b', 'LineWidth', 2);
legend('raw', 'filtered')
xlabel('Frequency [Hz]')
ylabel('Norm PSD')     

% tsa
subplot(2, 2, 3)
plot(tsaSig)
xlabel('sample')
ylabel('amplitude')
title('average interference')
axis tight

% histogram of cycle duration
subplot(2, 2, 4)
h = histogram(diff(revIdx / fs * 1000), 100);
xlabel('cycle duration [ms]')
ylabel('counts')
h.EdgeColor = 'none';
h.FaceColor = 'k';
box off
axis tight
set(gca, 'TickLength', [0 0])
