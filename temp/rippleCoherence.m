

% coherence 
% https://www.mathworks.com/help/signal/ug/cross-spectrum-and-magnitude-squared-coherence.html


basepath = 'F:\Data\lh123\lh123_221212_091133';
cd(basepath)
[~, basename] = fileparts(basepath);

session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'forceDef', true, 'forceL', true, 'saveVar', true);      
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
[~, basename] = fileparts(basepath);


ch = [spkgrp{:}];

lfpfile = fullfile(basepath, [basename, '.lfp']);
precision = 'int16';
nbytes = class2bytes(precision);


info = dir(lfpfile);
nsamps = info.bytes / nbytes / nchans;
m = memmapfile(lfpfile, 'Format', {precision, [nchans, nsamps] 'mapped'});
raw = m.Data;


epochs = round(ripp.epochs * fs);
nepochs = size(epochs, 1)
stamps = epochs(:, 1);
win = [zeros(1, nepochs); diff(epochs')]';

for iepoch = 10 : 20

[rippData, ~] = snipFromBinary('stamps', stamps(iepoch), 'fname', lfpfile,...
    'win', win(iepoch, :), 'nchans', nchans, 'ch', ch, 'align_peak', 'no',...
    'precision', precision, 'rmv_trend', [], 'saveVar', false,...
    'l2norm', false, 'raw', raw);

[cxy, f] = mscohere(rippData(2, :)', rippData(6, :)', [], [], [], 1250);

fh = figure;
subplot(1, 2, 1)
plot(f, cxy)
title('Magnitude-Squared Coherence')
xlabel('Frequency (Hz)')
grid

pxy = cpsd(rippData(2, :), rippData(6, :), [], [], [], 1250);

subplot(1, 2, 2)
plot(f, angle(pxy) / pi)
title('Cross Spectrum Phase')
xlabel('Frequency (Hz)')
ylabel('Lag (\times\pi rad)')
grid

end




