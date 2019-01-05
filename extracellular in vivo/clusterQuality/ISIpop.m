function ISIpop(basepath, spikes)

% INPUT
%   basepath    path to recording
%   spikes      struct (see getSpikes)
%
% CALLS
%   plotWaveform
%
% 04 jan 19 18 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nargs = nargin;
if nargs < 1 || isempty(basepath)
    basepath = pwd;
end
if nargs < 2 || isempty(spikes)
    warning('spikes will be loaded from %s', basepath)
    spikes = getSpikes('basepath', basepath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc ISI contamination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nunits = length(spikes.UID);
isi = zeros(nunits, 1);
for i = 1 : length(spikes.UID)
    isi(i, 1) = sum(diff(spikes.times{i}) < 0.002) / (length(spikes.times{i}) - 1) * 100;
    isi(i, 2) = sum(diff(spikes.times{i}) < 0.003) / (length(spikes.times{i}) - 1) * 100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot hist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;

binsize = 0.05;
bins = [0 : binsize : 5];

subplot(2, 2, 1)
h = histogram(isi(:, 1), bins);
h.EdgeColor = 'none';
axis tight
ax = gca;
line([0.1 0.1], ax.YLim, 'Color', 'k', 'LineWidth', 2)
line([0.5 0.5], ax.YLim, 'Color', 'k', 'LineWidth', 2)
title('ISI within 2 ms')
ylabel('No. clusters')

subplot(2, 2, 2)
h = histogram(isi(:, 2), bins);
h.EdgeColor = 'none';
axis tight
ax = gca;
line([1 1], ax.YLim, 'Color', 'k', 'LineWidth', 2)
title('ISI within 3 ms')
ylabel('No. clusters')
xlabel('ISI')

subplot(2, 2, 3)
h = histogram(spikes.L, bins);
h.EdgeColor = 'none';
axis tight
ax = gca;
line([0.05 0.05], ax.YLim, 'Color', 'k', 'LineWidth', 2)
title('L ratio')
ylabel('No. clusters')
xlabel('L ratio')

subplot(2, 2, 4)
binsize = 0.5;
bins = [0 : binsize : 100];
h = histogram(spikes.iDist, bins);
h.EdgeColor = 'none';
axis tight
ax = gca;
line([20 20], ax.YLim, 'Color', 'k', 'LineWidth', 2)
title('Isolation Distance')
ylabel('No. clusters')
xlabel('L ratio')

savePdf('ISIpop', basepath, f)

end

% EOF


