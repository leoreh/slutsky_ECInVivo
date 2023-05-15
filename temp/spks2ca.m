

% INPUT
%   spktimes    numeric vector of spktimes [s].
%   sig         calcium. 
%   fs          sampling frequency of calcium.

% manually align spktimes and calcium
manShift = manAlign('spktimes', spktimes, 'sig', sig, 'fs', fs);

% apply shift to spktimes.
spktimes = spktimes - manShift; 

% recording duration from calcium fs
tstamps = [1 : length(sig)] / fs;
recDur = tstamps(end);

% count spikes in bins corresponding to calcium
bins = [0 : 1 / fs : spktimes(end)];
fr = histcounts(spktimes, bins, 'Normalization', 'count');

fh = figure;
plot(tstamps, sig);
yyaxis right
plot(bins(2 : end), fr)