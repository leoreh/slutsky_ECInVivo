session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'forceDef', true, 'forceL', true, 'saveVar', true);      
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
[~, basename] = fileparts(basepath);

% params
sbands = [2, 3, 4];
ch = 1;
nbins = 96;

% load
load([basename, '.spec.mat'])
load([basename, '.sleep_sig.mat'], 'emg_rms');

specData = squeeze(spec.bands.db(:, :, ch));
specNorm = specData ./ specData(1, :);

% downsample emg
npnts = length(spec.bands.db);
tstamps = [1 : npnts] / (1 / spec.info.winstep);
emgData = [interp1([1 : length(emg_rms)], emg_rms, tstamps,...
    'spline')]';

% separate low and high emg
prct = 50;
highemg = emgData > prctile(emgData, prct);
lowemg = emgData <= prctile(emgData, 100 - prct);

clear avgdb_low avgdb_high
cnt = 1;
for iband = 1 : length(sbands)
    highdb = specData(cnt, highemg);
    lowdb = specData(cnt, lowemg);

    % separate to bins 
    [~, ~, locs] = histcounts([1 : sum(highemg)], nbins);
    avgdb_high(:, cnt) = accumarray(locs(:), highdb(:)) ./ accumarray(locs(:), 1);
    [~, ~, locs] = histcounts([1 : sum(lowemg)], nbins);
    avgdb_low(:, cnt) = accumarray(locs(:), lowdb(:)) ./ accumarray(locs(:), 1);
    cnt = cnt + 1;
end

fh = figure;
th = tiledlayout(2, 1, "TileSpacing", "compact");
nexttile
plot(avgdb_low)
legend(spec.bands.bandNames{sbands})

nexttile
plot(avgdb_high)
legend(spec.bands.bandNames{sbands})

