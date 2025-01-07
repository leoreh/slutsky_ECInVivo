
basepaths = {'F:\Data\lh96\lh96_220120_090157',...
    'F:\Data\lh96\lh96_220121_090213',...
    'F:\Data\lh96\lh96_220122_090154',...
    'F:\Data\lh96\lh96_220123_090009',...
    'F:\Data\lh96\lh96_220124_090127',...
    'F:\Data\lh96\lh96_220125_090041',...
    'F:\Data\lh96\lh96_220126_085016',...
    };


% load vars from each session
varsFile = ["datInfo"; "session"];
varsName = ["datInfo"; "session"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

% analyze
for isession = 1 : nsessions

    % file
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % params
    session = v(isession).session;
    nchans = session.extracellular.nChannels;
    fs = session.extracellular.sr;
    spkgrp = session.extracellular.spikeGroups.channels;

    % add timebins to datInfo
    nbins = 4;
    reqPnt = [];
    [timebins, timepnt] = metaInfo_timebins('reqPnt', reqPnt,...
        'nbins', nbins);
    winCalc = mat2cell(timebins, ones(nbins, 1), 2);

    % spk lfp
%     frange = [0.5, 2; 1, 4; 5, 11; 12, 18; 20, 35; 35, 50; 50, 70; 70, 100];
%     s = spklfp_wrapper('basepath', basepath, 'winCalc', winCalc,...
%         'ch', 9, 'frange', frange,...
%         'graphics', true, 'saveVar', false);
    
    % spike timing metrics
%     st = spktimesMetrics('winCalc', winCalc, 'forceA', true);
    
    % load lfp
    sig = double(bz_LoadBinary([basename, '.lfp'], 'duration', Inf,...
        'frequency', 1250, 'nchannels', nchans, 'start', 0,...
        'channels', 5 : 8, 'downsample', 1));
    sig = mean(sig, 2);

    % spec band
    spec = calc_spec('basepath', basepath, 'sig', sig,...
        'fs', 1250, 'graphics', true, 'saveVar', true,...
        'padfft', -1, 'winstep', 1, 'logfreq', false,...
        'ftarget', [0.4 : 0.2 : 120]);

end

% load vars from each session
varsFile = ["datInfo"; "session"; "spklfp"; "st_metrics"; "units"];
varsName = ["datInfo"; "session"; "s"; "st"; "units"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

brstVar = 'royer2';
unitType = 'rs';

nbins = 4;
cnt = 1;
clear mrld mrlt brst mrls
for isession = 1 : nsessions
    nunits = length(v(isession).s(ibin).phase.mrl(:, 1));
    
    su = v(isession).units.(unitType);
    
    for ibin = 1 : nbins
        mrls{cnt} = v(isession).s(ibin).phase.mrl(:, 1);
        mrld{cnt} = v(isession).s(ibin).phase.mrl(:, 2);
        mrlt{cnt} = v(isession).s(ibin).phase.mrl(:, 3);
        brst{cnt} = v(isession).st.(brstVar)(ibin, su);
        cnt = cnt + 1;
    end
end

mrls = cell2nanmat(mrls, 2);
mrld = cell2nanmat(mrld, 2);
mrlt = cell2nanmat(mrlt, 2);
brst = cell2nanmat(brst, 2);

xdata = [-5 * 6 : 6 : 136];

% graphics
setMatlabGraphics(true)
fh = figure;
th = tiledlayout(4, 1, 'TileSpacing', 'Compact');

nexttile
plot_boxMean('dataMat', mrls, 'clr', 'm', 'allPnts', true)
xticklabels(xdata);
ylabel('SWA [0.5-2 Hz] MRL')

nexttile
plot_boxMean('dataMat', mrld, 'clr', 'r', 'allPnts', true)
xticklabels(xdata);
ylabel('Delta [2-4 Hz] MRL')

nexttile
plot_boxMean('dataMat', mrlt, 'clr', 'b', 'allPnts', true)
xticklabels(xdata);
ylabel('Theta [5-12 Hz] MRL')

nexttile
plot_boxMean('dataMat', brst, 'clr', 'k', 'allPnts', true)
xticklabels(xdata);
xlabel('Time [h]')
ylabel('Burstiness')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short time reprasentative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = 'F:\Data\lh96\lh96_220125_090041';
cd(basepath)
[~, basename] = fileparts(basepath);
lfpfile = fullfile(basepath, [basename, '.lfp']);
datfile = fullfile(basepath, [basename, '.dat']);

% load vars 
varsFile = ["datInfo"; "session"; "spikes"; "units"];
varsName = ["datInfo"; "session"; "spikes"; "units"];
[v, ~] = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% params from session info
nchans = v.session.extracellular.nChannels;
fsLfp = v.session.extracellular.srLfp;
fs = v.session.extracellular.sr;

winCalc = [2 * 60, 2 * 60 + 30 * 60];
recDur = diff(winCalc);

% load
rawData = double(bz_LoadBinary(datfile, 'duration', recDur,...
    'frequency', fs, 'nchannels', nchans, 'start', winCalc(1),...
    'channels', 7, 'downsample', 1));

recDur = diff(winCalc);
lfpData = double(bz_LoadBinary(lfpfile, 'duration', recDur,...
    'frequency', fsLfp, 'nchannels', nchans, 'start', winCalc(1),...
    'channels', [1 : 15], 'downsample', 1));
lfpData = mean(lfpData, 2);

% filter
swa = filterLFP(lfpData, 'fs', fsLfp, 'type', 'butter', 'dataOnly', true,...
    'order', 3, 'passband', [1, 4], 'graphics', false);
theta = filterLFP(lfpData, 'fs', fsLfp, 'type', 'butter', 'dataOnly', true,...
    'order', 3, 'passband', [5, 11], 'graphics', false);


winPlot = 10 * 60;
winPlot = [winPlot, winPlot + 300];

% arrange spike times
spktimes = cellfun(@(x) [x(InIntervals(x, winPlot))]',...
    v.spikes.times, 'uni', false)';
[rs, sidx] = sort(v.units.rs, 'descend');
rs = logical(rs);
spktimes = spktimes(sidx);


% graphics
setMatlabGraphics(false)
fh = figure;
th = tiledlayout(3, 1, 'TileSpacing', 'Compact');

% raw data
nh1 = nexttile;
xWin = winPlot * fs;
xval = [xWin(1) : xWin(2)] / fs;
plot(xval, rawData(xWin(1) : xWin(2)))
ylabel('Amplitude [uV]')
yticks([])
title('Raw Data')

% sw
nh2 = nexttile;
xWin = winPlot * fsLfp;
xval = [xWin(1) : xWin(2)] / fsLfp;
plot(xval, swa(xWin(1) : xWin(2)))
ylabel('Amplitude [uV]')
yticks([])
title('Delta [1-4 Hz]')

% theta
% nh3 = nexttile;
% plot(xval, theta(xWin(1) : xWin(2)))
% yticks([])
% ylabel('Amplitude [uV]')
% title('Theta [5-11 Hz]')

% raster plot
LineFormat = struct();
LineFormat.Color = [0.1 0.1 0.8];
lineFormat.LineWidth = 3;
nh4 = nexttile;
plotSpikeRaster(spktimes(rs),...
    'PlotType', 'vertline', 'LineFormat', LineFormat);
hold on
LineFormat.Color = [0.8 0.1 0.1];
plotSpikeRaster(spktimes(~rs),...
    'PlotType', 'vertline', 'LineFormat', LineFormat,...
    'VertSpikePosition', sum(rs));
ylim([0, length(rs)])
yticks([])
ylabel('Cell Number')
title('RasterPlot]')
xlabel('Time [s]')

linkaxes([nh1, nh2, nh3, nh4], 'x')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lfp power in states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
varsFile = ["datInfo"; "session"; "spec"];
varsName = ["datInfo"; "session"; "spec"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

% analyze
for isession = 1 : nsessions

    % file
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    
    % params
    session = v(isession).session;
    nchans = session.extracellular.nChannels;
    fs = session.extracellular.sr;
    spkgrp = session.extracellular.spikeGroups.channels;
    
    % add timebins to datInfo
    timebins = v(isession).session.general.timebins;
    winCalc = mat2cell(timebins, ones(size(timebins, 1), 1), 2);
    
    % files
    sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

    % load emg
    load(sleepfile, 'emg_rms');    

    % spec bands
    frange = [0.5, 1; 1, 4; 5, 12];
    fnames = ["sw"; "delta"; "theta"];

    for iwin = 1 : size(winCalc, 1)

        % divide emg to percentiles. assumes winCalc is integers in seconds
        prct = 70;
        highemg = find(emg_rms > prctile(emg_rms, prct));
        lowemg = find(emg_rms < prctile(emg_rms, 100 - prct));
        sleepIdx = lowemg(lowemg < winCalc{iwin}(2) &...
            lowemg > winCalc{iwin}(1));
        wakeIdx = highemg(highemg < winCalc{iwin}(2) &...
            highemg > winCalc{iwin}(1));

        for ifreq = 1 : size(frange, 1)

            % index to frequency
            [~, fidx(1)] = min(abs(v(isession).spec.freq - frange(ifreq, 1)));
            [~, fidx(2)] = min(abs(v(isession).spec.freq - frange(ifreq, 2)));

            % mean power of spectrogram during wake or sleep
            specBand.(fnames{ifreq}).sleepPow{iwin} =...
                mean(v(isession).spec.s(sleepIdx, fidx(1) : fidx(2)), 2, 'omitnan');
            specBand.(fnames{ifreq}).wakePow{iwin} =...
                mean(v(isession).spec.s(wakeIdx, fidx(1) : fidx(2)), 2, 'omitnan');

            % organize in struct
            specBand.(fnames{ifreq}).fidx = fidx;            
        end
    end
    
    % save
    save([basename, '.specBand.mat'], 'specBand')
end




% concat
nbins = 4;
cnt = 1;
clear swSleep swWake deltaSleep deltaWake thetaSleep thetaWake
for isession = 1 : nsessions
    
    cd(basepaths{isession})
    [~, basename] = fileparts(basepaths{isession});
    load([basename, '.specBand.mat'], 'specBand')

    for ibin = 1 : nbins
        swSleep{cnt} = specBand.sw.sleepPow{ibin};
        swWake{cnt} = specBand.sw.wakePow{ibin};
        deltaSleep{cnt} = specBand.delta.sleepPow{ibin};
        deltaWake{cnt} = specBand.delta.wakePow{ibin};
        thetaSleep{cnt} = specBand.theta.sleepPow{ibin};
        thetaWake{cnt} = specBand.theta.wakePow{ibin};
        cnt = cnt + 1;
    end
end

% cat to mat
swSleep = cell2nanmat(swSleep, 2);
swWake = cell2nanmat(swWake, 2);
deltaSleep = cell2nanmat(deltaSleep, 2);
deltaWake = cell2nanmat(deltaWake, 2);
thetaSleep = cell2nanmat(thetaSleep, 2);
thetaWake = cell2nanmat(thetaWake, 2);


% graphics
xdata = [-5 * 6 : 6 : 136];
setMatlabGraphics(false)
fh = figure;
th = tiledlayout(2, 3, 'TileSpacing', 'Compact');

nexttile
bar(xdata, median(swSleep, 1, 'omitnan'))
hold on
errorbar(xdata, mean(swSleep, 1, 'omitnan'),...
    std(swSleep, [], 1, 'omitnan') / sqrt(length(swSleep)))

nexttile
bar(xdata, median(deltaSleep, 1, 'omitnan'))
hold on
errorbar(xdata, mean(deltaSleep, 1, 'omitnan'),...
    std(deltaSleep, [], 1, 'omitnan') / sqrt(length(deltaSleep)))

nexttile
bar(xdata, median(thetaSleep, 1, 'omitnan'))
hold on
errorbar(xdata, mean(thetaSleep, 1, 'omitnan'),...
    std(thetaSleep, [], 1, 'omitnan') / sqrt(length(thetaSleep)))

nexttile
bar(xdata, median(swWake, 1, 'omitnan'))
hold on
errorbar(xdata, mean(swWake, 1, 'omitnan'),...
    std(swWake, [], 1, 'omitnan') / sqrt(length(swWake)))

nexttile
bar(xdata, median(deltaWake, 1, 'omitnan'))
hold on
errorbar(xdata, mean(deltaWake, 1, 'omitnan'),...
    std(deltaWake, [], 1, 'omitnan') / sqrt(length(deltaWake)))

nexttile
bar(xdata, median(thetaWake, 1, 'omitnan'))
hold on
errorbar(xdata, mean(thetaWake, 1, 'omitnan'),...
    std(thetaWake, [], 1, 'omitnan') / sqrt(length(thetaWake)))


% normalized
xdata = [-5 * 6 : 6 : 136];
setMatlabGraphics(false)
fh = figure;
th = tiledlayout(2, 2, 'TileSpacing', 'Compact');

nexttile
ydata = deltaSleep / mean(mean(deltaSleep(:, 1 : 5), 1, 'omitnan'));
errorbar(xdata, mean(ydata, 1, 'omitnan'), std(ydata, [], 1, 'omitnan') ./ sum(~isnan(ydata)))

nexttile
ydata = thetaSleep / mean(mean(thetaSleep(:, 1 : 5), 1, 'omitnan'));
errorbar(xdata, mean(ydata, 1, 'omitnan'), std(ydata, [], 1, 'omitnan') ./ sum(~isnan(ydata)))

nexttile
ydata = deltaWake / mean(mean(deltaWake(:, 1 : 5), 1, 'omitnan'));
errorbar(xdata, mean(ydata, 1, 'omitnan'), std(ydata, [], 1, 'omitnan') ./ sum(~isnan(ydata)))

nexttile
ydata = thetaWake / mean(mean(thetaWake(:, 1 : 5), 1, 'omitnan'));
errorbar(xdata, mean(ydata, 1, 'omitnan'), std(ydata, [], 1, 'omitnan') ./ sum(~isnan(ydata)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% firing rate during low / high emg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
varsFile = ["datInfo"; "session"; "spikes"; "units"];
varsName = ["datInfo"; "session"; "spikes"; "units"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

% analyze
tlow = 0; thigh = 0;
clear frhigh frlow rshigh fshigh rslow fslow
for isession = 1 : nsessions

    % file
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    
    % params
    session = v(isession).session;
    nchans = session.extracellular.nChannels;
    fs = session.extracellular.sr;
    spkgrp = session.extracellular.spikeGroups.channels;
    
    % files
    sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

    % load emg
    load(sleepfile, 'emg_rms');    
    
    % divide emg to percentiles and create bouts
    prct = 50;
    highemg = binary2bouts('vec', [emg_rms > prctile(emg_rms, prct)]',...
        'minDur', 2, 'maxDur', [], 'interDur', 2, 'exclude', false);
    lowemg = binary2bouts('vec', [emg_rms < prctile(emg_rms, 100 - prct)]',...
        'minDur', 2, 'maxDur', [], 'interDur', 2, 'exclude', false);
    
    % clip spktimes to bouts
    spktimes = v(isession).spikes.times;
    spkhigh = cellfun(@(x) x(InIntervals(x, highemg)), spktimes, 'uni', false);
    spklow = cellfun(@(x) x(InIntervals(x, lowemg)), spktimes, 'uni', false);

    % calc firing rate
    [frhigh, ~, tstamps] = times2rate(spkhigh, 'binsize', 60,...
        'winCalc', [0 Inf], 'c2r', true);
    thigh = [thigh, tstamps + max(thigh)];
    
    [frlow, ~, tstamps] = times2rate(spklow, 'binsize', 60,...
        'winCalc', [0 Inf], 'c2r', true);
    tlow = [tlow, tstamps + max(tlow)];
    
    % separate rs and fs
    rshigh{isession} = frhigh(v(isession).units.rs, :);
    fshigh{isession} = frhigh(v(isession).units.fs, :);
    rslow{isession} = frlow(v(isession).units.rs, :);
    fslow{isession} = frlow(v(isession).units.fs, :);

end

% correct tstamps
thigh = thigh(2 : end);
tlow = tlow(2 : end);

% cat
rshighmat = cell2nanmat(rshigh);
fshighmat = cell2nanmat(fshigh);
rslowmat = cell2nanmat(rslow);
fslowmat = cell2nanmat(fslow);


% graphics
fh = figure;
th = tiledlayout(2, 1, 'TileSpacing', 'Compact');

nexttile
plot(thigh / 60 / 60, mean(rshighmat, 1, 'omitnan'))
hold on
plot(tlow / 60 / 60, mean(rslowmat, 1, 'omitnan'))
title('rs')

nexttile
plot(thigh / 60 / 60, mean(fshighmat, 1, 'omitnan'))
hold on
plot(tlow / 60 / 60, mean(fslowmat, 1, 'omitnan'))
title('fs')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% firing rate during low / high emg - alt 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
varsFile = ["datInfo"; "session"; "spikes"; "units"; "fr"];
varsName = ["datInfo"; "session"; "spikes"; "units"; "fr"];
[v, ~] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

% analyze
tlow = 0; thigh = 0;
clear frhigh frlow rshigh fshigh rslow fslow
for isession = 1 : nsessions

    % file
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    
    % params
    session = v(isession).session;
    nchans = session.extracellular.nChannels;
    fs = session.extracellular.sr;
    spkgrp = session.extracellular.spikeGroups.channels;
    
    % files
    sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

    % load emg
    load(sleepfile, 'emg_rms');    
    
    % calculate mean emg rms in 1 min bins
    binsize = v(isession).fr.info.binsize;
    nsec = length(emg_rms);
    m = nsec - mod(nsec, binsize);
    emg_mat = reshape(emg_rms(1 : m), binsize, []);
    emg_bins = mean(emg_mat, 1, 'omitnan');
    
    % divide emg to percentiles and create bouts
    prct = 50;
    highemg = emg_bins > prctile(emg_bins, prct);
    lowemg = emg_bins < prctile(emg_bins, 100 - prct);   

    % separate rs and fs
    rs = v(isession).units.rs;
    fs = v(isession).units.fs;

    % fr  
    rshigh{isession} = v(isession).fr.strd(rs, highemg);
    fshigh{isession} = v(isession).fr.strd(fs, highemg);
    rslow{isession} = v(isession).fr.strd(rs, lowemg);
    fslow{isession} = v(isession).fr.strd(fs, lowemg);

    % timestamps
    thigh = [thigh, v(isession).fr.tstamps(highemg) + max(thigh)];
    tlow = [tlow, v(isession).fr.tstamps(lowemg) + max(tlow)];

end

% correct tstamps
thigh = thigh(2 : end) / 60 / 60 - 32;
tlow = tlow(2 : end) / 60 / 60 - 32;

% cat
rshighmat = cell2nanmat(rshigh);
fshighmat = cell2nanmat(fshigh);
rslowmat = cell2nanmat(rslow);
fslowmat = cell2nanmat(fslow);


% graphics
fh = figure;
th = tiledlayout(2, 1, 'TileSpacing', 'Compact');

nexttile
plot(thigh, mean(rshighmat, 1, 'omitnan'))
hold on
plot(tlow, mean(rslowmat, 1, 'omitnan'))
title('rs')

nexttile
plot(thigh, mean(fshighmat, 1, 'omitnan'))
hold on
plot(tlow, mean(fslowmat, 1, 'omitnan'))
title('fs')