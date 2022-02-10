
for isession = 1 : nsessions
    
    % load
    % basepath = pwd;
    basepath = basepaths{isession};
    [~, basename] = fileparts(basepath);
    cd(basepath)
    
    load([basename, '.spktimes.mat'])
    load([basename, '.spikes.cellinfo.mat'])
    load([basename, '.units.mat'])
    load([basename, '.cell_metrics.cellinfo.mat'])
    
    fsLfp = 1250;
    
    % clip spikes times to reasonable window
    winCalc = [0, 180 * 60];
    stimes = cellfun(@(x) x(InIntervals(x, winCalc))', spikes.times, 'uni', false);
    
    % sort by fr
    [~, sidx] = sort(cellfun(@length, stimes, 'uni', true), 'descend');
    stimes = stimes(sidx);
    
    % choose pyr
    stimes = stimes(units.idx(1, :));
    
    % load lfp
    recDur = diff(winCalc);
    sig = double(bz_LoadBinary([basename, '.lfp'], 'duration', recDur,...
        'frequency', fsLfp, 'nchannels', 19, 'start', winCalc(1),...
        'channels', [9 : 11], 'downsample', 1));
    sig = mean(sig, 2);
    
    % filter
    % sig_filt = filterLFP(sig, 'fs', fsLfp, 'type', 'butter', 'dataOnly', true,...
    %     'order', 3, 'passband', [0.5 8], 'graphics', false);
    
    % buzcode
    subpops = cell(1, length(units.rs));
    for iunit = 1 : length(subpops)
        if units.fs(iunit)
            subpops{iunit} = 'fs';
        elseif units.rs(iunit)
            subpops{iunit} = 'rs';
        end
    end
    
    lfp.data = sig;
    lfp.timestamps = [winCalc(1) : 1 / fsLfp : winCalc(2) - 1 / fsLfp];
    lfp.samplingRate = fsLfp;
    lfp.channels = 1;
    [SpikeLFPCoupling] = bz_GenSpikeLFPCoupling(spikes, lfp,...
        'spikeLim', 500000, 'frange', [1 100], 'nfreqs', [20],...
        'cellclass', subpops, 'sorttype', 'rate', 'saveMat', true);
    
    
end



varsFile = ["spkLFPcoupling"; "fr"; "units"; "datInfo"; "session"];
varsName = ["SpikeLFPCoupling"; "fr"; "units"; "datInfo"; "session"];

mname = 'lh96';
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""]);
[~, basenames] = cellfun(@fileparts, basepaths, 'uni', false);
nsessions = length(basepaths);

xLabels = ["Baseline"; "Baseline"; "Baclofen"; "Baclofen"; "Baclofen";...
    "Baclofen"; "Washout"];

fh = figure;
freqs = v(isession).SpikeLFPCoupling.freqs;
for isession = 1 : nsessions
    sunits = v(isession).units.idx(1, :)' & v(isession).fr.mfr > 0.8;
    data = v(isession).SpikeLFPCoupling.cell.spikephasemag(sunits, :);
    subplot(2, 4, isession)
    plot(freqs, mean(data, 'omitnan'), 'k')
    hold on
    patch([freqs, flip(freqs)], [mean(data, 'omitnan') + std(data, 'omitnan'),...
        flip(mean(data, 'omitnan') - std(data, 'omitnan'))],...
        'k', 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    title(xLabels{isession})
    xlabel('Frequency [Hz]')
    ylabel('Spike-Phase Magnitude')
    set(gca, 'xscale', 'log')
    ylim([0 0.2])
end



% get phase
h = hilbert(sig_filt);
sig_phase = angle(h);
sig_amp = abs(h);

windowSize = [0.5];
pflag = 0;
clear ccurve
for iunit = 1 : 5
    [sta, ccurve(iunit, :), ~, ~, ~, ~, freq] =...
        getCoherence(stimes{iunit} - winCalc(1), sig, fsLfp, windowSize,...
        0, pflag, 0.5, 50);
end

fh = figure;
plot(freq, mean(ccurve, 1, 'omitnan'))
xlim([0 50])
ylim([0 0.5])
title(basename)






% plot raster
fh = figure;
sb1 = subplot(3, 1, 1);
plot([winCalc(1) : 1 / fsLfp : winCalc(2) - 1 / fsLfp],...
    sig, 'LineWidth', 1, 'Color', 'k')
title(basename)

sb2 = subplot(3, 1, 2);
plot([winCalc(1) : 1 / fsLfp : winCalc(2) - 1 / fsLfp],...
    sig_filt, 'LineWidth', 1, 'Color', 'k')
title(basename)

sb3 = subplot(3, 1, 3);
plotSpikeRaster(stimes(rs), 'PlotType', 'vertline', 'VertSpikeHeight', 1,...
    'XLimForCell', [winCalc(1), winCalc(1) + 30]);

linkaxes([sb1, sb2, sb3], 'x')
xlim([winCalc(1), winCalc(1) + 5])






% cat mu across tetrodes
stimes = sort(vertcat(spktimes{:}) / fs);

% recalculate firing rate in 1 s bins with a moving window of 25 ms
hbins = [0 : 1 : max(stimes)];
hcnts = histcounts(stimes, hbins);

fh = figure;
plot(hbins(2 : end) / 2 / 60 / 60, hcnts)