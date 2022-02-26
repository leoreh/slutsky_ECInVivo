function [spklfp] = spkLfpCoupling(varargin)

% investigate spk - lfp coupling as point to fields.

% should calculate across frequencies and for specific bands

% current metrices include (1) phase coupling per cell, including significance
% test (Rayleigh or jittering), (2) rate map for each cell by lfp phase and
% power of a single band, (3) correlation between spk counts and lfp
% magnitude across lfp frequencies, and correlation between population
% synchrony (by spk counts) and magnitude across frequencies.

% currently i think only (3) is suitable for wide band investigation.
% should think how and when to separate rs vs fs cells

% based on bz_GenSpikeLFPCoupling, bz_PhaseModulation, and
% bz_PowerPhaseRatemap. main differences from bz_GenSpikeLFPCoupling: (1)
% original also had a subrouting for calculating mutual information between
% power and ISI distribution which i ignored. could be relavent for
% theta-driving interneurons (Li, Sci. Reports, 2017). (2) original loopes
% across lfp channels, probably for linear probes.

% alternatives: (1) coherencypt from chronux which calculates the
% cross-spectrom of the multitaper specs (point-process) spktimes and
% (contineous) lfp. 
% see http://www-users.med.cornell.edu/~jdvicto/pdfs/pubo08.pdf
% (2) the sta method of Vinck 2012 which is implemented in fieldtrip. see 
% https://www.fieldtriptoolbox.org/tutorial/spikefield/
% however, i think both of these methods are mainly suitable for tasks with
% repeated  trials

% currently only a single lfp channel is used. need to think how to improve
% this. in fieldtrip they mention bleeding of a spike waveform energy into
% the lfp recorded on the same channel. could be problamatic 

% INPUT
%   basepath    char. fullpath to recording folder {pwd}
%   winCalc     n x 2 mat
%   ch          numeric vec of lfp channels to load. if vec than lfp will
%               be averaged across channels
%   frange      
%   nfreq       
%   jitterSig   logical {false} 
%   saveVar     logical {true}
%   graphics    logical {true}
%   srtUnits    logical {true}. currently not implemented
%
% CALLS
%   bz_LoadBinary
%   bz_SpktToSpkmat
%   bz_Counter
%   CircularDistribution (fmat)
%   NormToInt
%
% TO DO LIST
%       # allow user to input data (spktimes, lfp, units)
%
% 25 feb 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'winCalc', [0 Inf], @isnumeric)
addParameter(p, 'ch', 1, @isnumeric)
addParameter(p, 'frange', [1.5 100], @isnumeric)
addParameter(p, 'nfreq', 40, @isnumeric)
addParameter(p, 'jitterSig', false, @isnumeric)
addParameter(p, 'srtUnits', true, @islogical)
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
winCalc         = p.Results.winCalc;
ch              = p.Results.ch;
frange          = p.Results.frange;
nfreq           = p.Results.nfreq;
jitterSig       = p.Results.jitterSig;
srtUnits        = p.Results.srtUnits;
graphics        = p.Results.graphics;
saveVar         = p.Results.saveVar;

% direct run
% winCalc = [0, 180 * 60];
% basepath = pwd;
% ch = 9 : 11;
% frange = [10 20];
% nfreq = 1;
% jittersig
% saveVar
% graphics
% srtUnits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
[~, basename] = fileparts(basepath);
cd(basepath)
spksfile = fullfile(basepath, [basename, '.spikes.cellinfo.mat']);
unitsfile = fullfile(basepath, [basename, '.units.mat']);
lfpfile = fullfile(basepath, [basename, '.lfp']);
sessionfile = fullfile(basepath, [basename, '.session.mat']);
spklfpfile = fullfile(basepath, [basename, '.spklfp.mat']);

% load session info
if ~exist(sessionfile, 'file')
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'forceDef', false, 'forceL', false, 'saveVar', false);
else
    load(sessionfile)
end

% params from session info
nchans = session.extracellular.nChannels;
fs = session.extracellular.srLfp;

% load spikes
if exist(spksfile, 'file')
    load(spksfile, 'spikes')
end
nunits = length(spikes.times);

% get subpopulations of cells (e.g. rs and fs)
pop.subpops = cell(1, nunits);
pop.subpops(:) = {'unwn'};
if exist(unitsfile, 'file')
    load(unitsfile, 'units')
    for iunit = 1 : nunits
        if units.fs(iunit)
            pop.subpops{iunit} = 'fs';
        elseif units.rs(iunit)
            pop.subpops{iunit} = 'rs';
        end
    end
end
pop.names = unique(pop.subpops);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange lfp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load lfp
sig = double(bz_LoadBinary(lfpfile, 'duration', diff(winCalc),...
    'frequency', fs, 'nchannels', nchans, 'start', winCalc(1),...
    'channels', ch, 'downsample', 1));
sig = mean(sig, 2);
tstamps = [winCalc(1) : 1 / fs : winCalc(2) - 1 / fs];

% filter, norm, and hilbert. for multiple frequencies, bz uses the wavelet
% transform i assume to allow log-spaced frequencies and to catch obrupt
% changes. still, perhaps consider switching to multitaper spec
switch nfreq
    case 1
        sig_filt = filterLFP(sig, 'fs', fs, 'type', 'fir1', 'dataOnly', true,...
            'order', 7, 'passband', frange, 'graphics', false);
        sig_z = hilbert(sig_filt);
        sig_z = sig_z ./ mean(abs(sig_z));
        
    otherwise
        tmp = bz_WaveSpec(sig, 'showprogress', true,...
            'ncyc', 7, 'nfreqs', nfreq, 'frange', frange,...
            'samplingRate', 1250);
        freqs = tmp.freqs;
        sig_z = tmp.data ./ mean(abs(tmp.data), 1, 'omitnan');
        clear tmp
end

% must remove copies of the signal due to memory (raw, filt, hilbert)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restrict anaylsis to epochs with enough power in band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if working on multiple frequencies simultaneously, must run this in a
% loop

% approximate power in frequency band
power = fastrms(sig_filt, ceil(fs ./ frange(1)), 1);
powThr = 2;
minpow = mean(power) + std(power) * powThr;

% set the minimum duration to two cycles
mindur = (fs ./ frange(2)) * 2; 
    
% find epochs with power > low threshold. correct for durations. NOT SURE
% THIS IS THE CORRECT APPROACH FOR BAC EXPERIMENTS WHEN DELTA IS ALWAYS
% HIGH
bad_epochs = binary2epochs('vec', power < minpow, 'minDur', mindur,...
    'maxDur', Inf, 'interDur', 0);
nepochs = size(bad_epochs, 1);

% subtract out low power intervals
lfp_epochs = SubtractIntervals(winCalc, bad_epochs / fs);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% restrict spktimes to lfp epochs
spktimes = cellfun(@(x) x(InIntervals(x, lfp_epochs)),...
    spikes.times, 'uni', false);



% if computation takes too long, can add a limit to the number of spikes
% here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bz_PhaseModulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbins = 180;

% initialize
spkphase        = cell(1,length(spikes.times));
phase.dist      = nan(nbins, nunits);
phase.kappa     = nan(1, nunits);
phase.theta     = nan(1, nunits);
phase.r         = nan(1, nunits);
phase.p         = nan(1, nunits);
h = [];
for iunit = 1 : nunits
    
    if isempty(spktimes{iunit})
        continue
    end
    
    % find phase of lfp for each spike. notes: (1) in bz_PhaseModulation
    % this is done by indexing angle(sig_z) with spktimes * fs which is
    % faster but less accurate. (2) converting angles to 0 : 2pi rather
    % than -pi : pi is neccassary because that is how cicrcularDistribution
    % divides nbins. (3) in bz_genSpikeLfpCoupling they take
    % the angle and magnitude of the mean(hilbert) instead of averaging the
    % angle and magnitude separately which i think is a mistake
    spkphase{iunit} = interp1(tstamps, mod(angle(sig_z), 2 * pi), spktimes{iunit}, 'nearest');
    
    % use zugaro for circular statstics. mean phase - theta; mean resultant
    % length - r; Von Mises concentraion - kappa; Rayleigh significance of
    % uniformity - p; see also:
    % https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Analysis.pdf
    [phase.dist(:, iunit), phase.bins, tmp] = CircularDistribution(spkphase{iunit}, 'nBins', nbins);
    phase.kappa(iunit) = tmp.k;
    phase.theta(iunit) = tmp.m;
    phase.r(iunit) = tmp.r;
    phase.p(iunit) = tmp.p;

    % plotting   
    h(end+1) = figure;
    hax = subplot(1,2,1);
    rose(spkphase{iunit})
    title(sprintf('Cell %d; Rayleigh = %.2f', iunit, phase.p(iunit)))
    
    hax = subplot(1,2,2);
    bar(phase.bins * 180 / pi, phase.dist(:, iunit))
    xlim([0 360])
    set(hax,'XTick',[0 90 180 270 360])
    hold on;
    plot([0:360],cos(pi/180*[0:360])*0.05*max(phase.dist(:, iunit))+0.95*max(phase.dist(:, iunit)),'color',[.7 .7 .7])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rate map as a function of power and phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here, a bivariate histogram is created for each cell by dividing the lfp
% power and phase to nbins x nbins and counting the number of spikes in each
% bin. the counts or normalized to rate by occupance - the duration of each
% bin. based on bz_PowerPhaseRatemap
nbins_rate = 20;

% convert amp to power and normalize by z-scoring. dicided to include
% entire signal and not just epochs of high power
sig_pow = NormToInt(log10(abs(sig_z)), 'Z', [0 Inf], fs);
sig_phase = mod(angle(sig_z), 2 * pi);

% get power and phase of each spike for each cell
spikes.pow = cellfun(@(x) interp1(tstamps, sig_pow, x, 'nearest'),...
    spktimes, 'uni', false);
spikes.phase = cellfun(@(x) interp1(tstamps, sig_phase, x, 'nearest'),...
    spktimes, 'uni', false);

% create bins for power and phase
phase_edges = linspace(0, 2 * pi, nbins_rate + 1);
ratemap.phase_bins = phase_edges(1 : end - 1) + 0.5 .* diff(phase_edges(1 : 2));
pow_edges = linspace(-1.8, 1.8, nbins_rate + 1);
ratemap.power_bins = pow_edges(1 : end - 1) + 0.5 .* diff(pow_edges(1 : 2));
pow_edges(1) = -Inf; pow_edges(end) = Inf;

% calculate occupance; duration [sec] of the signal for each pair of phase
% and power
ratemap.occupancy = ...
    histcounts2(sig_pow, sig_phase, pow_edges, phase_edges) / fs;

% for each cell, count spikes in each bin of power and phase
ratemap.counts = cellfun(@(x, y) histcounts2(x, y, pow_edges, phase_edges),...
    spikes.pow, spikes.phase, 'uni', false);

% normalize counts to rate by dividing with occupancy
ratemap.rate = cellfun(@(x) x ./ ratemap.occupancy,...
    ratemap.counts, 'uni', false);

% convert to 3d mat (power x phase x cell)
ratemap.counts = cat(3, ratemap.counts{:});
ratemap.rate = cat(3, ratemap.rate{:});

% graphics. THIS SHOULD BE DIVIDED PER POPULATION
figure
    
% mean rate across all cells
subplot(2, 2, 1)
imagesc(ratemap.phase_bins, ratemap.power_bins,...
    mean(ratemap.rate, 3, 'omitnan'))
hold on
imagesc(ratemap.phasebins + 2 * pi, ratemap.power_bins,...
    mean(ratemap.rate, 3, 'omitnan'))
plot(linspace(0, 2 * pi, 100), cos(linspace(-pi, 2 * pi, 100)), 'k')
xlim([0 2 * pi])
axis xy
colorbar
xlabel('Phase');
ylabel('Norm. Power')
title('Mean Rate')

% histogram of power occupancy
subplot(2, 2, 3)
bar(ratemap.power_bins, sum(ratemap.occupancy, 2))
xlabel('Norm. Power')
ylabel('Time [sec]')
box off
axis tight
title('Occupancy')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get significance of spkphase by jittering spktimes. this takes an insane
% amount of time. however, currently it seems that with Rayleigh all cells
% have significant modulation of phase. need to dive deeper
if jitterSig
    
    njitt = 30;
    jitterwin = 2 / frange(1);
    jitterbuffer = zeros(nunits, nfreq, njitt);
    
    for ijitt = 1 : njitt
        if mod(ijitt, 10) == 1
            display(['Jitter ', num2str(ijitt), ' of ', num2str(njitt)])
        end
        jitterspikes = bz_JitterSpiketimes(spktimes, jitterwin);
        jitt_lfp = cellfun(@(X) interp1(tstamps, sig_z, X, 'nearest'),...
            jitterspikes, 'UniformOutput', false);
        
        for iunit = 1 : nunits
            spkz = mean(jitt_lfp{iunit}, 1, 'omitnan');
            jitterbuffer(iunit, :, ijitt) = abs(spkz);
        end
        
    end
    jitt_mean = mean(jitterbuffer, 3);
    jitt_std = std(jitterbuffer, [], 3);
    spkphase_jitt = (spkmag - jitt_mean) ./ jitt_std;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wide band metrices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filter lfp at nfreq frequencies (log spaced) and normalize to mean
frange = [0.5 120];
nfreq = 20;
tmp = bz_WaveSpec(sig, 'showprogress', true,...
    'ncyc', 7, 'nfreqs', nfreq, 'frange', frange,...
    'samplingRate', 1250);
freqs = tmp.freqs;
sig_z = tmp.data ./ mean(abs(tmp.data), 1, 'omitnan');
clear tmp

% -------------------------------------------------------------------------
% correlations 
% rate modulation by counting the number of spikes in 5 ms
% bins with a moving average of 4 bins, then calculating the spearman
% correlation between counts and lfp magnitude. seems to me that stepsize
% for spkcounts should be adjusted according to the frequency of interest.

% mat of spike counts
window = 0.02;
stepsize = 0.005;
spkcnts = bz_SpktToSpkmat(spktimes, 'binsize', window, 'dt', stepsize);

% get lfp mag/phase of each frequency during each bin of spkcounts
spkcnts.lfp_filt = interp1(tstamps, sig_z, spkcnts.timestamps, 'nearest');

% rate-mag modulation. 
[wideband.ratemag_r, wideband.ratemag_p] = corr(spkcnts.data, abs(spkcnts.lfp_filt),...
    'type', 'spearman', 'rows', 'complete');

% -------------------------------------------------------------------------
% population synchrony
% must be computed separately for rs and fs cells
for ipop = 1 : length(pop.names)
    
    pop.popidx(:, ipop) = strcmp(pop.subpops, pop.names{ipop});
    
    % calc population synchrony by averaging for each time bin the number
    % of cells that were active. for this, spkcounts is converted to binary.
    popsyn = sum(spkcnts.data(:, pop.popidx(:, ipop)) > 0, 2) ./...
        length(pop.popidx(:, ipop));
    
    % standardize to the mean
    popsyn = popsyn ./ mean(popsyn);
    
    % correlation between synchrony and lfp mag for each frequency
    [pop.synmag_r(:, ipop), pop.synmag_p(:, ipop)] =...
        corr(popsyn, abs(spkcnts.lfp_filt), 'type', 'spearman',...
        'rows', 'complete');
    
    % synchrony phase coupling. not sure i understand why incorporating the
    % synchrony in the complex number provides the coupling, and could not
    % find a reference to this
    popz = mean(abs(spkcnts.lfp_filt) .*...
        (popsyn .* exp(1i .* angle(spkcnts.lfp_filt))), 1, 'omitnan');
    pop.synmag(:, ipop) = abs(popz);
    pop.synphase(:, ipop) = angle(popz);
    
end
wideband.pop = pop;

% -------------------------------------------------------------------------
% spike triggered average and pairwise phase consistency. this is perhaps
% best done by fieldtrip, but their data structure is impossible at the
% moment. here i examined building an sta via fmatoolbox. very slow
% durWin = [-0.5, 0.5];   % [sec]
% nbinsMap = floor(fs * diff(durWin) / 2) * 2 + 1; % must be odd
% [r, i] = Sync([tstamps' sig], spktimes{iunit}, 'durations', durWin);
% tmpmap = SyncMap(r, i, 'durations', durWin,...
%     'nbins', nbinsMap, 'smooth', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort units by the number of spikes or the strength of the spkphase
% coupling. only the figure and not the struct output will be sorted
if srtUnits
    [~, srtOrder] = sort(spkphase(:, 1));
    srtOrder = [1 : nunits];
end

% select frequency with most coupling
[~, ifreq] = max(sum(spkphase, 1, 'omitnan'));

% organize struct
spklfp.info.runtime     = datetime(now, 'ConvertFrom', 'datenum');
spklfp.info.winCalc     = winCalc;
spklfp.info.ch          = ch;
spklfp.freqs            = freqs;
spklfp.pop              = pop;
spklfp.spkmag           = spkmag;
spklfp.spkphase         = spkphase;
spklfp.ratemag_r        = ratemag_r;
spklfp.ratemag_p        = ratemag_p;

if saveVar
    save(spklfpfile, 'spklfp')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    setMatlabGraphics(false)
    fh = figure;
    posnegcolor = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
    
    sb1 = subplot(2, 3, 1);     % syn mag correlation
    sb2 = subplot(2, 3, 2);     % cell map of rate mag coupling
    sb3 = subplot(2, 3, 3);     % polar of spk phase/mag coupling
    sb4 = subplot(2, 3, 4);     % syn phase coupling
    sb5 = subplot(2, 3, 5);     % cell map of spk phase coupling
    sb6 = subplot(2, 3, 6);     % rose histogram of spk phse coupling
    
    % syn mag correlation
    fh.CurrentAxes = sb1;
    hold on
    plot(sb1, freqs, spklfp.pop.synmag_r)
    plot(sb1, freqs([1 end]), [0 0], 'k--')
    axis tight
    box off
    set(gca, 'xscale', 'log')
    title('Syn - Mag Correlation')
    xlabel('Frequency [Hz]');
    ylabel('Rho')
    
    % syn phase coupling
    fh.CurrentAxes = sb4;
    hold on
    plot(freqs, spklfp.pop.synphase)
    axis tight
    box off
    set(gca, 'xscale', 'log')
    xlabel('Frequency [Hz]');
    ylabel('Mean Phase')
    ylabel('Syn - Phase Coupling')
    
    % cell map of rate mag coupling
    fh.CurrentAxes = sb2;
    imagesc(freqs, 1 : nunits, real(spklfp.ratemag_r(srtOrder, :)))
    xlabel('Frequency [Hz]');
    ylabel('Cell')
    axis tight
    set(gca, 'xscale', 'log')
    colormap(gca, posnegcolor)
    title('Spike - Mag Correlation')
    
    % cell map of spk phase coupling
    fh.CurrentAxes = sb5;
    imagesc(freqs, 1 : nunits, spklfp.spkphase(srtOrder, :))
    xlabel('Frequency [Hz]');
    ylabel('Cell')
    set(gca, 'xscale', 'log')
    title('Spike - Phase Coupling')
    
    % polar of spk phase/mag coupling
    fh.CurrentAxes = sb3;
    for ipop = 1 : length(pop.names)
        ph = polar(spklfp.spkphase(spklfp.pop.popidx(:, ipop), ifreq),...
            spklfp.spkmag(spklfp.pop.popidx(:, ipop), ifreq), '.');
        hold on
    end
    title(sprintf('Spk Coupling at %.1f [Hz]', freqs(ifreq)))
    
    % rose histogram of spk phse coupling
    fh.CurrentAxes = sb6;
    for ipop = 1 : length(pop.names)
        rh{ipop} = rose(spklfp.spkphase(spklfp.pop.popidx(:, ipop), ifreq));
        rh{ipop}.Color(4) = 0.5;
        hold on
    end
    % delete1 grid lines
    h = findall(gca, 'type', 'line');
    for ipop = 1 : length(pop.names)
        h(h == rh{ipop}) = [];
    end
    delete(h);   
   
    % save
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%.spk_lfp', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF

