function [spklfp] = spklfp_singleband(varargin)

% investigate spk - lfp coupling as point to fields in a single frequency
% band. metrices include (1) phase coupling per cell including significance
% test (Rayleigh or jittering), and rate map for each cell by lfp phase and
% power.

% based on bz_GenSpikeLFPCoupling, bz_PhaseModulation, and
% bz_PowerPhaseRatemap. 

% currently only a single lfp channel is used. need to think how to improve
% this. in fieldtrip they mention bleeding of a spike waveform energy into
% the lfp recorded on the same channel which could be problamatic 

% INPUT
%   basepath    char. fullpath to recording folder {pwd}
%   sig         lfp filtered in the frequency band
%   spktimes    cell array of spike times
%   fs          sampling frequency of lfp data
%   frange      2 x 1 numeric of passband frequency range
%   saveVar     logical {true}
%   graphics    logical {true}
%   srtUnits    logical {true}. currently not implemented
%
% CALLS
%   CircularDistribution (fmat)
%   NormToInt (bz)
%
% TO DO LIST
%   # jitter in a reasonable way
%
% 25 feb 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'spktimes', {}, @iscell)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'frange', [1.5 100], @isnumeric)
addParameter(p, 'srtUnits', true, @islogical)
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
sig             = p.Results.sig;
spktimes        = p.Results.spktimes;
fs              = p.Results.fs;
frange          = p.Results.frange;
srtUnits        = p.Results.srtUnits;
graphics        = p.Results.graphics;
saveVar         = p.Results.saveVar;

% params 
powThr = 2;                         % stds above mean power in band
mindur = (fs ./ frange(2)) * 2;     % minimum duration of epoch is two cycles
nbins_phase = 180;
nbins_rate = 20;

% params for spk count mat
movWin = 0.02;
stepsize = 0.005;
nunits = length(spktimes);

% files
[~, basename] = fileparts(basepath);
cd(basepath)
unitsfile = fullfile(basepath, [basename, '.units.mat']);
spklfpfile = fullfile(basepath, [basename, '.spklfp.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fit timestamps to lfp data
tstamps = 1 / fs : 1 / fs : length(sig) / fs;

% hilbert 
sig_z = hilbert(sig);

% convert amp to power and normalize by z-scoring. this is for the rate map
% which includes the entire signal and not just high power epochs
sig_pow = NormToInt(log10(abs(sig_z)), 'Z', [0 Inf], fs);

% norm to mean
% sig_z = sig_z ./ mean(abs(sig_z));

% get angles in range 0 : 2pi rather than -pi : pi. this is neccassary
% because of the way circularDistribution bins the data. 
sig_phase = mod(angle(sig_z), 2 * pi);

% restrict anaylsis to epochs with enough power in band. not sure why use
% the rms and not the normalized power. this is only for phase coupling
sig_rms = fastrms(sig, ceil(fs ./ frange(1)), 1);
minrms = mean(sig_rms) + std(sig_rms) * powThr;

% find epochs with power > low threshold. correct for durations
bad_epochs = binary2epochs('vec', sig_rms < minrms, 'minDur', mindur,...
    'maxDur', Inf, 'interDur', 0);
clear sig_rms

% remove low power intervals
lfp_epochs = SubtractIntervals([0, length(sig)], bad_epochs) / fs;  
lfp_epochs(end) = length(sig) / fs;

% high power occupancy
lfp_occupancy = sum(lfp_epochs(:, 2) - lfp_epochs(:, 1)) ./...
    (length(sig) / fs) * 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% restrict spktimes to lfp epochs. if computation takes too long, can add a
% limit to the number of spikes here

% get power and phase of each spike for each cell. in bz_PhaseModulation
% this is done by indexing angle(sig_z) with spktimes * fs which is
% faster but less accurate.
spk_pow = cellfun(@(x) interp1(tstamps, sig_pow, x, 'nearest'),...
    spktimes, 'uni', false);
spk_phase = cellfun(@(x) interp1(tstamps, sig_phase, x, 'nearest'),...
    spktimes, 'uni', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase modulation per cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
phase.dist      = nan(nbins_phase, nunits);
phase.kappa     = nan(1, nunits);
phase.theta     = nan(1, nunits);
phase.r         = nan(1, nunits);
phase.p         = nan(1, nunits);
for iunit = 1 : nunits
    
    % restrict spikes to epochs of high power. for other metrices all
    % spikes are used
    spkidx = InIntervals(spktimes{iunit}, lfp_epochs);
    if sum(spkidx) == 0 || isempty(spkidx)
        continue
    end

    % use zugaro for circular statstics. mean phase - theta; mean resultant
    % length - r; Von Mises concentraion - kappa; Rayleigh significance of
    % uniformity - p; see also:
    % https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Analysis.pdf
    % in bz_genSpikeLfpCoupling they take the angle and magnitude of
    % the mean(hilbert) instead of averaging the angle and magnitude
    % separately which i think is a mistake
    [phase.dist(:, iunit), phase.bins, tmp] =...
        CircularDistribution(spk_phase{iunit}(spkidx), 'nBins', nbins_phase);
    phase.kappa(iunit) = tmp.k;
    phase.theta(iunit) = tmp.m;
    phase.mrl(iunit) = tmp.r;
    phase.p(iunit) = tmp.p;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rate map as a function of power and phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here, a bivariate histogram is created for each cell by dividing the lfp
% power and phase to nbins x nbins and counting the number of spikes in each
% bin. the counts or normalized to rate by occupance - the duration of each
% bin. based on bz_PowerPhaseRatemap

% create bins for power and phase
phase_edges = linspace(0, 2 * pi, nbins_rate + 1);
pow_edges = linspace(-1.8, 1.8, nbins_rate + 1);
ratemap.phase_bins = phase_edges(1 : end - 1) + 0.5 .* diff(phase_edges(1 : 2));
ratemap.power_bins = pow_edges(1 : end - 1) + 0.5 .* diff(pow_edges(1 : 2));
pow_edges(1) = -Inf; pow_edges(end) = Inf;

% calculate occupance; duration [sec] of the signal for each phase-power
% bin
ratemap.occupancy = ...
    histcounts2(sig_pow, sig_phase, pow_edges, phase_edges) / fs;

% for each cell, count spikes in each bin of power and phase
ratemap.counts = cellfun(@(x, y) histcounts2(x, y, pow_edges, phase_edges),...
    spk_pow, spk_phase, 'uni', false);

% normalize counts to rate by dividing with occupancy
ratemap.rate = cellfun(@(x) x ./ ratemap.occupancy,...
    ratemap.counts, 'uni', false);

% convert to 3d mat (power x phase x cell)
ratemap.counts = cat(3, ratemap.counts{:});
ratemap.rate = cat(3, ratemap.rate{:});

% mean rate
ratemap.ratemean = squeeze(mean(ratemap.rate, 1, 'omitnan'));

% clear memory
% clear spk_pow 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rate - magnitude correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rate modulation by counting the number of spikes in 5 ms
% bins with a moving average of 4 bins, then calculating the spearman
% correlation between counts and lfp magnitude. seems to me that stepsize
% for spkcounts should be adjusted according to the frequency of interest.

% mat of spike counts
spkcnts = bz_SpktToSpkmat(spktimes, 'binsize', movWin, 'dt', stepsize);

% get lfp mag/phase of each frequency during each bin of spkcounts
lfp_z = interp1(tstamps, sig_z, spkcnts.timestamps, 'nearest');

% rate mag correlation. 
[ratemag.r, ratemag.p] = corr(spkcnts.data, abs(lfp_z),...
    'type', 'spearman', 'rows', 'complete');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% population synchrony
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get subpopulations of cells (e.g. rs and fs)
nunits = length(spktimes);
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

% initialize
pop.phase.dist      = nan(nbins_phase, length(pop.names));
pop.phase.kappa     = nan(1, length(pop.names));
pop.phase.theta     = nan(1, length(pop.names));
pop.phase.mrl         = nan(1, length(pop.names));
pop.phase.p         = nan(1, length(pop.names));

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
        corr(popsyn, abs(lfp_z), 'type', 'spearman',...
        'rows', 'complete');
    
    % synchrony phase coupling. in bz this is doen by averaging the complex
    % number but decided to adapt the same procedure as
    % for single cells
    pop_z = abs(lfp_z) .* (popsyn .* exp(1i .* angle(lfp_z)));
    pop_angle = mod(angle(pop_z), 2 * pi);
    [pop.phase.dist(:, ipop), pop.phase.bins, tmp] =...
        CircularDistribution(pop_angle, 'nBins', nbins_phase);
    pop.phase.kappa(ipop) = tmp.k;
    pop.phase.theta(ipop) = tmp.m;
    pop.phase.mrl(ipop) = tmp.r;
    pop.phase.p(ipop) = tmp.p;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort units by the number of spikes or the strength of the spkphase
% coupling. only the figure and not the struct output will be sorted
if srtUnits
%     [~, srtOrder] = sort(spk_phase);
    srtOrder = [1 : nunits];
end

% organize struct
spklfp.info.runtime         = datetime(now, 'ConvertFrom', 'datenum');
spklfp.info.frange          = frange;
spklfp.info.epochs          = lfp_epochs;
spklfp.info.occupancy       = lfp_occupancy;
spklfp.ratemap              = ratemap;
spklfp.phase                = phase;
spklfp.pop                  = pop;
spklfp.ratemag              = ratemag;

% save
if saveVar
    save(spklfpfile, 'spklfp')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    setMatlabGraphics(false)
     
    fh = figure;
    
    % mean rate across all cells
    sb1 = subplot(2, 3, 1);
    hold on
    imagesc(sb1, spklfp.ratemap.phase_bins, spklfp.ratemap.power_bins,...
        mean(spklfp.ratemap.rate, 3, 'omitnan'))
    imagesc(sb1, spklfp.ratemap.phase_bins + 2 * pi, spklfp.ratemap.power_bins,...
        mean(spklfp.ratemap.rate, 3, 'omitnan'))
    plot(sb1, linspace(0, 4 * pi, 100), cos(linspace(0, 4 * pi, 100)), 'k')
    axis tight
    ax = gca;
    xticks([0 : pi : 4 * pi])
    xticklabels(string(0 : 4) + "pi")
    axis xy
    colorbar
    xlabel('Phase [rad]');
    ylabel('Norm. Power')
    title('Mean Rate')

    % polar plot of mean phase and resultant length per cell
    subplot(2, 3, 2)
    polarplot(spklfp.phase.theta, spklfp.phase.mrl, '.')
    rlim([0 1])
    rticks([])
    thetaticks([0 : 90 : 270])
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    pax.GridAlpha = 0.2;
    title('phase and MRL per cell')   

    % polar histogram of phase distribution between cells
    subplot(2, 3, 3)
    polarhistogram(spklfp.phase.theta, 10)
    title('phase distribution')
    rticks([])
    thetaticks([0 : 90 : 270])
    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    pax.GridAlpha = 0.2;

    % histogram of power occupancy
    subplot(2, 3, 4)
    bar(spklfp.ratemap.power_bins, sum(spklfp.ratemap.occupancy, 2))
    xlabel('Norm. Power')
    ylabel('Time [sec]')
    box off
    axis tight
    title('Occupancy')   
    
    % histogram of phase modulation of population synchrony
    subplot(2, 3, 5)
    bar(spklfp.pop.phase.bins(2 : end), spklfp.pop.phase.dist(2 : end, :))
    xlabel('Phase [rad]')
    ylabel('Time [sec]')
    ax = gca;
    xlim([0 2 * pi])
    xticks([0 : pi : 2 * pi])
    xticklabels(string(0 : 2) + "pi")
    box off
    title('Synchrony phase')
    
    sgtitle(sprintf('%.1f - %.1f Hz', frange))
    
    % save
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%.spk_lfp', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
    
    % ---------------------------------------------------------------------
    % per unit
    %     sunits = randperm(nunits, 10);
    %     for iunit = 1 : 10
    %         su = sunits(iunit);
    %
    %         fh = figure;
    %
    %         subplot(1, 2, 1);
    %         rose(spklfp.phase.spk_phase{su})
    %         title(sprintf('Cell %d; Rayleigh = %.2f', iunit, phase.p(su)))
    %
    %         subplot(1, 2, 2);
    %         bar(spklfp.phase.bins * 180 / pi, spklfp.phase.dist(:, su))
    %         xlim([0 360])
    %         set(gca, 'XTick', [0 90 180 270 360])
    %         hold on;
    %         plot([0 : 360], cos(pi / 180 * [0 : 360]) * 0.05 *...
    %             max(spklfp.phase.dist(:, su)) + 0.95 * max(spklfp.phase.dist(:, su)),...
    %             'color', [.7 .7 .7])
    %     end
end

end

% EOF

