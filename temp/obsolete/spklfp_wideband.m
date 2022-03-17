function [spklfp] = spklfp_wideband(varargin)

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

% alternatives: 
% (1) coherencypt from chronux which calculates the cross
% multitaper spectroms of the (point-process) spktimes and (contineous)
% lfp. see http://www-users.med.cornell.edu/~jdvicto/pdfs/pubo08.pdf 
% (2) the sta method of Vinck 2012 which is implemented in fieldtrip. see
% https://www.fieldtriptoolbox.org/tutorial/spikefield/ however, i think
% both of these methods are mainly suitable for tasks with repeated trials
% i think


% currently only a single lfp channel is used. need to think how to improve
% this. in fieldtrip they mention bleeding of a spike waveform energy into
% the lfp recorded on the same channel. could be problamatic 

% INPUT
%   basepath    char. fullpath to recording folder {pwd}
%   sig         lfp filtered in the frequency band
%   spktimes    cell array of spike times
%   fs          sampling frequency of lfp data
%   frange      2 x 1 numeric of passband frequency range. note that since
%               spike are counted in 5 ms bins, the max frequency should
%               be <= 100 Hz (i think)
%   nfreq       number of frequencies for calculation
%   saveVar     logical {true}
%   graphics    logical {true}
%   srtUnits    logical {true}. currently not implemented
%
% CALLS
%   bz_SpktToSpkmat
%
% TO DO LIST
%
% 25 feb 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'spktimes', {}, @iscell)
addParameter(p, 'frange', [1 100], @isnumeric)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'nfreq', 30, @isnumeric)
addParameter(p, 'srtUnits', true, @islogical)
addParameter(p, 'graphics', true, @islogical)
addParameter(p, 'saveVar', true, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
sig             = p.Results.sig;
spktimes        = p.Results.spktimes;
frange          = p.Results.frange;
fs              = p.Results.fs;
nfreq           = p.Results.nfreq;
srtUnits        = p.Results.srtUnits;
graphics        = p.Results.graphics;
saveVar         = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params for wavelet filtering
ncyc = 7;

% params for spk count mat
window = 0.02;
stepsize = 0.005;

% files
[~, basename] = fileparts(basepath);
cd(basepath)
unitsfile = fullfile(basepath, [basename, '.units.mat']);
spklfpfile = fullfile(basepath, [basename, '.spklfp.mat']);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fit timestamps to lfp data
tstamps = 1 / fs : 1 / fs : length(sig) / fs;

% filter, norm, and hilbert. for multiple frequencies bz uses the wavelet
% transform i assume to allow log-spaced frequencies and to catch abrupt
% changes. still, perhaps consider switching to multitaper spec
tmp = bz_WaveSpec(sig, 'showprogress', true,...
    'ncyc', ncyc, 'nfreqs', nfreq, 'frange', frange,...
    'samplingRate', 1250);
freq = tmp.freqs;
sig_z = tmp.data ./ mean(abs(tmp.data), 1, 'omitnan');
clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rate - magnitude correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rate modulation by counting the number of spikes in 5 ms
% bins with a moving average of 4 bins, then calculating the spearman
% correlation between counts and lfp magnitude. seems to me that stepsize
% for spkcounts should be adjusted according to the frequency of interest.

% mat of spike counts
spkcnts = bz_SpktToSpkmat(spktimes, 'binsize', window, 'dt', stepsize);

% get lfp mag/phase of each frequency during each bin of spkcounts
lfp_z = interp1(tstamps, sig_z, spkcnts.timestamps, 'nearest');

% rate mag correlation. 
[ratemag.r, ratemag.p] = corr(spkcnts.data, abs(lfp_z),...
    'type', 'spearman', 'rows', 'complete');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% population synchrony
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
nbins_phase = 180;
phase.dist      = nan(nbins_phase, nfreq, length(pop.names));
phase.kappa     = nan(nfreq, length(pop.names));
phase.theta     = nan(nfreq, length(pop.names));
phase.r         = nan(nfreq, length(pop.names));
phase.p         = nan(nfreq, length(pop.names));

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
    % number but i think this is a mistake. adapted the same procedure as
    % for a single cell
    pop_z = abs(lfp_z) .* (popsyn .* exp(1i .* angle(lfp_z)));
    pop_angle = mod(angle(pop_z), 2 * pi);
    for ifreq = 1 : length(freq)
        [phase.dist(:, ifreq, ipop), phase.bins, tmp] =...
            CircularDistribution(pop_angle(:, ifreq), 'nBins', nbins_phase);
        phase.kappa(ifreq, ipop) = tmp.k;
        phase.theta(ifreq, ipop) = tmp.m;
        phase.r(ifreq, ipop) = tmp.r;
        phase.p(ifreq, ipop) = tmp.p;
    end   
end
pop.phase = phase;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ppc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% organize struct
spklfp.info.runtime     = datetime(now, 'ConvertFrom', 'datenum');
spklfp.freqs            = freq;
spklfp.pop              = pop;
spklfp.ratemag          = ratemag;

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
    sb2 = subplot(2, 3, 2);     % cell map of rate mag correlation
    sb3 = subplot(2, 3, 3);     % polar of spk phase/mag coupling
    sb4 = subplot(2, 3, 4);     % syn phase coupling
    sb5 = subplot(2, 3, 5);     % cell map of spk phase coupling
    sb6 = subplot(2, 3, 6);     % rose histogram of spk phse coupling
    
    % syn mag correlation
    fh.CurrentAxes = sb1;
    hold on
    plot(sb1, freq, spklfp.pop.synmag_r)
    plot(sb1, freq([1 end]), [0 0], 'k--')
    axis tight
    box off
    set(gca, 'xscale', 'log')
    title('Syn - Mag Correlation')
    xlabel('Frequency [Hz]');
    ylabel('Rho')
    
    % syn phase coupling
    fh.CurrentAxes = sb4;
    hold on
    plot(freq, pop.phase.r)
    axis tight
    box off
    set(gca, 'xscale', 'log')
    xlabel('Frequency [Hz]');
    ylabel('Mean Resultant Length')
    ylabel('Syn - Phase Coupling')
    
    % cell map of rate mag correlation
    fh.CurrentAxes = sb2;
    imagesc(freq, 1 : nunits, spklfp.ratemag.r)
    xlabel('Frequency [Hz]');
    ylabel('Cell')
    axis tight
    colormap(gca, posnegcolor)
    title('Spike - Mag Correlation')
    
    % cell map of rate mag correlation
    fh.CurrentAxes = sb5;
    data = bz_NormToRange(squeeze(pop.phase.dist(2 : end, 1 : end, 2)), [0 1]);
    imagesc(pop.phase.bins, freq, data)
    xlabel('Phase');
    ylabel('Frequency')
    axis tight
    colormap
    title('Synchrony - Phase Correlation')    
    
    % save
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%.spk_lfp', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF

