
% spklfp_wrapper

% investigate spk - lfp coupling.

% spklfp_singleband calculates for a specifiec band (e.g. theta) the phase coupling and rate
% map as a function of phase and lfp power. spklfp_wideband calculates the
% correlation between spike counts in 5 ms bins and the lfp magnitude and
% the correlation between population synchrony (in 5 ms bins) and lfp
% magnitude.

% in addition to working on multiple / single frequency, single band
% filters the data by fir whereas wide band filters each frequency band
% with by wavelet convolution. 

% based on bz_GenSpikeLFPCoupling, bz_PhaseModulation, and
% bz_PowerPhaseRatemap. alternatives include:
% (1) coherencypt from chronux which calculates the cross
% multitaper spectroms of the (point-process) spktimes and (contineous)
% lfp. see http://www-users.med.cornell.edu/~jdvicto/pdfs/pubo08.pdf 
% (2) the sta method of Vinck 2012 which is implemented in fieldtrip. see
% https://www.fieldtriptoolbox.org/tutorial/spikefield/ however, i think
% both of these methods are mainly suitable for tasks with repeated trials
% i think

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[~, basename] = fileparts(basepath);
cd(basepath)

% files
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

% load lfp
winCalc = [0 3 * 60 * 60];
recDur = diff(winCalc);
ch = 9 : 11;
sig = double(bz_LoadBinary(lfpfile, 'duration', diff(winCalc),...
    'frequency', fs, 'nchannels', nchans, 'start', winCalc(1),...
    'channels', ch, 'downsample', 1));
sig = mean(sig, 2);

% clip spike times
spktimes = cellfun(@(x) x(InIntervals(x, winCalc)),...
    spikes.times, 'uni', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single band analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frange = [0.5 1.5; 1.5, 4; 4, 6; 6, 10;];
frange = [frange; [10 : 10 : 90; 20 : 10 : 100]'];

clear spklfp
for ifreq = 1 : size(frange, 1)
    
    fprintf('working on frequency %d of %d', ifreq, size(frange, 1))
    
    % filter lfp
    sig_filt = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
        'order', 3, 'passband', frange(ifreq, :), 'graphics', false);
    
    [spklfp(ifreq)] = spklfp_singleband('basepath', basepath, 'fs', fs,...
        'sig', sig_filt, 'spktimes', spktimes, 'frange', frange(ifreq, :),...
        'srtUnits', true, 'graphics', true,...
        'saveVar', false);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wide band graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cat struct
sl = catfields(spklfp, 'catdef', 'addim');

% center frequencies 
freq = sl.info.frange(1, :) + diff(sl.info.frange) / 2;

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
plot(sb1, freq, sl.pop.synmag_r)
plot(sb1, freq([1 end]), [0 0], 'k--')
xticks([0.1, 1, 10, 100])
box off
set(gca, 'xscale', 'log')
title('Syn - Mag Correlation')
xlabel('Frequency [Hz]');
ylabel('Rho')

% syn phase coupling
fh.CurrentAxes = sb4;
hold on
plot(freq, sl.pop.phase.mrl)
xticks([0.1, 1, 10, 100])
box off
set(gca, 'xscale', 'log')
xlabel('Frequency [Hz]');
ylabel('Mean Resultant Length')

% cell map of rate mag correlation
fh.CurrentAxes = sb2;
imagesc(freq, 1 : nunits, sl.ratemag.r)
xlabel('Frequency [Hz]');
ylabel('Cell')
axis tight
axis xy
colormap(gca, posnegcolor)
title('Spike - Mag Correlation')

% map of mean rate across phases and frequencies
fh.CurrentAxes = sb5;
data = squeeze(mean(mean(sl.ratemap.rate, 1, 'omitnan'), 3, 'omitnan'));
imagesc(sl.ratemap.phase_bins(:, 1), freq, data)
ylabel('Frequency [Hz]')
xlabel('Phase [rad]')
xticks([0 : pi : 2 * pi])
xticklabels(string(0 : 2) + "pi")
title('Mean Rate Across Cells')

% mrl of all cells vs. frequency
fh.CurrentAxes = sb6;
plot_boxMean(sl.phase.mrl, 'allPnts', true)
xticklabels(string(freq))
xlabel('Frequency [Hz]')
ylabel('Mean Resultant Length')
title('MRL vs. Frequency')

% save
figpath = fullfile(basepath, 'graphics');
mkdir(figpath)
figname = fullfile(figpath, sprintf('%.spk_lfp', basename));
export_fig(figname, '-tif', '-transparent', '-r300')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% buzsaki
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wide band analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frange = [1 100];
nfreq = 20;
[spklfp] = spklfp_wideband('basepath', basepath, 'fs', fs,...
    'sig', sig, 'spktimes', spktimes, 'frange', frange,...
    'nfreq', nfreq, 'srtUnits', true, 'graphics', true,...
    'saveVar', true);

