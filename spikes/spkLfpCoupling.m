function [spklfp] = spkLfpCoupling(varargin)

% investigate spk - lfp coupling. based on bz_GenSpikeLFPCoupling.
% main differences from bz_GenSpikeLFPCoupling: (1) original also had a
% subrouting for calculating mutual information between power and ISI
% distribution which i ignored. could be relavent for theta-driving
% interneurons (Li, Sci. Reports, 2017). (2) original loopes across lfp
% channels, probably for linear probes. otherwise this version is simply
% more organized. 

% should also implement coherencypt from chronux which calculates the
% cross-spectrom of the multitaper specs (point-process) spktimes and
% (contineous) lfp. also should inspect the sta method of Vinck 2012 which is
% implemented in fieldtrip

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

% restrict to winCalc
spktimes = cellfun(@(x) x(InIntervals(x, winCalc)), spikes.times, 'uni', false);

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
% spike count matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a matrix of spike counts (each row is a differnet cell). the time
% resolution is binsize and the overlap for moving average is window /
% stepsize.
window = 0.02;
stepsize = 0.005;
spkcnts = bz_SpktToSpkmat(spktimes, 'binsize', window, 'dt', stepsize);

% if computation takes too long, can add a limit to the number of spikes
% here

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
        sig_filt = 	hilbert(sig_filt);
        sig_filt = sig_filt ./ mean(abs(sig_filt));
        
    otherwise
        tmp = bz_WaveSpec(sig, 'showprogress', true,...
            'ncyc', 7, 'nfreqs', nfreq, 'frange', frange,...
            'samplingRate', 1250);
        freqs = tmp.freqs;
        sig_filt = tmp.data ./ mean(abs(tmp.data), 1, 'omitnan');
        clear tmp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restrict anaylsis to epochs with enough power in frange
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% approximate power in frequency band
power = fastrms(sig_filt, ceil(fs ./ frange(1)), 1);
powThr = 2;
minpow = mean(power) + std(power) * powThr;

% set the minimum duration to two cycles
mindur = (fs ./ frange(2)) * 2; 
    
% if working on multiple frequencies simultaneously, must run this in a
% loop
% find epochs with power > low threshold. correct for durations
bad_epochs = binary2epochs('vec', power < minpow, 'minDur', mindur,...
    'maxDur', 0, 'interDur', 0);
nepochs = size(bad_epochs, 1);

% subtract out low power intervals
intervals = SubtractIntervals(winCalc, bad_epochs);  


angle(sig_filt) - mod(angle(sig_filt),2*pi)
        lfpphase = mod(angle(hilb),2*pi);

% get the lfp mag/phase of each frequency during each bin of spkcounts
spkcnts.lfp_filt = interp1(tstamps, sig_filt, spkcnts.timestamps, 'nearest');

% for each cell, get the lfp mag/phase of each frequency during each spike
for iunit = 1 : nunits
    bz_Counter(iunit, nunits, 'Interpolating Cell')
    spikes.lfp_filt{iunit} = interp1(tstamps, sig_filt, spktimes{iunit},'nearest');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rate-magnitue modulation. correlation between spkcounts in each bin to
% the magnitude of the lfp in each frequecny band. it seems to me that
% stepsize for spkcounts should be adjusted according to the frequency of
% interest.
[ratemag_r, ratemag_p] = corr(spkcnts.data, abs(spkcnts.lfp_filt),...
    'type', 'spearman', 'rows', 'complete');

% spk-phase coupling
spkz = cellfun(@(x) mean(x, 1, 'omitnan'), spikes.lfp_filt, 'uni', false);
spkz = cell2mat(spkz');
spkmag = abs(spkz);
spkphase = angle(spkz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get significance of spkphase by jitterign spktimes. this takes an insane
% amount of time
if jitterSig
    
    njitt = 30;
    jitterwin = 2 / frange(1);
    jitterbuffer = zeros(nunits, nfreq, njitt);
    
    for ijitt = 1 : njitt
        if mod(ijitt, 10) == 1
            display(['Jitter ', num2str(ijitt), ' of ', num2str(njitt)])
        end
        jitterspikes = bz_JitterSpiketimes(spktimes, jitterwin);
        jitt_lfp = cellfun(@(X) interp1(tstamps, sig_filt, X, 'nearest'),...
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
% population synchrony
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Population Synchrony: Phase Coupling and Rate Modulation
% here is where the separation of rs and fs is useful
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
    figname = fullfile(figpath, sprintf('%spk_lfp', basename));
    export_fig(figname, '-tif', '-transparent', '-r300')
end

end

% EOF

