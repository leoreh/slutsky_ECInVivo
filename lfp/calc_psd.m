function [psd, faxis, psd_bouts] = calc_psd(varargin)

% calculates the psd of a signal with the multitaper estimate by chronux.
% the spectrum is calculated in specific timebins, or the timebins are
% concatenated and the psd is calculated once (see alternative methods
% below). timebins are provided as a cell array such that the psd
% will be calculated for each cell separately.
% -------------------------------------------------------------------------
% notes (see accompanying presentation).
% (1) with mtspectrumsegc, each timebin is divided into segments and the
% psd is averaged across segments. this permits the frequency resolution
% to be independent of the timebin length (as opposed to mspectrumc).
% (2) with mtspectrumc there is much higher variability in high frequencies
% (3) averaging a tetrode and calculating the psd provides a similar result
% to calculating the psd per electrode and averaging the psds. the results
% are different then calculating the psd based on a single channel.
% -------------------------------------------------------------------------
% see also Savolainen, scientific reports, 2021 and K. van Vugt et al., J.
% neurosci. met., 2007
% -------------------------------------------------------------------------
%
% INPUT:
%   sig             numeric. lfp / eeg data (1 x n)
%   fs              numeric. sampling frequency {1250} of the signal.
%   ftarget         numeric. requested frequencies of psd estimate
%   bins            cell of n x 2 mats. the psd will be calculated for each
%                   cell separately according to the time windows specified
%                   in each cell [sec]. e.g., bins = ss.boutTimes.
%   graphics        logical. plot figure {true}
% 
% OUTPUT
%   psd             numeric nbins x ftarget. raw averaged psd for each bin 
%   faxis           numeric frequency values
%   psd_bouts       cell of nwins with the psd per bin
%
% DEPENDENCIES
%   mtspectrumc
% 
% TO DO LIST
%   # implement fCWT https://www.nature.com/articles/s43588-021-00183-z
%   see also https://www.nature.com/articles/s41598-021-86525-3#Sec4
%   # separate the oscillations from the power law
%   https://link.springer.com/article/10.1007/s12021-022-09581-8    
%
% 26 jul 21 LH      updates:
% 05 jan 21         removed downsampling
% 12 jan 22         bins. removed state dependancy
% 07 sep 22         applied mtspectrumc by chronux 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'sig', [], @isnumeric);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'ftarget', [], @isnumeric);
addOptional(p, 'bins', []);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
sig             = p.Results.sig;
fs              = p.Results.fs;
ftarget         = p.Results.ftarget;
bins            = p.Results.bins;
graphics        = p.Results.graphics;

if isempty(ftarget)
    ftarget = [1 : 0.1 : 100];
end

% margin to remove from the start and end of an bout 
sigMrgin = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assert time bins
if isempty(bins)
    bins = {[0 Inf]};
end
if ~iscell(bins)
    if size(bins, 2) ~= 2
        error('check the input bins')
    end
    bins = {bins};
end

maxsec = floor(length(sig) / fs);
if ~isempty(bins{end})
    if bins{end}(end) > maxsec
        bins{end}(end) = maxsec;
    end
end

% get number of bins
ncells = size(bins, 2);
nbins = cellfun(@(x) size(x, 1), bins, 'uni', true);
for icell = 1 : ncells
    if nbins(icell) < 1
        fprintf('\nWARNING: some timebins are empty\n')
        bins{icell} = [];
    end
end

% get minimum bin length
minbin = min(cellfun(@(x) min(diff(x')), bins(nbins > 0), 'uni', true));

% prep frequencies
frange = [ftarget(1), ftarget(end)];
if frange(2) > fs / 2
    error('requested max frequency greater than nyquist')
end

% chronux params
window = min([5, minbin]);     % [s]
segave = 1;                    % in mtspectrumsegc, average across segments
mtspec_params.pad = 0;
mtspec_params.Fs = fs;
mtspec_params.fpass = frange;
mtspec_params.tapers = [ceil(window / 2), ceil(window / 2) * 2 - 1];
mtspec_params.tapers = [3, 5];
mtspec_params.trialave = 1;
mtspec_params.err = [0, 0];

% frequency grid from chronux (for initializing). only used when sig is divided
% to bins (mtspectrumsegc)
N = window * fs;
nfft = max(2 ^ (nextpow2(N) + mtspec_params.pad), N);
[f, ~] = getfgrid(fs, nfft, frange); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc power for each bout separatly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alt = 1;

% initialize
if alt < 3
    psd = nan(ncells, length(f));
else
    psd = nan(ncells, length(ftarget));
end
psd_bouts = cell(ncells, 1);

for icell = 1 : ncells
    
    if nbins(icell) == 0
        continue
    end

    % ---------------------------------------------------------------------
    % ALT 1: calc psd per bout and average across bouts (mspectrumsegc)
    if alt == 1        
        
        % initialize
        psd_bouts{icell} = nan(nbins(icell), length(f));
        
        % loop through bouts
        for ibin = 1 : nbins(icell)

            % idx to signal w/o edges (sigMrgin)
            dataIdx = (bins{icell}(ibin, 1) + sigMrgin) * fs :...
                (bins{icell}(ibin, 2) - sigMrgin) * fs - 1;

            if length(dataIdx) < window * fs
                nbins(icell) = nbins(icell) - 1;
                warning('some bouts too short for fft window')
                continue
            end

            % calc power and sum across bins
            [psd_bouts{icell}(ibin, :), f] = mtspectrumsegc(sig(dataIdx), window, mtspec_params, segave);

        end
        
        % average power across bouts
        psd(icell, :) = mean(psd_bouts{icell}, 1, 'omitnan');

    % ---------------------------------------------------------------------
    % ALT 2: cat bouts and calc psd once (mspectrumsegc)
    elseif alt == 2
        dataIdx = Restrict([1 : length(sig)], bins{icell} * fs);
        [psd(icell, :), f] = mtspectrumsegc(sig(dataIdx), window, mtspec_params, segave);
    
    % ---------------------------------------------------------------------
    % ALT 3: cat bouts and calc psd once (mspectrumc). doesn't work good. 
    elseif alt == 3
        dataIdx = Restrict([1 : length(sig)], bins{icell} * fs);
        [tmp, f] = mtspectrumc(sig(dataIdx), mtspec_params);
    
        % adjust the fequency domain
        fidx = zeros(1, length(ftarget)); % find closest indices in f
        for ifreq = 1 : length(ftarget)
            [fdev(ifreq), fidx(ifreq)] = min(abs(f - ftarget(ifreq)));
        end
        psd(icell, :) = tmp(fidx);
        faxis = f(fidx);

    % ---------------------------------------------------------------------
    % ALT 4: calc psd per bout and average across bouts (mspectrumc)
    elseif alt == 4
        
        psd_bouts{icell} = nan(nbins(icell), length(ftarget));
        for ibin = 1 : nbins(icell)

            % idx to signal w/o edges (sigMrgin)
            dataIdx = (bins{icell}(ibin, 1) + sigMrgin) * fs :...
                (bins{icell}(ibin,     2) - sigMrgin) * fs - 1;

            if length(dataIdx) < window * fs
                nbins(icell) = nbins(icell) - 1;
                warning('some bouts too short for fft window')
                continue
            end

            % calc power and sum across bins
            [tmp, f] = mtspectrumc(sig(dataIdx, :), mtspec_params);
            
            % adjust the fequency domain
            fidx = zeros(1, length(ftarget)); % find closest indices in f
            for ifreq = 1 : length(ftarget)
                [fdev(ifreq), fidx(ifreq)] = min(abs(f - ftarget(ifreq)));
            end
            psd_bouts{icell}(ibin, :) = tmp(fidx);
            faxis = f(fidx);

        end
        % average power
        psd(icell, :) = mean(psd_bouts{icell}, 1, 'omitnan');
    end
end

faxis = f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectral power per window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use 10*log10(psdStates) for [dB]

if graphics
    xLimit = [0 ftarget(end)];
    fh = figure;
    
    % raw psd
    sb1 = subplot(1, 2, 1);
    ph = plot(faxis, psd, 'LineWidth', 3);
    xlim(xLimit)
    xlabel('Frequency [Hz]')
    ylabel('PSD [mV^2/Hz]')
    set(gca, 'YScale', 'log')
    
    % norm psd
    sb2 = subplot(1, 2, 2);
    ph = plot(faxis, psd ./ sum(psd, 2), 'LineWidth', 3);
    xlim(xLimit)
    xlabel('Frequency [Hz]')
    ylabel('norm PSD')
    set(gca, 'YScale', 'log')
end

end

% EOF