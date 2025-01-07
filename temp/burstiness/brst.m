function   all_su_bursts = brst(data, isiThr, minSpksBrst) 

% calculates the multitaper spectrogram via chronux with refinements from
% accusleep. data can be loaded from .lfp or provided as input.
% multichannel support includes averaging channels (e.g. electrodes from
% the same tetrode) or avergaing spectrograms across channels and/or
% calculating the spectrogram separately for groups of channels. these are
% determined by the input 'ch'
%
% INPUT
%   basepath    char. fullpath to recording folder {pwd}
%   saveVar     logical. organize and save struct {true}
% 
% OUTPUT
%
% CALLS
%   binary2bouts
% 
% TO DO LIST
%
% 05 feb 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'saveVar', true, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
saveVar         = p.Results.saveVar;

spktimes = spikes.times;
isiThr = 0.02;  % units must be as spktimes
binsize = 3600;  
minSpksBrst = 2;    % minimum no. spikes in a burst

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nunits = length(spktimes);

% separate recording to bins
recLen = max(cellfun(@max, spktimes, 'uni', true));
bins = n2chunks('n', recLen, 'chunksize', binsize, 'overlap', [1 0]);
bins(1) = 0;
nbins = size(bins, 1);

% initialize
clear brst
brst.avg.nspks = zeros(nbins, nunits);
brst.avg.dur = zeros(nbins, nunits);
brst.avg.freq = zeros(nbins, nunits);
brst.avg.spkprct = zeros(nbins, nunits);
brst.avg.ibi = zeros(nbins, nunits);    
brst.nbrsts.detect = zeros(nbins, nunits);
brst.nbrsts.dur = zeros(nbins, nunits);
brst.nbrsts.freq = zeros(nbins, nunits);
brst.nbrsts.freqNorm = zeros(nbins, nunits);
brst.nbrsts.shortPrct = zeros(nbins, nunits);

% make cell array of spikes per unit per time
for ibin = 1 : nbins
    for iunit = 1 : nunits
        
        spks = spktimes{iunit}(InIntervals(spktimes{iunit}, bins(ibin, :)));
        isi = diff(spks);
        nspks = length(spks);

        % get indices of first and last spike in each burst and force minSpksBrst
        [b.idx, nbouts] = binary2bouts('vec', isi <= isiThr,...
            'minDur', minSpksBrst - 1, 'printFlag', false);

        if isempty(b.idx)
            continue
        end

        % ratio of bursts that pass the min spks criterion
        nbrsts.shortPrct = (nbrsts.detect - nbrsts.dur) / nbrsts.detect * 100;
        
        % number of bursts relative to binsize and ~firing rate
        nbrsts.freq = nbrsts.dur / binsize;
        nbrsts.freqNorm = nbrsts.freq / (nspks / binsize);

        % number of spikes per burst
        b.nspks = diff(b.idx, 1, 2) + 1;
        
        % start and stop times per burst
        b.times = spks(b.idx(:, 1));
        b.times(:, 2) = spks(b.idx(:, 2));
        
        % burst duration [s]
        b.dur = diff(b.times, 1, 2); 
        
        % average frequency within burst
        b.freq = b.nspks ./ (b.dur);
        
        % percent of spikes that participate in bursts
        b.spkprct = sum(b.nspks) / nspks;
        
        % inter-burst interval
        b.ibi = b.times(2 : end, 1) - b.times(1 : end - 1, 2);
        
        % organize per bin stats
        brst.all(ibin, iunit) = b;

        % average stats
        brst.avg.nspks(ibin, iunit) = mean(b.nspks);
        brst.avg.dur(ibin, iunit) = mean(b.dur);
        brst.avg.freq(ibin, iunit) = mean(b.freq);
        brst.avg.spkprct(ibin, iunit) = mean(b.spkprct);
        brst.avg.ibi(ibin, iunit) = mean(b.ibi);
        
        % number of bursts stats
        brst.nbrsts.detect(ibin, iunit) = nbrsts.detect;
        brst.nbrsts.dur(ibin, iunit) = nbrsts.dur;
        brst.nbrsts.freq(ibin, iunit) = nbrsts.freq;
        brst.nbrsts.freqNorm(ibin, iunit) = nbrsts.freqNorm;
        brst.nbrsts.shortPrct(ibin, iunit) = nbrsts.shortPrct;

    end
end

end

% EOF