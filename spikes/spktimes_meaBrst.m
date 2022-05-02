function brst = spktimes_meaBrst(spktimes, varargin) 

% calculates burstiness stats per unit in timebins. based on analysis done
% for mea recordings
%
% INPUT
%   spktimes    cell of spike times per unit. typically in [s]. can also be
%               samples, ms, etc. but then binsize, and isiThr must be the
%               same units
%   binsize     numeric. duration of bins {3600 [s]}
%   bins        cell array of n x 2 mats of intervals.
%               metrices will be calculated for each cell by limiting
%               spktimes to the intervals. can be for example
%               ss.stateEpochs. must be the same units as spikes.times
%               (e.g. [s]). will override binsize
%   isiThr      numeric. minimum inter-spike interval for defining a burst
%               {0.02 [s]}
%   minSpks     numeric. minimum number of spikes for defining a burst {2}
%   
%   basepath    char. fullpath to recording folder {pwd}
%   saveVar     logical. save struct {true}
%   force       logical. force analyze even if file exists {false}
% 
% OUTPUT
%   brst        struct
%
% CALLS
%   binary2epochs
%   n2chunks
% 
% TO DO LIST
%
% 05 feb 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'isiThr', 0.02, @isnumeric)
addParameter(p, 'binsize', 3600, @isnumeric)
addParameter(p, 'bins', [])
addParameter(p, 'minSpks', 2, @isnumeric)
addParameter(p, 'saveVar', true, @islogical)
addParameter(p, 'force', false, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
isiThr          = p.Results.isiThr;
binsize         = p.Results.binsize;
bins            = p.Results.bins;
minSpks         = p.Results.minSpks;
saveVar         = p.Results.saveVar;
force           = p.Results.force;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% separate recording to bins
recLen = max(cellfun(@max, spktimes, 'uni', true));
if isempty(bins)
    if binsize > recLen
        binsize = recLen;
        bins = [0, binsize];
    else
        bins = n2chunks('n', recLen, 'chunksize', binsize, 'overlap', [1 0]);
        bins(1) = 0; bins(end - 1, 2) = recLen; bins(end, :) = [];
    end
end
if ~iscell(bins)
    bins = mat2cell(bins, ones(size(bins, 1), 1), size(bins, 2));
end
nbins = size(bins, 1);

% load if exists
[~, basename] = fileparts(basepath);
brstfile = fullfile(basepath, [basename, '.brst.mat']);
if exist(brstfile, 'file') && ~force
    load(brstfile, 'brst')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc burstiness per unit per timebin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nunits = length(spktimes);

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

for ibin = 1 : nbins
    for iunit = 1 : nunits
        
        spks = spktimes{iunit}(InIntervals(spktimes{iunit}, bins{ibin}));
        isi = diff(spks);
        nspks = length(spks);

        % get indices of first and last spike in each burst and force minSpksBrst
        [b.idx, nbrsts] = binary2epochs('vec', isi <= isiThr,...
            'minDur', minSpks - 1, 'printFlag', false);

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
        b.spkprct = sum(b.nspks) / nspks * 100;
        
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if saveVar
    brst.info.runtime = datetime(now, 'ConvertFrom', 'datenum');
    brst.info.input = p.Results;
    brst.info.bins = bins;
    brst.info.recLen = recLen;

    save(brstfile, 'brst')
end

end

% EOF