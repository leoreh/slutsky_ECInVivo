function brst = spktimes_meaBrst(spktimes, varargin) 

% calculates burstiness stats per unit in timebins. based on analysis done
% for mea recordings
%
% INPUT
%   basepath    char. fullpath to recording folder {pwd}
%   spktimes    cell of spike times per unit. typically in [s]. can also be
%               samples, ms, etc. but then binsize, and isiThr must be the
%               same units
%   binsize     numeric. duration of bins {3600 [s]}
%   bins        cell array of n x 2 mats of intervals.
%               metrices will be calculated for each cell by limiting
%               spktimes to the intervals. can be for example
%               ss.boutTimes. must be the same units as spikes.times
%               (e.g. [s]). will override binsize
%   isiThr      numeric. minimum inter-spike interval for defining a burst
%               {0.02 [s]}
%   minSpks     numeric. minimum number of spikes for defining a burst {2}
%   
%   flg_all     logical. save data of each burst (true) or just mean across
%               bursts {fale}
%   flg_save     logical. save struct {true}
%   flg_force   logical. flg_force analyze even if file exists {false}
% 
% OUTPUT
%   brst        struct
%
% CALLS
%   binary2bouts
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
addParameter(p, 'flg_save', true, @islogical)
addParameter(p, 'flg_all', false, @islogical)
addParameter(p, 'flg_force', false, @islogical)

parse(p, varargin{:})
basepath        = p.Results.basepath;
isiThr          = p.Results.isiThr;
binsize         = p.Results.binsize;
bins            = p.Results.bins;
minSpks         = p.Results.minSpks;
flg_save         = p.Results.flg_save;
flg_all         = p.Results.flg_all;
flg_force       = p.Results.flg_force;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% separate recording to bins
recLen = max(cellfun(@max, spktimes, 'uni', true));
if isempty(binsize) || binsize > recLen
    binsize = recLen;
end
if isempty(bins)
    bins = n2chunks('n', recLen, 'chunksize', binsize, 'overlap', [1 0]);
    bins(1) = 0; bins(end - 1, 2) = recLen; bins(end, :) = [];
end

% convert bins to cell
if ~iscell(bins)
    bins = mat2cell(bins, ones(size(bins, 1), 1), size(bins, 2));
end
nbins = length(bins);

for ibin = 1 : nbins
    if bins{ibin}(end) > recLen
        bins{ibin}(end) = recLen;
    end
end

% get binsize
for ibin = 1 : nbins
    binsize(ibin) = sum(bins{ibin}(:, 2) - bins{ibin}(:, 1));
end

% load if exists
[~, basename] = fileparts(basepath);
brstfile = fullfile(basepath, [basename, '.st_brst.mat']);
if exist(brstfile, 'file') && ~flg_force
    load(brstfile, 'brst')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc burstiness per unit per timebin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nunits = length(spktimes);

% initialize
clear brst
brst.info.runtime       = datetime("now");
brst.info.input         = p.Results;
brst.info.bins          = bins;
brst.info.binsize       = binsize;
brst.info.recLen        = recLen;
brst.detect             = nan(nbins, nunits);
brst.nspks              = nan(nbins, nunits);
brst.spkprct            = nan(nbins, nunits);
brst.brstDur            = nan(nbins, nunits);
brst.freq               = nan(nbins, nunits);
brst.rate               = nan(nbins, nunits);
brst.rateNorm           = nan(nbins, nunits);
brst.ibi                = nan(nbins, nunits);
brst.shortPrct          = nan(nbins, nunits);

for ibin = 1 : nbins
    for iunit = 1 : nunits
        
        spks = spktimes{iunit}(InIntervals(spktimes{iunit}, bins{ibin}));
        isi = diff(spks);
        nspks = length(spks);

        % get indices of first and last spike in each burst and flg_force minSpksBrst
        [b.idx, nbrsts] = binary2bouts('vec', isi <= isiThr,...
            'minDur', minSpks - 1, 'flgPrnt', false);

        if isempty(b.idx)
            continue
        end

        % number of spikes per burst
        b.nspks = diff(b.idx, 1, 2) + 1;
        
        % start and stop times per burst
        b.times = spks(b.idx(:, 1));
        b.times(:, 2) = spks(b.idx(:, 2));
        
        % burst duration [s]
        b.dur = diff(b.times, 1, 2); 
        
        % frequency within burst
        b.freq = b.nspks ./ (b.dur);
        
        % percent of spikes that participate in bursts
        b.spkprct = sum(b.nspks) / nspks * 100;
        
        % inter-burst interval
        b.ibi = b.times(2 : end, 1) - b.times(1 : end - 1, 2);
        
        % percent bursts that pass the min spks criterion
        brst.shortPrct(ibin, iunit) =  (nbrsts.dur / nbrsts.detect) * 100;

        % number of bursts relative to binsize and ~firing rate
        brst.rate(ibin, iunit) = nbrsts.dur / binsize(ibin);
        brst.rateNorm(ibin, iunit) = brst.rate(ibin, iunit) / (nspks / binsize(ibin));

        % save all burst data
        if flg_all
            brst.all(ibin, iunit) = b;
        else
            brst.all = [];
        end

        % organize and average stats
        brst.detect(ibin, iunit) = nbrsts.detect;
        brst.nspks(ibin, iunit) = mean(b.nspks);
        brst.brstDur(ibin, iunit) = mean(b.dur);
        brst.freq(ibin, iunit) = mean(b.freq);
        brst.spkprct(ibin, iunit) = mean(b.spkprct);
        brst.ibi(ibin, iunit) = mean(b.ibi);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flg_save
    save(brstfile, 'brst')
end

end

% EOF