function [bursty] = isbursty(varargin)

% this function defines a cell as bursty if (1) there is a local maximum in
% the isi log distribution between blim {[1.5 8]} and (2) if the ratio
% between isis within blim and total isi is greater than bratio {0.05}.
% based on Medrihan et al., 2017.
%
% INPUT
%   spktimes    a cell array of vectors. each vector (unit) contains the
%
% OUTPUT
%   bursty      logical {1 = bursty}. 
%
% DEPENDENCIES
%   gaussFilt       from TSToolbox
%
% TO DO LIST
%   # implement Buzsaki ACG critirion
%   # use prominance of peak
% 
% 14 may 20 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'spktimes', @iscell);

parse(p, varargin{:})
spktimes = p.Results.spktimes;

nunits = length(spktimes);

% params
bwin = [1.5 8];     % win for peak detection in isi [ms]
bratio = 0.05;      % ratio of observed over expected isi histogram peak.
peaksep = 0.5;      % minimum time [ms] between peaks.
binsize = 0.05;     % binsize for iis distribution [log10(ms)].
isiwin = [-1 : binsize : 4];   % win for calc isi distribution [ms]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
bpeak = nan(nunits, 1);
bursty = zeros(nunits, 1);

for i = 1 : nunits
    s = spktimes{i} * 1000;   % convert to ms
    
    % isi distribution
    isi = diff(s);
    isi = isi(isi < 10 ^ isiwin(end));
    h = histcounts(log10(isi), isiwin(1 : end - 1));
    hs = gaussFilt(h(:), 2, 0) / (length(s) - 1);
    
    % find local maximum in range
    [~, locs, ~, ~] = findpeaks(hs, 'MinPeakDistance', peaksep / binsize);
    locs(isiwin(locs) < log10(bwin(1))) = [];
    
    % critirion for burstiness
    if ~isempty(locs)
        n = sum(isi >= bwin(1) & isi <= bwin(2)) / length(isi);
        bursty(i) = 10 ^ isiwin(locs(1)) < bwin(2) & n > bratio;
        bpeak(i) = 10 ^ isiwin(locs(1));
    end
    
end

end

% EOF
