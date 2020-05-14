function spikes = getSPKfromDat(varargin)

% neurosuite and other utilities high pass filter the data before
% extracting spikes but this typically distorts the waveform. here, spikes
% are directly extracted from the raw dat file and detrended (i.e. the mean
% is substracted and the slope removed). this eliminates the effects of,
% e.g., a slow oscillation "underneath" the spike. another option in this
% function is to divide each waveform with its norm (L2) such that the sum
% of squares of each waveform is one. This effectively neglects the
% amplitude and allows, e.g., for spikes in a burst to be classified
% together. However, since this function is designed to extract spikes
% after sorting I do not see any reason to implement this.
% -------------------------------------------------------------------------
% originally, the intention was to use snipFromDat.m. However, this
% recreates a memory map for each unit. Thus, for 69 units and 2608583
% spikes this took 24.5 minutes w/ snipFromDat and 42.5 min with one time
% memmap. Thus, memmap was embedded within this function.
% -------------------------------------------------------------------------
% since detrend is essentialy matrix multiplication with a regressor that
% only depends on the waveform length, the matrices (w, r, c) are
% calculated once in this function and used for all spikes. this
% significantly decreased computation time. Note I did not thoroughly
% understand the mathematics here but rather "stole" the lines from matlab
% detrend.m
% -------------------------------------------------------------------------
% compared detrending each spike and averaging vs. detrending the average
% waveform. for 10 units, detrending each spike took 240 s vs 78 s for
% detrending the average. The maximum difference was ~1e-13. Still, it
% seems to me risky to detrend the average thus the default will be to
% detrend each spike separately.
% -------------------------------------------------------------------------
%
% INPUT:
%   basepath    string. path to .dat file (not including dat file itself)
%               {pwd}.
%   fname       string. name of dat file. can be empty if only one dat file
%               exists in basepath or if fname can be extracted from basepath
%   fs          sampling frequency
%   spktimes    a cell array of vectors. each vector (unit) contains the
%               timestamps of spikes for that unit. 
%   win         vec of 2 elements {[-20 19]}. determines length of snip.
%               for example, win = [5 405] than each snip will be 401
%               samples, starting 5 samples after the corresponding stamp
%               and ending 405 samples after stamp. if win = [-16 16] than
%               snip will be of 33 samples symmetrically centered around
%               stamp.
%   nchans      numeric. number of channels in dat file {35}
%   grp         vec. spike group to which each cell belongs. 
%   ch          cell array. each cell contains the channels of a spk group.
%   precision   char. sample precision of dat file {'int16'}
%   b2uv        numeric. conversion of bits to uV {0.195}
%   dtrend      logical. detrend snippets {false}.
%   l2norm      logical. L2 normalize snippets {false}.
%
% OUTPUT
%   spikes      same struct as getSpikes with the following fields:
%       avgwv   mat (unit x ch x sample) of avg waveform .
%       stdwv   mat (unit x ch x sample) of std waveform .
%       maxwv   mat (unit x sample) of waveform on channel with maximum amp 
%       maxch   vec of channel with maximum amplitude for each unit
% 
% CALLS:
%   class2bytes
%
% TO DO LIST:   
%       # fix output for when each group has a different number of channels
%       (done)
%       # calculate trend based on x samples before and after the spike
%       (not including the waveform)
%
% 13 may 20 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'fs', 20000, @isnumeric);
addOptional(p, 'spktimes', {}, @iscell);
addOptional(p, 'win', [-20 19], @isnumeric);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'grp', [], @isnumeric);
addOptional(p, 'ch', {}, @iscell);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'b2uv', 0.195, @isnumeric);
addOptional(p, 'dtrend', false, @islogical);
addOptional(p, 'l2norm', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
fname = p.Results.fname;
fs = p.Results.fs;
spktimes = p.Results.spktimes;
win = p.Results.win;
nchans = p.Results.nchans;
grp = p.Results.grp;
ch = p.Results.ch;
precision = p.Results.precision;
b2uv = p.Results.b2uv;
dtrend = p.Results.dtrend;
l2norm = p.Results.l2norm;

% size of one data point in bytes
nbytes = class2bytes(precision);

if isempty(ch)
    ch = {1 : nchans};
end
if isempty(grp)
    grp = ones(length(spktimes), 1);
end

nunits = length(spktimes);
basename = bz_BasenameFromBasepath(basepath);

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle dat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
cd(basepath)
datFiles = dir([basepath filesep '**' filesep '*dat']);
if isempty(datFiles)
    error('no .dat files found in %s', basepath)
end
if isempty(fname)
    if length(datFiles) == 1
        fname = datFiles.name;
    else
        fname = [bz_BasenameFromBasepath(basepath) '.dat'];
        if ~contains({datFiles.name}, fname)
            error('please specify which dat file to process')
        end
    end
end
fprintf('\nextracting waveforms from %s\n\n', fname)

% memory map to dat file
info = dir(fname);
nsamps = info.bytes / nbytes / nchans;
m = memmapfile(fname, 'Format', {precision, [nchans, nsamps] 'mapped'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build regressor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
N = length(win(1) : win(2));
s = 0 : N - 1;
scaleS = s(end);
a = s./scaleS;
b = max(a, 0);
W = b(:);
W = [reshape(W, N, []), ones(N,1)];
[Q,R] = qr(W,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get waveforms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1 : nunits
    
    stamps = round(spktimes{j} * fs);
       
    % initialize
    sniplength = diff(win) + 1;
    nsnips = length(stamps);
    snips = zeros(length(ch), sniplength, nsnips);
    
    % go over stamps and snip data
    for i = 1 : length(stamps)
        if stamps(i) + win(1) < 1 || stamps(i) + win(2) > nsamps
            warning('skipping stamp %d because snip incomplete', i)
            snips(:, :, i) = nan(length(ch), sniplength);
            continue
        end
        v = double(m.Data.mapped(ch{grp(j)},...
            stamps(i) + win(1) : stamps(i) + win(2)));
        
        % convert bits to uV
        if b2uv
            v = double(v) * b2uv;
        end
        
        % L2 normalize
        if l2norm
            v = v ./ vecnorm(v, 2, 2);
        end
        
        % remove best fit
        if dtrend
            v = [v' - W*(R\Q'*v')]';
        end
        
        snips(:, :, i) = v;
    end

    avgwv{j} = mean(snips, 3);    
    stdwv{j} = std(snips, [], 3);    
    x = squeeze(avgwv{j});
    [~, maxi] = max(abs(min(x, [], 2) - max(x, [], 2)));
    maxch(j) = ch{grp(j)}(maxi);
    maxwv(j, :) = avgwv{j}(maxi, :);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spkname = [basename '.spikes.mat'];
load(spkname);

spikes.avgwv = avgwv;
spikes.stdwv = stdwv;
spikes.maxwv = maxwv;
spikes.maxch = maxch;

save([basepath filesep spkname], 'spikes')

fprintf('\nthat took %.1f minutes\n', toc / 60)

end

% EOF
