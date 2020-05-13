function spikes = getSPKfromDat(varargin)

% neurosuite and other utilities high pass filter the data before
% extracting spikes but this typically distorts the waveform. here, spikes
% are directly extracted from the raw dat file and detrended (i.e. the mean
% is substracted and the slope removed). this eliminates the effects of,
% e.g., a slow oscillation "underneath" the spike.
% memmap. Thus, memmap was embedded within this function. 
% -------------------------------------------------------------------------
% 13 may 20 - compared detrending each spike and averaging vs. detrending
% the average waveform. for 10 units, detrending each spike took 240 s vs
% 78 s for detrending the average. The maximum difference was ~1e-13.
% Still, it seems to me risky to detrend the average thus the default will
% be to detrend each spike separately.
% -------------------------------------------------------------------------
% originally, the intention was to use snipFromDat.m. However, this
% recreates a memory map for each unit. Thus, for 69 units and 2608583
% spikes this took 24.5 minutes w/ snipFromDat and 42.5 min with one time
%
% INPUT:
%   basepath    string. path to .dat file (not including dat file itself)
%               {pwd}.
%   fname       string. name of dat file. can be empty if only one dat file
%               exists in basepath or if fname can be extracted from basepath
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
%
% 13 may 20 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'spktimes', {}, @iscell);
addOptional(p, 'win', [-20 19], @isnumeric);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'grp', [], @isnumeric);
addOptional(p, 'ch', [], @isnumeric);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'b2uv', 0.195, @isnumeric);
addOptional(p, 'dtrend', false, @islogical);
addOptional(p, 'l2norm', false, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
fname = p.Results.fname;
spktimes = p.Results.spktimes;
win = p.Results.win;
nchans = p.Results.nchans;
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
spkname = [basename '.spikes.mat'];

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
fprintf('extracting waveforms from %s\n\n', fname)

% memory map to dat file
info = dir(fname);
nsamps = info.bytes / nbytes / nchans;
m = memmapfile(fname, 'Format', {precision, [nchans, nsamps] 'mapped'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get waveforms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 11 : 25
    
    stamps = round(spikes.times{j} * spikes.samplingRate);
       
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
        
        % detrend
        if dtrend
            v = [detrend(v')]';
        end
        
        snips(:, :, i) = v;
    end

    avgwv(j, :, :) = [detrend(mean(snips, 3)')]';    
    stdwv(j, :, :) = std(snips, [], 3);    
    x = squeeze(avgwv(j, :, :));
    [~, maxi] = max(abs(min(x, [], 2) - max(x, [], 2)));
    maxch(j) = ch{grp(j)}(maxi);
    maxwv(j, :) = avgwv(j, maxi, :);
    
end

spikes.avgwv = avgwv;
spikes.stdwv = stdwv;
spikes.maxwv = maxwv;
spikes.maxch = maxch;

save([basepath filesep basename '.spikes.mat'], 'spikes')

fprintf('\nthat took %d minutes\n', toc / 60)


figure
subplot(1, 2, 1)
plotWaveform('avgwv', avgwv, 'sbar', false)
subplot(1, 2, 2)
plotWaveform('avgwv', spikes.avgWaveform{i}, 'sbar', false)


