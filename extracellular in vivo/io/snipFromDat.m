function snips = snipFromDat(varargin)

% maps dat file to memory and snips segments sorrounding specific
% samples. user defines the window length sorrounding each snippet (does
% not have to be symmetrical) and can arrange the output according to
% channel groups
%
% INPUT:
%   datpath     string. path to .dat file (not including dat file itself)
%               {pwd}.
%   fname       string. name of dat file. can be empty if only one dat file
%               exists in datpath or if fname can be extracted from datpath
%   stamps      vec. pointers to dat files from which snipping will occur
%               [samples].
%   win         vec of 2 elements. determines length of snip. for example,
%               win = [5 405] than each snip will be 401 samples, starting
%               5 samples after the corresponding stamp and ending 405
%               samples after stamp. if win = [-16 16] than snip will be of
%               33 samples symmetrically centered around stamp.
%   nchans      numeric. number of channels in dat file {35}
%   ch          vec. channels to load from dat file {[]}. if empty than all will
%               be loaded
%   precision   char. sample precision of dat file {'int16'}
%   dtrend      logical. detrend and L2 normalize snippets {false}.
%   saveVar     logical. save snips {true} or not (false).
%
% OUTPUT
%   snips       matrix of ch x sampels x stamps. if detrend is false than
%               will be same precision as dat file. if true than will be
%               double. 
% 
% CALLS:
%   class2bytes
%
% TO DO LIST:
%   understand the porpuse of L2 normalization
%
% 10 apr 20 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'datpath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'stamps', [], @isnumeric);
addOptional(p, 'win', [-16 16], @isnumeric);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'ch', [], @isnumeric);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'dtrend', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
datpath = p.Results.datpath;
fname = p.Results.fname;
stamps = p.Results.stamps;
win = p.Results.win;
nchans = p.Results.nchans;
ch = p.Results.ch;
precision = p.Results.precision;
dtrend = p.Results.dtrend;
saveVar = p.Results.saveVar;

% size of one data point in bytes
nbytes = class2bytes(precision);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle dat file
cd(datpath)
datFiles = dir([datpath filesep '**' filesep '*dat']);
if isempty(datFiles)
    error('no .dat files found in %s', datpath)
end
if isempty(fname)
    if length(datFiles) == 1
        fname = datFiles.name;
    else
        fname = [bz_BasenameFromBasepath(datpath) '.dat'];
        if ~contains({datFiles.name}, fname)
            error('please specify which dat file to process')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% memory map to dat file
info = dir(fname);
nsamps = info.bytes / nbytes / nchans;
m = memmapfile(fname, 'Format', {precision, [nchans, nsamps] 'mapped'});

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
    v = m.Data.mapped(ch, stamps(i) + win(1) : stamps(i) + win(2));
     
    % L2 normalize and detrend
    if dtrend
        v = double(v) ./ vecnorm(double(v), 2, 2);
        v = [detrend(v')]';
    end
    snips(:, :, i) = v;
end

clear m

end

% EOF