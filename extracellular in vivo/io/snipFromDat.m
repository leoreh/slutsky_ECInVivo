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
%   win         vec.  
%   nchans      numeric. number of channels in dat file {35}
%   precision   char. sample precision of dat file {'int16'}
%   ch          vec. channels to load from dat file {[]}. if empty than all will
%               be loaded
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
%   re
%
% 10 apr 20 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'datpath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'chunksize', 5e6, @isnumeric);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'mapch', [], @isnumeric);
addOptional(p, 'rmvch', [], @isnumeric);
addOptional(p, 'pli', 0, @isnumeric);
addOptional(p, 'bkup', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
datpath = p.Results.datpath;
fname = p.Results.fname;
chunksize = p.Results.chunksize;
precision = p.Results.precision;
nchans = p.Results.nchans;
mapch = p.Results.mapch;
rmvch = p.Results.rmvch;
pli = p.Results.pli;
bkup = p.Results.bkup;
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
[~, basename, ~] = fileparts(fname);

info = dir(fname);
nsamps = info.bytes / nbytes / nchans;

datpath = 'E:\Data\Dat\lh50\lh50_200402\190448_e2r1-7';
nchans = 35;
ch = 1 : 4;
stamps = din.data
win = [1 : 0.002 * 20000];

nsamps = length(win) + 1;
nsnips = length(stamps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% memory map to dat file
m = memmapfile(fname, 'Format', {precision, [nchans, nsamps] 'mapped'});

idx = stamps(i) + 3020;

% go over stamps and snip data
snips = zeros(length(ch), nsamps, nsnips);
for i = 1 : length(stamps) / 4
    v = m.Data.mapped(ch, stamps(i) : stamps(i) + win(end));
    
    
    v = m.Data.mapped(ch, idx : idx + win(end));

    % L2 norm (euclidean distance)
    if dtrend
        for j = 1 : length(ch)
            nv = double(v(j, :));
            nv = nv / norm(nv);
            snips(j, :, i) = detrend(nv);
        end
    end
%     figure
%     plot(snip(j, :, i))
%     yyaxis right
%     plot(v(j, :))
%     legend('detrend', 'raw')
end

clear m

end


% EOF