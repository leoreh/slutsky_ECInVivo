function fepsp = getfEPSPfromOE(varargin)

% this is a wrapper to get fEPSP signals from OE. Assumes preprocOE has
% been called beforehand and that basepath contains both the raw .dat file
% and din.mat.
%  
% INPUT
%   basepath    string. path to .dat file (not including dat file itself)
%   fname       string. name of dat file. if empty and more than one dat in
%               path, will be extracted from basepath
%   nchans      numeric. original number of channels in dat file {35}.
%   spkgrp         cell array. each cell represents a spkgrprode and contains 
%               the channels in that spkgrprode
%   intens      vec describing stimulus intensity [uA]. must be equal in
%               length to number of recording files in experiment. 
%   win         vec of 2 elements. determines length of snip. for example,
%               win = [5 405] than each snip will be 401 samples, starting
%               5 samples after the corresponding stamp and ending 405
%               samples after stamp. if win = [-16 16] than snip will be of
%               33 samples symmetrically centered around stamp.
%   precision   char. sample precision of dat file {'int16'}
%   force       logical. force reload {false}.
%   concat      logical. concatenate blocks (true) or not {false}. 
%               used for e.g stability.
%   saveVar     logical. save variable {1}. 
%   graphics    logical. plot graphics {1}. 
%   
% CALLS
%   snipFromDat
%   
% OUTPUT
%   fepsp       struct
% 
% TO DO LIST
%   # code more efficient way to convert tstamps to idx
%   # add concatenation option for stability
%   # improve graphics
%   # add option to resample (see getAcc)
% 
% 22 apr 20 LH          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'intens', [], @isnumeric);
addOptional(p, 'win', [1 2000], @isnumeric);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'force', false, @islogical);
addOptional(p, 'concat', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
fname = p.Results.fname;
nchans = p.Results.nchans;
spkgrp = p.Results.spkgrp;
intens = p.Results.intens;
win = p.Results.win;
precision = p.Results.precision;
force = p.Results.force;
concat = p.Results.concat;
saveVar = p.Results.saveVar;
graphics = p.Results.graphics;

% params
if isempty(spkgrp)
    spkgrp = num2cell(1 : nchans, 2);
end
nspkgrp = size(spkgrp, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure dat and stim files exist
datfiles = dir([basepath filesep '**' filesep '*dat']);
if isempty(datfiles)
    error('no .dat files found in %s', basepath)
end
if isempty(fname)
    if length(datfiles) == 1
        fname = datfiles.name;
    else
        fname = [bz_BasenameFromBasepath(basepath) '.dat'];
        if ~contains({datfiles.name}, fname)
            error('please specify which dat file to process')
        end
    end
end
[~, basename, ~] = fileparts(fname);

% load fepsp if already exists
fepspname = [basename '.fepsp.mat'];
if exist(fepspname, 'file') && ~force
    fprintf('\n loading %s \n', fepspname)
    load(fepspname)
    return
end

% load digital input
stimname = fullfile(basepath, [basename, '.din.mat']);
if exist(stimname) 
    fprintf('\n loading %s \n', stimname)
    load(stimname)
else
    error('%s not found', stimname)
end

% load dat info
infoname = fullfile(basepath, [basename, '.datInfo.mat']);
if exist(infoname, 'file')
    fprintf('\n loading %s \n', infoname)
    load(infoname)
end
nfiles = length(datInfo.origFile);  % number of intensities 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snip data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert tstamps to idx of samples 
for i = 1 : length(din.data)
    stamps(i) = find(datInfo.tstamps == din.data(i));
end

% snip
snips = snipFromDat('basepath', basepath, 'fname', fname,...
    'stamps', stamps, 'win', win, 'nchans', nchans, 'ch', [],...
    'dtrend', false, 'precision', precision);
snips = snips / 1000;   % uV to mV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find indices according to intesities (different files)
csamps = [0 cumsum(datInfo.nsamps)];
maxstim = 1;
stimidx = cell(nfiles, 1);
for i = 1 : nfiles
    stimidx{i} = find(stamps > csamps(i) &...
        stamps <  csamps(i + 1));
    maxstim = max([maxstim, length(stimidx{i})]);   % used for cell2mat
end

% rearrange snips according to intensities and spkgrprodes
% extract amplitude and waveform
for j = 1 : nspkgrp
    for i = 1 : nfiles
        wv{j, i} = snips(spkgrp{j}, :, stimidx{i});
        wvavg(j, i, :) = mean(mean(wv{j, i}, 3), 1);
        ampcell{j, i} = mean(squeeze(abs(min(wv{j, i}, [], 2) - max(wv{j, i}, [], 2))), 1);
        amp(j, i) = mean(ampcell{j, i});
        stimcell{i} = stamps(stimidx{i});
    end
end

if concat 
    for i = 1 : length(blocks)
        amp = [amp; ampcell{i}];
    end
else
    mat = cellfun(@(x)[x(:); NaN(maxstim - length(x), 1)], stimcell,...
        'UniformOutput', false);
    stimidx = cell2mat(mat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    fh = figure;
    k = 1;
    for j = 1 : nspkgrp
        for i = 1 : nfiles
            subplot(nspkgrp, nfiles, k)
            plot(squeeze(wvavg(j, i, :)))
            ylim([min(wvavg(:)) max(wvavg(:))])
            set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
                'Color', 'none')
            if j == nspkgrp
                xlabel([num2str(intens(i)) ' uA'])
            end
            if i == 1
                set(gca, 'YTickLabel', [ceil(min(wvavg(:))), floor(max(wvavg(:)))])
                ylabel(['T' num2str(j)])
            end
            k = k + 1;
        end
    end
    
    fh = figure;
    k = 1;
    for j = 1 : nspkgrp
        subplot(nspkgrp, 1, k)
        plot(amp(j, :))
        k = k + 1;
        ylim([min(amp(:)) max(amp(:))])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct
fepsp.wv = wv;
fepsp.wvavg = wvavg;
fepsp.ampcell = ampcell;
fepsp.amp = amp;
fepsp.stim = stimidx;
fepsp.intens = intens;
fepsp.t = win(1) : win(2);
fepsp.spkgrp = spkgrp;

if saveVar
    save(fepspname, 'fepsp');
end

end
