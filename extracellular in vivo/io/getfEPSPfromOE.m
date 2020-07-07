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
%   spkgrp      cell array. each cell represents a spkgrprode and contains 
%               the channels in that spkgrprode
%   intens      vec describing stimulus intensity [uA]. must be equal in
%               length to number of recording files in experiment. 
%   win         vec of 2 elements. determines length of snip. for example,
%               win = [5 405] than each snip will be 401 samples, starting
%               5 samples after the corresponding stamp and ending 405
%               samples after stamp. if win = [-16 16] than snip will be of
%               33 samples symmetrically centered around stamp.
%   dt          numeric. dead time for calculating amplitude. 
%               important for exclusion of stim artifact. {40}[samples]
%   precision   char. sample precision of dat file {'int16'}
%   force       logical. force reload {false}.
%   concat      logical. concatenate blocks (true) or not {false}. 
%               used for e.g stability.
%   saveVar     logical. save variable {1}. 
%   saveFig     logical. save graphics {1}. 
%   graphics    logical. plot graphics {1}. 
%   vis         char. figure visible {'on'} or not ('off')
%   
% CALLS
%   snipFromDat
%   
% OUTPUT
%   fepsp       struct
% 
% TO DO LIST
%   # more efficient way to convert tstamps to idx
%   # add concatenation option for stability
%   # improve graphics
%   # add option to resample (see getAcc)
%   # sort mats according to intensity
% 
% 22 apr 20 LH   UPDATES:
% 28 jun 20 LH      first average electrodes and then calculate range
%                   dead time to exclude stim artifact

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
addOptional(p, 'dt', 40, @isnumeric);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'force', false, @islogical);
addOptional(p, 'concat', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'saveFig', true, @islogical);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'vis', 'on', @ischar);

parse(p, varargin{:})
basepath = p.Results.basepath;
fname = p.Results.fname;
nchans = p.Results.nchans;
spkgrp = p.Results.spkgrp;
intens = p.Results.intens;
win = p.Results.win;
dt = p.Results.dt;
precision = p.Results.precision;
force = p.Results.force;
concat = p.Results.concat;
saveVar = p.Results.saveVar;
saveFig = p.Results.saveFig;
graphics = p.Results.graphics;
vis = p.Results.vis;

% params
if isempty(spkgrp)
    spkgrp = num2cell(1 : nchans, 2);
end
nspkgrp = length(spkgrp);

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
    fprintf('\nloading %s \n', stimname)
    load(stimname)
else
    error('%s not found', stimname)
end

% load dat info
infoname = fullfile(basepath, [basename, '.datInfo.mat']);
if exist(infoname, 'file')
    fprintf('loading %s \n', infoname)
    load(infoname)
end
nfiles = length(datInfo.origFile);  % number of intensities 

% load timestamps
load(fullfile(basepath, [basename, '.tstamps.mat']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snip data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert tstamps to idx of samples 
stamps = zeros(1, length(din.data));
for i = 1 : length(din.data)
    stamps(i) = find(tstamps == din.data(i));
end

% snip
snips = snipFromDat('basepath', basepath, 'fname', fname,...
    'stamps', stamps, 'win', win, 'nchans', nchans, 'ch', [],...
    'dtrend', false, 'precision', precision);
snips = snips / 1000;   % uV to mV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rearrange indices according to intesities (files)
csamps = [0 cumsum(datInfo.nsamps)];
maxstim = 1;
stimidx = cell(nfiles, 1);
for i = 1 : nfiles
    stimidx{i} = find(stamps > csamps(i) &...
        stamps <  csamps(i + 1));
    maxstim = max([maxstim, length(stimidx{i})]);   % used for cell2mat
end

% rearrange snips according to intensities and tetrodes. extract amplitude
% and waveform.
for j = 1 : nspkgrp
    for i = 1 : nfiles
        wv{j, i} = snips(spkgrp{j}, :, stimidx{i});
        wvavg(j, i, :) = mean(mean(wv{j, i}, 3), 1);
        wvamp = mean(wv{j, i}(:, dt : end, :), 1);
        ampcell{j, i} = squeeze(range(wvamp));
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
% arrange struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[intens, ia] = sort(intens);

% arrange struct
fepsp.wv = wv(:, ia);
fepsp.wvavg = wvavg(:, ia, :);
fepsp.ampcell = ampcell(:, ia);
fepsp.amp = amp(:, ia);
fepsp.stim = stimidx(:, ia);
fepsp.intens = intens;
fepsp.t = win(1) : win(2);
fepsp.spkgrp = spkgrp;
fepsp.dt = dt;

if saveVar
    save(fepspname, 'fepsp');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    if saveFig
        figpath = fullfile(basepath, 'graphics');
        mkdir(figpath)
    end
    for i = 1 : nspkgrp
        fh = figure('Visible', vis);
        suptitle(sprintf('T#%d', i))
        subplot(1, 2, 1)
        plot(squeeze(fepsp.wvavg(i, :, :))')
        axis tight
        y = [min([fepsp.wvavg(:)]) max([fepsp.wvavg(:)])];
        hold on
        plot([dt dt], y, '--r')
        ylim(y)
        xlabel('Time [samples]')
        ylabel('Voltage [mV]')
        legend(split(num2str(intens)))
        box off
        
        subplot(1, 2, 2)
        ampmat = cell2nanmat(fepsp.ampcell(i, :));
        boxplot(ampmat, 'PlotStyle', 'traditional')
        ylim([min(vertcat(fepsp.ampcell{:})) max(vertcat(fepsp.ampcell{:}))])
        xticklabels(split(num2str(intens)))
        xlabel('Intensity [uA]')
        ylabel('Amplidute [mV]')
        box off
        
        if saveFig
            figname = [figpath '\fepsp_t' num2str(i)];
            export_fig(figname, '-tif', '-r300', '-transparent')
        end
    end
end

end

% EOF