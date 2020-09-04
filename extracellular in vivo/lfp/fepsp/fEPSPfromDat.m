function fepsp = fEPSPfromDat(varargin)

% this is a wrapper to get fEPSP signals snipped from .dat or .lfp files
% based on stim indices. compatible for both OE and TDT. if OE, assumes
% preprocOE has been called beforehand and that basepath conatins din.mat.
% if TDT, assumes tdt2dat was called beforehand and basepath is the data
% tank with all relavent blocks.
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
%   dur         numeric. duration of snip {0.15}[s]
%   dt          numeric. dead time for calculating amplitude.
%               important for exclusion of stim artifact. {2}[ms]
%   precision   char. sample precision of dat file {'int16'}
%   extension   string. load from 'dat' or {'lfp'} file
%   recSystem   string. data from 'tdt' or {'oe'}
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
%   tdtbin2mat
%
% OUTPUT
%   fepsp       struct
%
% TO DO LIST
%   # more efficient way to convert tstamps to idx
%   # add concatenation option for stability
%   # improve graphics
%   # add option to resample (see getAcc) (done)
%   # sort mats according to intensity (done)
%
% 22 apr 20 LH   UPDATES:
% 28 jun 20 LH      first average electrodes and then calculate range
%                   dead time to exclude stim artifact
% 03 sep 20 LH      snip from lfp
% 04 sep 20 LH      tdt and oe

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'intens', [], @isnumeric);
addOptional(p, 'dur', 0.15, @isnumeric);
addOptional(p, 'dt', 2, @isnumeric);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'extension', 'lfp', @ischar);
addOptional(p, 'recSystem', 'oe', @ischar);
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
dur = p.Results.dur;
dt = p.Results.dt;
precision = p.Results.precision;
extension = p.Results.extension;
recSystem = p.Results.recSystem;
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

switch recSystem
    case 'oe'
        b2uv = 0.195;
    case 'tdt'
        b2uv = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure dat file exist
datfiles = dir([basepath filesep '**' filesep '*' extension]);
if isempty(datfiles)
    error('no .dat files found in %s', basepath)
end
if isempty(fname)
    if length(datfiles) == 1
        fname = datfiles.name;
    else
        fname = [bz_BasenameFromBasepath(basepath) '.' extension];
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

% load dat info
infofile = dir('*datInfo*');
infoname = infofile.name;
if exist(infoname, 'file')
    fprintf('loading %s \n', infoname)
    load(infoname)
end
nfiles = length(datInfo.nsamps);  % number of intensities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get stim indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch recSystem
    case 'oe'
        fsIn = 20000;
        
        % load digital input
        stimname = fullfile(basepath, [basename, '.din.mat']);
        if exist(stimname)
            fprintf('\nloading %s \n', stimname)
            load(stimname)
        else
            error('%s not found', stimname)
        end
        
        % load timestamps
        load(fullfile(basepath, [basename, '.tstamps.mat']));
        
        % convert tstamps to idx of samples
        stamps = zeros(1, length(din.data));
        for i = 1 : length(din.data)
            stamps(i) = find(tstamps == din.data(i));
        end
        
    case 'tdt'
        % for some blocks the length of two stores (e.g. Raw1 and stim) defer even
        % if they were sampled at the same frequency. this appears random (i.e.
        % does not depend on block duration and can go either way. the workaround
        % here is to shorten / elongate stim according to Raw1. this assumes the
        % missing / additional samples are at the end of the block and are not
        % important (so far has proven to be true)
        fsIn = 24414.0625;
        blockfiles = dir('block*');
        blocknames = natsort({blockfiles.name});
        stim = [];
        for i = 1 : length(datInfo.blocks)
            blockpath = fullfile(basepath, datInfo.blocks{i});
            store = datInfo.store;
            dat = TDTbin2mat(blockpath, 'TYPE', {'streams'},...
                'STORE', store, 'T1', 0, 'T2', 0);
            dat = dat.streams.(store).data;
            store = 'Stim';
            s = TDTbin2mat(blockpath, 'TYPE', {'streams'},...
                'STORE', store, 'T1', 0, 'T2', 0);
            s = s.streams.(store).data;
            
            d = length(s) - length(dat);
            if d ~= 0
                warning('length stim and raw not equal')
                if d > 0
                    s(end : -1 : end - d + 1) = [];
                else
                    s(length(s) : length(dat)) = 0;
                end
            end
            stim = [stim, s];
        end   
        % stim onset from diff
        stamps = find(diff(stim) > max(diff(stim)) / 2);
end

if strcmp(extension, 'lfp')
    fsRatio = fsIn / 1250;     % dat / lfp
    stamps = round(stamps / fsRatio);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snip data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set window of snip (150 ms as in wcp).
% assumes dat recorded at 20 kHz and lfp at 1.25 kHz.
if strcmp(extension, 'lfp')
    fsOut = 1250;
else
    fsOut = fsIn;
end
win = round([1 dur * fsOut]);
dt = round(dt / 1000 * fsOut);

% snip
snips = snipFromDat('basepath', basepath, 'fname', fname,...
    'stamps', stamps, 'win', win, 'nchans', nchans, 'ch', [],...
    'dtrend', false, 'precision', precision, 'extension', extension,...
    'b2uv', b2uv);
snips = snips / 1000;   % uV to mV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rearrange indices according to intesities (files)
switch extension
    case 'lfp'
        nsamps = datInfo.nsamps / fsRatio;
    case 'dat'
        nsamps = datInfo.nsamps;
end
csamps = [0 cumsum(nsamps)];
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