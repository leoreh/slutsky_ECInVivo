function fepsp = fEPSPfromDat(varargin)

% gets fEPSP signals snipped from .dat or .lfp files based on stim indices.
% compatible for both OE and TDT. if OE, assumes preprocOE has been called
% beforehand and that basepath conatins din.mat. if TDT, assumes tdt2dat
% was called beforehand and basepath is the data tank with all relavent
% blocks.
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
%   protocol    stimulus protocol. can be a string ('io' or 'stp') or a
%               numeric vector describing the times [ms] of stimulations
%               for each trace (currently not implemented)
%   dc          logical. remove DC component {1}
%   cf          numeric. cut-off frequency for low-pass filter {450}
%   precision   char. sample precision of dat file {'int16'}
%   extension   string. load from 'dat' or {'lfp'} file
%   recSystem   string. data from 'tdt' or {'oe'}
%   force       logical. force reload {false}
%   saveVar     logical. save variable {1}
%   inspect     logical. inspect traces {0}
%   anaflag     logical. send to analysis {1}
%   saveFig     logical. save graphics {1}. Only relevant If anaflag is
%               true
% CALLS
%   snipFromDat
%   tdtbin2mat
%   fepsp_analysis
%
% OUTPUT
%   fepsp           struct with fields:
%       Info        struct, info about the recourding, with fields:
%           basename    cell, the name of the file/ files that the traces were
%                       exstracted from.
%           recSystem   char, type of recourding system (wcp/oe/tdt), in
%                       this function oe/tdt.
%           protocol    stimulus protocol. can be a string ('io' or 'stp') or a
%                       numeric vector describing the times [ms] of stimulations
%                       for each trace (numeric currently not implemented)
%           fsOrig      double, original sampling frequency of file.
%           fs          double, sampling frequency after downsampling
%           spkgrp      cell,  which channles are realated to which
%                       electrode. Each column == diffrent electrode.
%           stimIdx     double, the indexes of the stimuli in the long dat
%                       file, dims electrode X intensities.
%           stimTs      double, the index of the stimuli in the TS of the
%                       long dat file.
%           stamps      double, the index of the start and end of each
%                       trace in the long dat file (need to double check meaning for orginization).
%           intensOrig  double, the intensities as user inputed them, without 
%                       any sorting and merging.
%           rm          cell, for each intensity the number of the traces
%                       user removed. Column order matching intensOrig.
%           lowPass     double, the cut off frequency for the low pass filter used.
%           inspect     logical, true if user inspected the traces in each
%                       intensity in order to remove unwanted ones.
%       intens      double the intensities in this file after sorting and merging.
%       tstamps     the time stamps for all traces. 0 is the first stimulus.
%       traces      cell, all the traces extracted, dims electrode X intensities.
%
% TO DO LIST
%   # more efficient way to convert tstamps to idx
%   # add concatenation option for stability
%   # improve graphics (done)
%   # add option to resample (see getAcc) (done)
%   # sort mats according to intensity (done)
%
% 22 apr 20 LH   UPDATES:
% 28 jun 20        first average electrodes and then calculate range
%                  dead time to exclude stim artifact
% 03 sep 20        snip from lfp
% 04 sep 20        tdt and oe
% 16 oct 20        added inspect, dc, and cf
%                  compatible with wcp
%                  separated analysis
% 01 Nov 20 LD     change fepsp.info.basename to cell for compatiblity with
%                  wcp, and add description of base output
% 05 Dec 20 LD     Fliped tstamps vector on fepsp struct to match WCP
%                  Now passing saveVar & saveFig to fEPSP_analysis
%                  When anaflag is true, fepsp output + savedVar will be the
%                  output of fEPSP_analysis
%                  Comment taking intens from Old type existing fepsp file
%                  while force is true
%                  Addif force is false check for analysis before returning
%                  (and send to analyse if needed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'fname', '', @ischar);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'intens', [], @isnumeric);
addOptional(p, 'protocol', 'io', @ischar);
addOptional(p, 'precision', 'int16', @ischar);
addOptional(p, 'extension', 'lfp', @ischar);
addOptional(p, 'recSystem', 'oe', @ischar);
addOptional(p, 'cf', 450, @isnumeric);
addOptional(p, 'dc', true, @islogical);
addOptional(p, 'force', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'anaflag', true, @islogical);
addOptional(p, 'inspect', false, @islogical);
addOptional(p, 'saveFig', true, @islogical);


parse(p, varargin{:})
basepath    = p.Results.basepath;
fname       = p.Results.fname;
nchans      = p.Results.nchans;
spkgrp      = p.Results.spkgrp;
intens      = p.Results.intens;
protocol    = p.Results.protocol;
dc          = p.Results.dc;
cf          = p.Results.cf;
precision   = p.Results.precision;
extension   = p.Results.extension;
recSystem   = p.Results.recSystem;
inspect     = p.Results.inspect;
force       = p.Results.force;
saveVar     = p.Results.saveVar;
anaflag     = p.Results.anaflag;
saveFig     = p.Results.saveFig;

% params
if isempty(spkgrp)
    spkgrp = num2cell(1 : nchans, 2);
end
nspkgrp = length(spkgrp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure dat / lfp file exists
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
if exist(fepspname, 'file')
    fprintf('\nloading %s \n', fepspname)
    load(fepspname)
    if ~force && ~anaflag
        return
    elseif ~force && anaflag
        if isfield(fepsp,'amp')
            return
        else
            fprintf('Loaded File isn''t analysed, analaysing...\n')
            fepsp = fEPSP_analysis('fepsp', fepsp,'saveFig',saveFig,'saveVar',saveVar,'savename',fepspname);
        end
    else
%         if isfield(fepsp, 'origIntens')
            intens = fepsp.info.intensOrig;
%         elseif isfield(fepsp, 'intens')
%             intens = fepsp.intens;
%         end
    end
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
        b2uv = 0.195;

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
        b2uv = [];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% snip data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set window of snip.
% assumes lfp files sampled at 1.25 kHz.
if strcmp(extension, 'lfp')
    fsOut = 1250;
    fsRatio = fsIn / fsOut;     % dat / lfp
    stamps = round(stamps / fsRatio);
else
    fsOut = fsIn;
end

switch protocol
    case 'io'
        % single pulse of 500 us. recording length 150 ms.
        % repeated once every 15 s. negative peak of response typically
        % 10 ms after stim.
        nstim = 1;
        baseline = [1 floor(29 * fsOut / 1000)];        % samples of baseline period
        snipwin = round([-0.03 * fsOut 0.12 * fsOut]);
        tstamps = [snipwin(1) : snipwin(2)] / fsOut * 1000;
        ts = [];
    case 'stp'
        % 5 pulses of 500 us at 50 Hz. recording length 200 ms. repeated
        % once every 30 s.
        nstim = 5;
        baseline = [1 floor(9 * fsOut / 1000)];        % samples of baseline period
        snipwin = round([-0.01 * fsOut 0.19 * fsOut]);
        tstamps = [snipwin(1) : snipwin(2)] / fsOut * 1000;
        ts = diff(stamps);
        stamps = stamps(1 : 5 : end);
end

% snip
snips = snipFromDat('basepath', basepath, 'fname', fname,...
    'stamps', stamps, 'win', snipwin, 'nchans', nchans, 'ch', [],...
    'dtrend', false, 'precision', precision, 'extension', extension,...
    'b2uv', b2uv);
snips = snips / 1000;   % uV to mV

% clean snips - remove DC componenet and low pass filter
import iosr.dsp.*
filtRatio = cf / (fsOut / 2);
for i = 1 : size(snips, 1)
    x = squeeze(snips(i, :, :));
    if dc
        x = rmDC(x, 'dim', 1, 'win', baseline);
    end
    
    if cf && ~strcmp(extension, 'lfp')
        sz = size(x);
        x = x(:);
        x = iosr.dsp.sincFilter(x, filtRatio);
        x = reshape(x, sz);
    end
    snips(i, :, :) = x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rearrange stim indices according to intesities (files)
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

% combine files if belong to same intensity
[~, ind] = unique(intens);
ind = setdiff([1 : nfiles], ind);
if ~isempty(ind)
    for i = 1 : length(ind)
        stimidx{ind(i) - 1} = [stimidx{ind(i) - 1 : ind(i)}];
    end
    stimidx(ind) = [];
    nfiles = length(stimidx);
end

% rearrange snips according to intensities and tetrodes (spkgrp). averages
% across electrodes from the same spkgrp.
traces = cell(nspkgrp, nfiles);
for j = 1 : nspkgrp
    for i = 1 : nfiles
        traces{j, i} = squeeze(mean(snips(spkgrp{j}, :, stimidx{i}), 1));
        stimcell{i} = stamps(stimidx{i});
    end
end
stimidx = cell2nanmat(stimcell);

% correct orientation if only one trace
for i = 1 : nfiles
    for j = 1 : nspkgrp
        if size(traces{j, i}, 1) < size(traces{j, i}, 2)
            traces{j, i} = traces{j, i}';
        end
    end
end

% manually inspect and remove unwanted traces
sg = 1;                 % selected tetrode 
rm = cell(nfiles);
if inspect
    for i = 1 : nfiles
        [~, rm{i}] = rmTraces(traces{sg, i});
        for j = 1 : nspkgrp
            traces{j, i}(:, rm{i}) = [];
        end
    end
end
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create struct and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort according to intensity
intensOrig = intens;
[intens, ia] = unique(intens);

clear fepsp
fepsp.info.basename = {basename};
fepsp.info.recSystem = recSystem;
fepsp.info.protocol = protocol;
fepsp.info.fsOrig = fsIn;
fepsp.info.fs = fsOut;
fepsp.info.spkgrp = spkgrp;
fepsp.info.stimIdx = stimidx(:, ia);
fepsp.info.stimTs = ts;                 % for stp
fepsp.info.stamps = stamps;             % to recover tdt
fepsp.info.intensOrig = intensOrig;
fepsp.info.rm = rm;
fepsp.info.lowPass = cf;
fepsp.info.inspect = logical(inspect);
fepsp.intens = intens;
fepsp.tstamps = tstamps';
fepsp.traces = traces(:, ia);

if saveVar && ~anaflag
    save(fepspname, 'fepsp');
end

% send to analysis
if anaflag
    fepsp = fEPSP_analysis('fepsp', fepsp,'saveFig',saveFig,'saveVar',saveVar,'savename',fepspname);
end

return

% EOF

% EOF