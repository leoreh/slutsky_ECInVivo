function lfp = getfEPSPfromTDT(varargin)

% this is a wrapper to get lfp signals and stim from tdt
%  
% INPUT
%   basename    string. filename of lfp file. if empty retrieved from
%               basepath. if .lfp should not include extension, if .wcp
%               should include extension
%   basepath    string. path to load filename and save output {pwd}
%   extension   load from {'lfp'} (neurosuite), 'abf', or 'wcp'.
%   forceL      logical. force reload {false}.
%   fs          numeric. requested sampling frequency {1250}
%   mapch       new order of channels {[]}.
%   ch          numeric. channel for lfp signal. if more than one
%               specificed then signal will be the mean of the selected
%               channels. select according to new order (mapch).
%   clip        array of mats indicating times to diregard from recording.
%               each cell corresponds to a block. for example:
%               clip{3} = [0 50; 700 Inf] will remove the first 50 s from
%               block-3 and the time between 700 s and the end of Block-3
%   fdur        duration of fepsp waveform after stim {0.2} [s]
%   saveVar     save variable {1}. 

%   
% OUTPUT
%   lfp         structure with the following fields:
%   fs
%   fs_orig
%   extension
%   interval    
%   duration    
%   chans
%   timestamps 
%   data  
% 
% TO DO LIST
%   # compare resampling output with that of ndmanager plugin
%   # for lh49_200327, the length of 'Stim' stream was shorter than 'Raw1'
%   when loaded from tdtbin2mat
% 
% 27 mar 20 LH          


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'andname', '');
addOptional(p, 'blocks', 1, @isnumeric);
addOptional(p, 'mapch', [1 : 16], @isnumeric);
addOptional(p, 'ch', [1 : 16], @isnumeric);
addOptional(p, 'clip', {}, @iscell);
addOptional(p, 'fdur', 0.2, @isnumeric);
addOptional(p, 'fs', 1250, @isnumeric);
addOptional(p, 'saveVar', true, @islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
andname = p.Results.andname;
blocks = p.Results.blocks;
mapch = p.Results.mapch;
ch = p.Results.ch;
clip = p.Results.clip;
fdur = p.Results.fdur;
fs = p.Results.fs;
saveVar = p.Results.saveVar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % load stim
% [info, stim] = tdt2dat('basepath', basepath, 'store', 'Stim',...
%     'blocks', blocks, 'chunksize', [], 'mapch', [], 'rmvch', [],...
%     'clip', clip, 'saveVar', false);
% 
% load lfp
[info, raw] = tdt2dat('basepath', basepath, 'store', 'Raw1',...
    'blocks', blocks, 'chunksize', [], 'mapch', mapch, 'rmvch', [],...
    'clip', clip, 'saveVar', false);

% workaround for f***ing tdt
blockfiles = dir('block*');
blocknames = {blockfiles.name};
blocknames = natsort(blocknames);
stim = [];
raw = [];
for i = 1 : length(blocks)
    blockpath = fullfile(basepath, blocknames{blocks(i)});
    store = 'Raw1';
    dat = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store,...
        'T1', 0, 'T2', 0);
    dat = dat.streams.(store).data;
    store = 'Stim';
    s = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store,...
        'T1', 0, 'T2', 0);
    s = s.streams.(store).data;
    
    d = length(s) - length(dat);
    if d ~= 0
        warning('length stim and raw not equal')
        if d > 0
            s(end : -1 : end - d + 1) = [];
        else
            dat(:, end : -1 : end - abs(d) + 1) = [];
        end
    end
    stim = [stim, s];
    raw = [raw, dat];
end

% find stim onset from diff
stim = find(diff(stim) > max(diff(stim)) / 2);
% convert stim onset times to new fs
stim = round(stim / info.fs * fs);

% resample
[p, q] = rat(fs / info.fs);
n = 5; beta = 20;
raw = [resample(double(raw)', p, q, n, beta)]';

% remove dc, average and convert to mV
raw = [rmDC(raw')]';
raw = mean(raw(ch, :)) / 1000;

% clip fepsp
for i = 1 : length(stim)
    fepsp(i, :) = raw(stim(i) : stim(i) + fdur * fs);
end

% calc amp and arrange in matrix according to blocks
stimidx = stim / fs;
cumdur = [0 cumsum(info.blockduration)];
maxstim = 1;
for i = 1 : length(blocks)
    idx{i} = find(stimidx > cumdur(i) & stimidx < cumdur(i + 1));
    maxstim = max([maxstim, length(idx{i})]);
end
for i = 1 : length(blocks)
    wv{i} = fepsp(idx{i}, :);
    wvavg(i, :) = [mean(wv{i})]';
    amp{i} = abs(min(wv{i}, [], 2) - max(wv{i}, [], 2));
    stimcell{i} = stim(idx{i});
end
mat = cellfun(@(x)[x(:); NaN(maxstim - length(x), 1)], amp,...
    'UniformOutput', false);
amp = cell2mat(mat);
mat = cellfun(@(x)[x(:); NaN(maxstim - length(x), 1)], stimcell,...
    'UniformOutput', false);
stim = cell2mat(mat);

% arrange struct
lfp.wv = wv;
lfp.wvavg = wvavg';
lfp.stim = stim;
lfp.t = 0 : 1 / fs : fdur;
lfp.amp = amp;
lfp.fs = fs;
lfp.fs_orig = info.fs;
lfp.blocks = blocks;
lfp.blockduration = info.blockduration;
lfp.mapch = mapch;
lfp.ch = ch;
lfp.clip = clip;
lfp.filename = [info.filename '.' andname '.mat'];

% save variable
if saveVar   
    save(lfp.filename, 'lfp')
end

end
