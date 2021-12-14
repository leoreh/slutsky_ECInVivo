
% gets the field (stimulus) response from continuous lfp recordings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load / save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'E:\Data\Dat\lh47\lh47_200307';
[~, filename] = fileparts(basepath);
load([filename '_io1.mat'])

% remove DC component
data.raw = rmDC(data.raw);
% convert to mV
data.raw = data.raw / 1000;
% save([filename '_io1.mat'], 'data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traceL = 0.2;         % trace length [s]
nblocks = length(data.blockduration);
ch = 4;
inspect = false;
fs = data.fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arrange data in mat according to stimulus response
[~, stimidx] = findpeaks(diff(data.stim), 'MinPeakHeight', 1000);
stimidx = stimidx(1 : 5 : end);
% stimidx = find(diff(data.stim) > 1000);
nstim = length(stimidx);

for i = 1 : nstim
    datamat(:, i, :) = data.raw(stimidx(i) : stimidx(i) + data.fs * traceL, :);
end

% inspect traces
if inspect
    for i = 1 : ch
        [datamat(:, :, i), rm{i}] = rmTraces(squeeze(datamat(:, :, i)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arrange mat according to blocks
dur = cumsum(data.blockduration * fs);
x = 0;
for i = 1 : length(data.blocks)
    idx = find(stimidx < dur(i) & stimidx > x);
    data.stimmat{i} = datamat(:, idx, :);
    x = dur(i);
end

% % rearrange if needed to sort intensity
% blockidx = [1, 7, 8, 2 : 6];
% data.stimmat = data.stimmat(blockidx);

% amplitude and trace for prism
tracet = [1 : size(datamat, 1)] / data.fs;
for i = 1 : length(data.stimmat)
    for j = 1 : ch
        mat = squeeze(data.stimmat{i}(:, :, j));
        amp{j, i} = max(mat, [], 1);
        trace{j}(i, :) = mean(mat, 2);
    end
end

maxtraces = 15;
for i = 1 : ch
    ampmat{i} = cellfun(@(x)[x(:); NaN(maxtraces-length(x), 1)], amp(i, :),...
        'UniformOutput', false);
    ampmat{i} = cell2mat(ampmat{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amp = squeeze(max(datamat));
ntraces = size(datamat, 2);

t = 1 : max(cumsum(data.blockduration));
t = interp1(linspace(0, 1, length(t)), t, linspace(0, 1, ntraces))';

tracet = 1 : 15 : ntraces * 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trace = datamat(:, :, 3);


