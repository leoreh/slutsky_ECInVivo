function emg = getEMG(basepath, store, blocks, rmvch, graphics)

% loads EMG from tank (TDT) and finds linear envelope.
%
% INPUT:
%   basepath    path to recording folder {pwd}.
%   store       stream. typically {'Raw1'} or 'Raw2'
%   blocks      vector. blocks to convert {all}. e.g. [1 2 4 5];
%   rmvch       channels to remove (according to original order) {[]}
%
% OUTPUT
%   emg         struct with fields:
%       raw         raw EMG
%       data        linear envelope
%       <params>    as in input + tdt params
%
% CALLS:
%   TDTbin2mat
%
% TO DO LIST:
%   handle arguments
%
% 29 apr 19 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(basepath)
    basepath = pwd;
end
if nargin < 2 || isempty(store)
    store = 'EMG1';
end
if nargin < 3 || isempty(blocks)
    blocks = [];
end
if nargin < 4 || isempty(rmvch)
    rmvch = [];
end
if nargin < 5 || isempty(graphics)
    graphics = true;
end

[~, basename] = fileparts(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get tank blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(basepath)
blockfiles = dir('block*');
blocknames = {blockfiles.name};
fprintf(1, '\nFound %d blocks in %s\n\n', length(blocknames), basepath);

if isempty(blocknames)
    error('no blocks in dir %s.', basepath)
end
if ~isempty(blocks)
    blocknames = blocknames(blocks);
end
nblocks = length(blocknames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load EMG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : nblocks
    blockpath = fullfile(basepath, blocknames{i});
    fprintf(1, 'Working on %s\n', blocknames{i});
    
    heads = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store, 'HEADERS', 1);
    nsec(i) = heads.stores.(store).ts(end);
    fs = heads.stores.(store).fs;
    
    raw = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store,...
        'T1', 0, 'T2', 0);
    raw = raw.streams.(store).data;
end

if ~isempty(rmvch)
    raw(rmvch, :) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find linear envelope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% OPTION 1:
% [EMG, ~] = envelope(double(EMG2), win, 'rms');

%%% OPTION 2:
% bandpass filter
% data = filterLFP(double(raw'), 'fs', fs,...
%     'type', 'butter', 'passband', [10 500], 'graphics', false);
% 
% % rectify
% data = abs(data - mean(data));
% 
% % low-pass filter
% win = round(0.5 * fs);    % 500 ms moving average
% data = movmean(data, win);
% 
% % normalize
% data = data ./ max(data);
data = double(raw');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    nchans = size(data, 2);
    t = (1 : size(data, 1)) / fs / 60;
    
    figure
    for i = 1 : nchans
        subplot(nchans, 1, i)
        plot(t, data(:, i))
        yticks([0 1])
        set(gca,'TickLength',[0 0])
        box off
        axis tight
        ylim([0 1])
        if i == 1
            title('EMG')
        end
    end
    xlabel('time [m]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save var
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emg.data = data;
emg.raw = raw';
emg.filename = basename;
emg.blocks = blocks;
emg.blockduration = nsec;
emg.rmvch = rmvch;
emg.fs = heads.stores.(store).fs;

save([basename, '.emg.mat'], 'emg');

end

% EOF