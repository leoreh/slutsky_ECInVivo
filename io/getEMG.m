function emg = getEMG(varargin)

% loads EMG from tank (TDT) and finds linear envelope.
%
% INPUT:
%   basepath    path to recording folder {pwd}.
%   store       stream. typically {'EMG1'} 
%   blocks      vector. blocks to convert {all}. e.g. [1 2 4 5];
%   ch          channels to extract
%   graphics    plot figure {1}.
%   saveVar     save variable {1}.
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
%   adapt chunks from tdt2dat
%
% 29 apr 19 LH
% 15 sep 19 LH channel handeling and arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'store', 'EMG1', @ischar);
addOptional(p, 'blocks', [], @isnumeric);
addOptional(p, 'ch', 1, @isnumeric);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'saveVar', true, @islogical);

parse(p, varargin{:})
basepath = p.Results.basepath;
store = p.Results.store;
blocks = p.Results.blocks;
ch = p.Results.ch;
graphics = p.Results.graphics;
saveVar = p.Results.saveVar;

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
        'CHANNEL', ch, 'T1', 0, 'T2', 0);
    raw = raw.streams.(store).data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find linear envelope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% OPTION 1:
% [EMG, ~] = envelope(double(EMG2), win, 'rms');

%%% OPTION 2:
% bandpass filter
data = filterLFP(double(raw'), 'fs', fs,...
    'type', 'butter', 'passband', [10 500], 'graphics', false);

% rectify
data = abs(data - mean(data));

% % low-pass filter
win = round(0.5 * fs);    % 500 ms moving average
data = movmean(data, win);

% normalize
data = data ./ max(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    nchans = size(data, 2);
%     idx = 1 : min(size(data, 1), fs * 60);
    idx = 1 : size(data, 1);
    t = idx / fs / 60;
    
    for i = 1 : nchans
        figure
        yyaxis left
        plot(t, raw(i, idx))
        set(gca,'TickLength',[0 0])
        box off
        axis tight
        
        yyaxis right
        plot(t, data(idx, i))
        yticks([0 1])
        set(gca,'TickLength',[0 0])
        box off
        axis tight
        ylim([0 1])
    end
    xlabel('time [m]')
    legend('Raw', 'Envelope')
    title('EMG')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save var
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar
    emg.data = data;
    emg.raw = raw';
    emg.filename = basename;
    emg.blocks = blocks;
    emg.blockduration = nsec;
    emg.ch = ch;
    emg.fs = heads.stores.(store).fs;
    
    save([basename, '.emg.mat'], 'emg');
end

end

% EOF