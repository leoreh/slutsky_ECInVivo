% tdt_checkStates


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mname           = 'lh110';
basepath        = 'G:\Data\lh110\lh110_220801_092300';
blocks          = [2];
store           = 'Spar';
ch              = [1 : 3];

% load streams
[~, basename] = fileparts(basepath);
blockfiles = dir('block*');
blocknames = {blockfiles.name};
blocknames = natsort(blocknames);
blockpath = fullfile(basepath, blocknames{blocks});

streamStruct = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store,...
    'T1', 0, 'T2', 0, 'CHANNEL', ch);
data = streamStruct.streams.(store).data;
fs = streamStruct.streams.(store).fs;
   
% get states
statefile = [basename, '.sleep_states.mat'];
load(statefile, 'ss')
sstates = [1, 4];
for istate = 1 : length(sstates)
    stateIdx = round(ss.boutTimes{sstates(istate)} * fs);
    stateIdx(stateIdx < 0) = 1;
    stateIdx(stateIdx > length(data)) = length(data);
    stateData(:, istate) = median(data(:, InIntervals([1 : length(data)], stateIdx)), 2);
end






