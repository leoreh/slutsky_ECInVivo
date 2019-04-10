function [] = getEMG()

% extract EMG signal from TDT and convert it matlab file
% 
%
% INPUT:
%   basepath    path to recording folder {pwd}.
%   store       stream. typically {'Raw1'} or 'Raw2'
%   blocks      vector. blocks to convert {all}. e.g. [1 2 4 5];
%   interval
%
%
% CALLS:
%   TDT2mat
%
%
% 02 April 19 RA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'block', '');
addOptional(p, 'store', '');
addOptional(p, 'interval', [0 Inf], @isnumeric);
addOptional(p, 'chans', [0 : 4], @isnumeric);
addOptional(p, 'saveVar', true, @islogical);


parse(p,varargin{:})
basepath = p.Results.basepath;
block = p.Results.block;
Store = p.Results.store;
interval = p.Results.interval;
chans = p.Results.chans;
saveVar = p.Results.saveVar;


nchans = length(chans);

cd(basepath)
if isempty(filename)
    [~, filename] = fileparts(basepath);
    filename = [filename '.lfp'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataEMG = TDTbin2mat(blockpath, 'TYPE', {'streams'}, 'STORE', store, 'T1', chunks(j, 1), 'T2', chunks(j, 2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EMG.data = dataEMG.streams.EMG1.data;
EMG.fs = dataEMG.streams.EMG1.fs;
EMG.channels = length(chans) ;
EMG.time = interval(1) : 1/fs : ((interval(1) + (length(EMG.data)-1)) /fs);
EMG.interval = interval;



end