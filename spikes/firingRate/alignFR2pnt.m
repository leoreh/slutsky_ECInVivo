function [frMat, timeIdx] = alignFR2pnt(varargin)

% reveives a cell of basepaths, loads the fr / sr struct, and combines
% firing rate vs. time for all units / tetrodes. also gets a vec of time
% points [s] (one for each session) which will be converted to indices to
% the fr mat. the mats will be combined such that these points
% are aligned. EXAMPLE x1 = [1 0 1 0]; x2 = [2 2 3 2 4]; pnts = [2 3];
% output: mat = [nan 1 0 1 0; 2 2 3 2 4];

% INPUT
%   basepaths       cell of chars describing the full path to sessions
%   suFlag          logical. if true will concat fr, if false sr
%   dataType        char. use 'strd' or 'norm'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepaths', {});
addOptional(p, 'suFlag', true, @islogical);
addOptional(p, 'dataType', 'strd', @ischar);

parse(p, varargin{:})
basepaths   = p.Results.basepaths;
suFlag      = p.Results.suFlag;
dataType    = p.Results.dataType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if suFlag
    dataStruct = 'fr';
else
    dataStruct = 'sr';
end

% load vars from each session
varsFile = ["fr"; "sr"; "datInfo"; "session"];
varsName = ["fr"; "sr"; "datInfo"; "session"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

% go over each session and find the index to the fr.strd mat based on the
% time point saved in the session struct
for isession = 1 : nsessions
    timepnt = v(isession).session.general.timepnt;
    [~, timeIdx(isession)] = min(abs(v(isession).(dataStruct).tstamps - timepnt));
    tLen(isession) = length(v(isession).(dataStruct).tstamps);
    nunits(isession) = size(v(isession).(dataStruct).(dataType), 1);
end

% initialize
frMat = nan(sum(nunits), max(timeIdx) + max(tLen));
cnt = 1;
for isession = 1 : nsessions
    startIdx = max(timeIdx) - timeIdx(isession) + 1;
    frIdx = startIdx : tLen(isession) + startIdx - 1;
    frMat(cnt : cnt + nunits(isession) - 1, frIdx) = v(isession).(dataStruct).(dataType);
    cnt = cnt + nunits(isession);
end

end

% EOF
