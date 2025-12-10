function [frMat, alignIdx] = alignFR2pnt(varargin)

% reveives a cell of basepaths, loads the fr / sr struct, and combines
% firing rate vs. time for all units / tetrodes. also gets a vec of time
% points [s] (one for each session) which will be converted to indices to
% the fr mat. the mats will be combined such that these points
% are aligned. EXAMPLE x1 = [1 0 1 0]; x2 = [2 2 3 2 4]; pnts = [2 3];
% output: mat = [nan 1 0 1 0; 2 2 3 2 4];

% INPUT
%   basepaths       cell of chars describing the full path to sessions
%   timeIdx         numeric. points by which to align fr mats. if empty
%                   will be extracted from session struct
%   iunit           numeric. extract RS (1), FS (2), or all (empty). only
%                   relavent for fr mat
%   suFlag          logical. if true will concat fr, if false sr
%   dataType        char. use 'strd' or 'norm'.

% UPDATES
%   13 apr 23       add iunit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepaths', {});
addOptional(p, 'timeIdx', [], @isnumeric);
addOptional(p, 'iunit', [], @isnumeric);
addOptional(p, 'suFlag', true, @islogical);
addOptional(p, 'dataType', 'strd', @ischar);

parse(p, varargin{:})
basepaths   = p.Results.basepaths;
timeIdx     = p.Results.timeIdx;
iunit       = p.Results.iunit;
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
if suFlag
    varsFile = ["fr"; "datInfo"; "session"; "units"];
else
    varsFile = ["sr"; "datInfo"; "session"; "units"];
end
v = basepaths2vars('basepaths', basepaths, 'vars', varsFile);

if ~suFlag && isfield(v, 'fr') && ~isfield(v, 'sr')
    [v.sr] = v.fr;
    v = rmfield(v, 'fr');
end
nsessions = length(basepaths);

% go over each session and find the index to fr.strd mat
for isession = 1 : nsessions

    % get time indices
    if isempty(timeIdx)
        timepnt = v(isession).session.general.timepnt;
    else
        timepnt = timeIdx(isession);
    end
    [~, alignIdx(isession)] = min(abs(v(isession).(dataStruct).tstamps - timepnt));
    tLen(isession) = length(v(isession).(dataStruct).tstamps);

    % get unit indices
    if ~isempty(iunit)
        units{isession} = v(isession).units.clean(iunit, :);
        nunits(isession) = sum(units{isession});
    else
        nunits(isession) = size(v(isession).(dataStruct).(dataType), 1);
        units{isession} = ones(1, nunits(isession));
    end
end

% get data mat
frMat = nan(sum(nunits), max(alignIdx) + max(tLen));        % initialize
cnt = 1;
for isession = 1 : nsessions
    startIdx = max(alignIdx) - alignIdx(isession) + 1;
    frIdx = startIdx : tLen(isession) + startIdx - 1;
    frMat(cnt : cnt + nunits(isession) - 1, frIdx) =...
        v(isession).(dataStruct).(dataType)(units{isession}, :);
    cnt = cnt + nunits(isession);
end
frMat = frMat';

end

% EOF
