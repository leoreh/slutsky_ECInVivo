function [bouts, nbouts] = binary2bouts(varargin)

% gets a binary vector and returns the start \ stop times of bouts. can be
% used for BS, ripples, anesthesia states etc. accounts for inter-bout
% interval and duration. note that bout defines start including idx and
% stop not including idx [). this is for legacy purposes. 
% 
% INPUT
%       vec         binary vector 
%       exclude     logical. exclude incomplete events (1) or not {0}.
%       minDur      minimum duration of bout
%       maxDur      maximum duration of bout
%       interDur    maximum inter-bout interval
%       flgPrnt     logical. print to screen the number of bouts
%
% OUTPUT
%       bouts      n x 2 mat where n is the number of events and the
%                   columns represent the start and end of each event [samps]
%       nbouts     struct with fields describing number of bouts detected
%                   in each step
%     
% EXAMPLE
%       see aneStates.m or getBS.m
%       bouts = binary2bouts('vec', binaryVec, 'minDur', minDur, 'maxDur', [],...
%           'interDur', interDur, 'exclude', false);
% 
% 14 jan 19 LH      updates:
% 05 feb 22         nbouts and flgPrnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'vec', [])
addParameter(p, 'exclude', false, @islogical)
addParameter(p, 'minDur', [], @isnumeric)
addParameter(p, 'maxDur', [], @isnumeric)
addParameter(p, 'interDur', [], @isnumeric)
addParameter(p, 'flgPrnt', true, @islogical)

parse(p, varargin{:})
vec         = p.Results.vec;
exclude     = p.Results.exclude;
minDur      = p.Results.minDur;
maxDur      = p.Results.maxDur;
interDur    = p.Results.interDur;
flgPrnt     = p.Results.flgPrnt;

% initialize output
bouts = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find start\stop times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vec = vec(:);
start = find([0; diff(vec)] > 0);
stop = find([0; diff(vec)] < 0);

if isempty(start) && isempty(stop)
    if unique(vec) == 1
        bouts = [1 length(vec)];
        nbouts.detect = 1;
        nbouts.dur = 1;
        nbouts.merge = 1;
    else
        nbouts.detect = 0;
        nbouts.dur = 0;
        nbouts.merge = 0;
        if flgPrnt
            fprintf('\nNo bouts detected\n');
        end
    end
    return
end

% exclude last event if it is incomplete
if length(stop) == length(start) - 1
    if exclude
        start = start(1 : end - 1);
    else
        stop = [stop; length(vec)];
    end
end

% complete first event if it is incomplete
if length(stop) - 1 == length(start)
    if exclude
        stop = stop(2 : end);
    else
        start = [1; start];
    end
end

% correct special case when both first and last events are incomplete
if start(1) > stop(1)
    if exclude
        stop(1) = [];
        start(end) = [];
    else
        start = [1; start];
        stop = [stop; length(vec)];
    end
end

bouts = [start(:), stop(:)];
nbouts.detect = size(bouts, 1);
if flgPrnt
    fprintf('\nAfter initial detection: %d events\n', nbouts.detect);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge events if inter-event period is too short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(interDur)
    iei = bouts(2 : end, 1) - bouts(1 : end - 1, 2);     
    while min(iei) < interDur
        temp = [];
        event = bouts(1, :);
        for i = 2 : nbouts.detect
            if bouts(i, 1) - event(2) < interDur
                event = [event(1) bouts(i, 2)];
            else
                temp = [temp; event];
                event = bouts(i, :);
            end
        end
        temp = [temp; event];
        bouts = temp;
        iei = bouts(2 : end, 1) - bouts(1 : end - 1, 2);
    end
    
    % correct last bout to end of recording
    if length(vec) - bouts(end, end) < interDur
        bouts(end, end) = length(vec);
    end    
end

nbouts.merge = size(bouts, 1);
if isempty(bouts)
    if flgPrnt
        fprintf('Event merge failed');
    end
    return
end
if flgPrnt
    fprintf('After merging: %d events\n', nbouts.merge);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard events by duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dur = bouts(:, 2) - bouts(:, 1);
if ~isempty(maxDur)
    idxMax = dur > maxDur;
else
    idxMax = zeros(nbouts.merge, 1);
end
if ~isempty(minDur)
    idxMin = dur < minDur;
else
    idxMin = zeros(nbouts.merge, 1);
end
idx = idxMax | idxMin;
bouts(idx, :) = [];

nbouts.dur = size(bouts, 1);
if isempty(bouts)
    if flgPrnt
        fprintf('duration test failed.');
    end
    return
end
if flgPrnt
    fprintf('After duration: %d events\n', nbouts.dur);
end

end

% EOF

