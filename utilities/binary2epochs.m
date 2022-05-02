function [epochs, nepochs] = binary2epochs(varargin)

% gets a binary vector and returns the start \ stop times of epochs. can be
% used for BS, ripples, anesthesia states etc. accounts for inter-epoch
% interval and duration. note that epoch defines start including idx and
% stop not including idx [). this is for legacy purposes. 
% 
% INPUT
%       vec         binary vector 
%       exclude     logical. exclude incomplete events (1) or not {0}.
%       minDur      minimum duration of epoch
%       maxDur      maximum duration of epoch
%       interDur    maximum inter-epoch interval
%       printFlag   logical. print to screen the number of epochs
%
% OUTPUT
%       epochs      n x 2 mat where n is the number of events and the
%                   columns represent the start and end of each event [samps]
%       nepochs     struct with fields describing number of epochs detected
%                   in each step
%     
% EXAMPLE
%       see aneStates.m or getBS.m
%       epochs = binary2epochs('vec', binaryVec, 'minDur', minDur, 'maxDur', [],...
%           'interDur', interDur, 'exclude', false);
% 
% 14 jan 19 LH      updates:
% 05 feb 22         nepochs and printFlag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'vec', [])
addParameter(p, 'exclude', false, @islogical)
addParameter(p, 'minDur', [], @isnumeric)
addParameter(p, 'maxDur', [], @isnumeric)
addParameter(p, 'interDur', [], @isnumeric)
addParameter(p, 'printFlag', true, @islogical)

parse(p, varargin{:})
vec         = p.Results.vec;
exclude     = p.Results.exclude;
minDur      = p.Results.minDur;
maxDur      = p.Results.maxDur;
interDur    = p.Results.interDur;
printFlag   = p.Results.printFlag;

% initialize output
epochs = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find start\stop times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start = find([0; diff(vec)] > 0);
stop = find([0; diff(vec)] < 0);

if isempty(start) && isempty(stop)
    if unique(vec) == 1
        epochs = [1 length(vec)];
        nepochs.detection = 1;
    else
        nepochs.detection = 0;
        if printFlag
            fprintf('\nNo epochs detected\n');
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

epochs = [start(:), stop(:)];
nepochs.detect = size(epochs, 1);
if printFlag
    fprintf('\nAfter initial detection: %d events\n', nepochs.detect);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge events if inter-event period is too short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(interDur)
    iei = epochs(2 : end, 1) - epochs(1 : end - 1, 2);     
    while min(iei) < interDur
        temp = [];
        event = epochs(1, :);
        for i = 2 : nepochs.thr
            if epochs(i, 1) - event(2) < interDur
                event = [event(1) epochs(i, 2)];
            else
                temp = [temp; event];
                event = epochs(i, :);
            end
        end
        temp = [temp; event];
        epochs = temp;
        iei = epochs(2 : end, 1) - epochs(1 : end - 1, 2);
    end
    
    % correct last epoch to end of recording
    if length(vec) - epochs(end, end) < interDur
        epochs(end, end) = length(vec);
    end    
end

nepochs.merge = size(epochs, 1);
if isempty(epochs)
    if printFlag
        fprintf('Event merge failed');
    end
    return
end
if printFlag
    fprintf('After merging: %d events\n', nepochs.merge);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard events by duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dur = epochs(:, 2) - epochs(:, 1);
if ~isempty(maxDur)
    idxMax = dur > maxDur;
else
    idxMax = zeros(nepochs.merge, 1);
end
if ~isempty(minDur)
    idxMin = dur < minDur;
else
    idxMin = zeros(nepochs.merge, 1);
end
idx = idxMax | idxMin;
epochs(idx, :) = [];

nepochs.dur = size(epochs, 1);
if isempty(epochs)
    if printFlag
        fprintf('duration test failed.');
    end
    return
end
if printFlag
    fprintf('After duration: %d events\n', nepochs.dur);
end

end

% EOF

