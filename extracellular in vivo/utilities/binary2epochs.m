function epochs = binary2epochs(varargin)

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
%
% OUTPUT
%       epochs      n x 2 mat where n is the number of events and the
%                   columns represent the start and end of each event [samps]
%     
% EXAMPLE
%       see aneStates.m or getBS.m
% 
% 14 jan 19 LH. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'vec', [])
addParameter(p, 'exclude', 0, @islogical)
addParameter(p, 'minDur', [], @isnumeric)
addParameter(p, 'maxDur', [], @isnumeric)
addParameter(p, 'interDur', [], @isnumeric)

parse(p, varargin{:})
vec = p.Results.vec;
exclude = p.Results.exclude;
minDur = p.Results.minDur;
maxDur = p.Results.maxDur;
interDur = p.Results.interDur;

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
    else
    fprintf('\nNo epochs detected\n');
    end
    return
else
    fprintf('\nAfter initial detection: %d events\n', length(start));
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
nevents = length(epochs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge events if inter-event period is too short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(interDur)
    iei = epochs(2 : end, 1) - epochs(1 : end - 1, 2);     
    while min(iei) < interDur
        temp = [];
        event = epochs(1, :);
        for i = 2 : nevents
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
    
    nevents = size(epochs, 1);
    if isempty(epochs)
        fprintf('Event merge failed');
        return
    else
        fprintf('After merging: %d events\n', nevents);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discard events by duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dur = epochs(:, 2) - epochs(:, 1);
if ~isempty(maxDur)
    idxMax = dur > maxDur;
else
    idxMax = zeros(nevents, 1);
end
if ~isempty(minDur)
    idxMin = dur < minDur;
else
    idxMin = zeros(nevents, 1);
end
idx = idxMax | idxMin;
epochs(idx, :) = [];

nevents = size(epochs, 1);
if isempty(epochs)
    fprintf('duration test failed.');
    return
else
    fprintf('After duration: %d events\n', nevents);
end

end

% EOF

