function [stateEpochs, epochStats] = as_epochs(varargin)

% recieves a labels vector and returns the state epochs. basically a
% wrapper for binary2epochs with params adjusted for sleep scoring
%
% INPUT:
%   labels          numeric. see as_classify 
%   minDur          numeric. minimum duration of an epoch. 
%                   if length(minDur) == nstates than a different minimum
%                   duration will be applied to each state
%   interDur        numeric. combine epochs separated by <= interDur
%   confMarg        2 x 1 numeric. confidance margin. for example, [1 2]
%                   will change the epoch [789 810] to [790 808]. 
%                   this is to assure that the previous and next state to
%                   not influence the current state
% 
% OUTPUT
%
% DEPENDENCIES
% 
% TO DO LIST
%
% 12 jan 22 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'labels', [], @isnumeric);
addOptional(p, 'minDur', [10, 5, 5, 10, 5, 5], @isnumeric);
addOptional(p, 'interDur', 4, @isnumeric);
addOptional(p, 'confMarg', [], @isnumeric);

parse(p, varargin{:})
labels          = p.Results.labels;
minDur          = p.Results.minDur;
interDur        = p.Results.interDur;
confMarg        = p.Results.confMarg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% state params
cfg = as_loadConfig();
nstates = cfg.nstates;

if length(minDur) == 1
    minDur = repmat(minDur, nstates, 1);
elseif length(minDur) ~= nstates
    error('minDur length is different than the number of states')
end
if length(interDur) == 1
    interDur = repmat(interDur, nstates, 1);
elseif length(interDur) ~= nstates
    error('interDur length is different than the number of states')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create state epochs
for istate = 1 : nstates
    binaryVec = zeros(length(labels), 1);
    binaryVec(labels == istate) = 1;
    stateEpochs{istate} = binary2epochs('vec', binaryVec, 'minDur', minDur(istate), 'maxDur', [],...
        'interDur', interDur(istate), 'exclude', false);
end

% remove confindance margins
if ~isempty(confMarg)
    funh = @(x) [x(:, 1) + confMarg(1) x(:, 2) - confMarg(2)];
    stateEpochs = cellfun(funh, stateEpochs, 'uni', false);
end

% epoch stats
epochStats.epLen = cellfun(@(x) (diff(x')'), stateEpochs, 'UniformOutput', false);
epochStats.nepochs = cellfun(@length, epochStats.epLen);
epochStats.totDur = cellfun(@sum, epochStats.epLen);
epochStats.confMarg = confMarg;

end

% EOF