function labels = bouts2labels(btimes, varargin)

% Converts bout intervals into a vector of labels, where each cell in
% btimes gets it's own number. for example, the indices in btimes{2} will
% be labeled 2. The function handles overlapping bouts by applying a
% precedence order specified by cellOrdr. Unassigned intervals are marked
% as the number of cells + 2.
%
% INPUT
%       btimes      1D cell array. Each cell contains n x 2 matrices
%                   representing start and stop times of bouts for each state.
%       cellOrdr   Numeric vector. Specifies the precedence order of states.
%                   Higher precedence states overwrite lower ones.
%       nlabels     scalar. length of labels vector. if empty will be taken
%                   to be max(btimes{:})
%
% OUTPUT
%       labels      Numeric vector. State labels reconstructed from bouts.
%                   Unassigned intervals are labeled as nstates + 2 (bin state).
%
% EXAMPLE
%       labels = bouts2labels('btimes', boutTimes, 'cellOrdr', [2, 1, 3], ...
%           'nstates', 6, 'labelsOrig', originalLabels);
%
% 10 Jan 25 LH      Initial version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'cellOrdr', [], @isnumeric)
addParameter(p, 'nlabels', [], @isnumeric)

parse(p, varargin{:})
cellOrdr    = p.Results.cellOrdr;
nlabels     = p.Results.nlabels;

% initialize output
if isempty(nlabels)
    nlabels = max([btimes{:}]);
end

% the use of ncell + 2 is for states, where this means unclassified
labels = ones(1, nlabels) * (length(btimes) + 2);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign state labels from bouts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for istate = 1 : length(cellOrdr)
    
    % State index for precedence
    state_idx = cellOrdr(istate);
    
    if isempty(btimes{state_idx})
        continue
    end

    % Extract indices from bout times and assign state labels
    labelsIdx = cellfun(@(x) arrayfun(@(start, stop) start : stop, ...
        x(:, 1), x(:, 2), 'UniformOutput', false), ...
        btimes(:, state_idx), 'UniformOutput', false);

    % Flatten indices and update labels
    labels(horzcat(labelsIdx{:}{:})) = state_idx;
    
end

end
