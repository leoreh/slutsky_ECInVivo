function tblOut = tbl_tNorm(tbl, varargin)
% TBL_TNORM Normalizes vector columns in a table based on a baseline window.
%
%   tblOut = TBL_TNORM(tbl, ...) normalizes numeric variables in 'tbl' using
%   statistics (mean, std, min, max) calculated from a specified window
%   ('winNorm').
%
%   This function serves as a wrapper for MATLAB's native NORMALIZE,
%   calculating the Center (C) and Scale (S) parameters from a baseline
%   window and applying them to the full data.
%
%   INPUTS:
%       tbl         - (table) Input table.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'varsInc'   - (cell/char) Variables to include. Default: All numeric.
%       'varsGrp'   - (cell/char) Variables to group by. If provided,
%                     baseline stats are calculated from the pooled window
%                     of the entire group and applied to all units in that group.
%                     This is useful for normalizing clusters to their
%                     group mean.
%       'winNorm'   - (1x2 double) [Start End] indices for the baseline window.
%                     Stats are calculated from this window. Default: Full vector.
%       'Method'    - (char) Normalization method {'percentage'}.
%                     Options match MATLAB's NORMALIZE where applicable, plus
%                     'percentage':
%                     'percentage' : x / mean(win)
%                     'zscore'     : (x - mean(win)) / std(win)
%                     'center'     : x - mean(win)
%                     'scale'      : x / std(win)
%                     'range'      : (x - min(win)) / (max(win) - min(win))
%
%       'flgGeom'   - (bool) Use Geometric statistics (default: false).
%                     If true, calculations are done in log-space:
%                     'percentage' : x / geomean(win)
%                     'zscore'     : (log(x) - mean(log(win))) / std(log(win))
%                     'center'     : log(x) - mean(log(win))
%                     'scale'      : log(x) / std(log(win))
%
%       'floorVal'  - (numeric) Min value for geometric calcs (avoid log(0)).
%                     If empty and flgGeom=true, defaults to half the minimum
%                     non-zero value of the entire variable.
%
%   EXAMPLE:
%       % Normalize 'LFP' column to baseline (indices 1-100) as percentage:
%       tbl = tbl_tNorm(tbl, 'varsInc', 'LFP', 'winNorm', [1, 100], 'Method', 'percentage');
%
%       % Z-score 'FR' by 'Genotype' group using baseline stats:
%       % Each unit is normalized by the MEAN and STD of the whole group's baseline.
%       tbl = tbl_tNorm(tbl, 'varsInc', 'FR', 'varsGrp', 'Genotype', ...
%                       'winNorm', [1, 500], 'Method', 'zscore');
%
%   See also: NORMALIZE

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'tbl', @istable);
addParameter(p, 'varsInc', [], @(x) iscell(x) || ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'varsGrp', [], @(x) iscell(x) || ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'winNorm', [], @(x) isnumeric(x) && (isempty(x) || numel(x)==2));
addParameter(p, 'Method', 'percentage', @(x) ischar(x) || isstring(x));
addParameter(p, 'flgGeom', false, @islogical);
addParameter(p, 'floorVal', [], @(x) isnumeric(x) && (isempty(x) || isscalar(x)));

parse(p, tbl, varargin{:});

varsInc  = p.Results.varsInc;
varsGrp  = p.Results.varsGrp;
winNorm  = p.Results.winNorm;
method   = p.Results.Method;
flgGeom  = p.Results.flgGeom;
floorVal = p.Results.floorVal;

% Pass unmatched arguments to normalize
unmatched = p.Unmatched;
normArgs = reshape([fieldnames(unmatched), struct2cell(unmatched)]', 1, []);

% Standardize cell arrays
if ischar(varsInc) || isstring(varsInc), varsInc = cellstr(varsInc); end
if ischar(varsGrp) || isstring(varsGrp), varsGrp = cellstr(varsGrp); end


%% ========================================================================
%  INITIALIZATION
%  ========================================================================

tblOut = tbl;
tblVars = tbl.Properties.VariableNames;

% Determine Variables
if isempty(varsInc)
    processVars = tblVars;
else
    % Filter vars that exist
    processVars = varsInc(ismember(varsInc, tblVars));
end

% Filter Numeric Variables Only
% Check first row to see if it is numeric
isNum = cellfun(@(x) isnumeric(tblOut.(x)), processVars);
processVars = processVars(isNum);

% Window Indices
wS = 1; wE = inf;
if ~isempty(winNorm)
    wS = max(1, winNorm(1));
    wE = winNorm(2);
end


%% ========================================================================
%  PRE-CALCULATE GROUPS
%  ========================================================================

if ~isempty(varsGrp)
    [uGrps, ~, ic] = unique(tblOut(:, varsGrp), 'rows');
    nGrps = height(uGrps);
    % ic contains the group index for each row, avoiding manual loop matching
    idxGrps = cell(nGrps, 1);
    for i = 1:nGrps
        idxGrps{i} = (ic == i);
    end
else
    idxGrps = {};
end


%% ========================================================================
%  TRANSFORM
%  ========================================================================

for iVar = 1:length(processVars)
    varName = processVars{iVar};
    varData = tblOut.(varName);

    % Determine Floor Value for Geometric Stats
    currFloor = 0;
    if flgGeom
        if isempty(floorVal)
            % Auto-calculate: Half of minimum non-zero value
            minVal = min(varData(varData > 0), [], 'all');
            if isempty(minVal)
                currFloor = eps;
            else
                currFloor = minVal / 2;
            end
        else
            currFloor = floorVal;
        end
    end

    % ---------------------------------------
    % Case 1: Row-wise Normalization (Default)
    % ---------------------------------------
    if isempty(idxGrps)

        % Extract Window
        winMat = varData(:, wS:wE);
        
        % Calculate Row-wise Params
        [C, S] = norm_params(winMat, method, flgGeom, currFloor);

        % Apply Normalize
        if flgGeom & ~strcmpi(method, 'percentage')
            % Log Output: (log(x) - muLog) / sigLog
            % Apply floor before log
            varData = log(max(varData, currFloor));
        end
        varData = normalize(varData, 2, 'center', C, 'scale', S, normArgs{:});

    % ---------------------------------------
    % Case 2: Group-wise Normalization
    % ---------------------------------------
    else
        for iGrp = 1:nGrps
            idx = idxGrps{iGrp};
            grpData = varData(idx, :);

            % Extract Pooled Window for Group
            subMat = grpData(:, wS:wE);
            
            % Calculate Params
            % Here we treat all samples in the window as part of the
            % baseline distribution. C and S will be scalars.
            % Note: We pass vector so norm_params sees it as one distribution
            [C, S] = norm_params([subMat(:)]', method, flgGeom, currFloor);

            % Apply Normalize to Group Rows
            if flgGeom & ~strcmpi(method, 'percentage')
                 % Log Output
                 grpData = log(max(grpData, currFloor));
            end
            varData(idx, :) = normalize(grpData, 2, 'center', C, 'scale', S, normArgs{:});
        end
    end

    tblOut.(varName) = varData;
end

end


%% ========================================================================
%  HELPER FUNCTION
%  ========================================================================

function [C, S] = norm_params(vals, method, flgGeom, floorVal)
% Calculates Center (C) and Scale (S) parameters.
% vals: Data to calculate stats on.

if nargin < 3, flgGeom = false; end
if nargin < 4, floorVal = eps; end

% Geometric Logic
if flgGeom
    vals = max(vals, floorVal); % Cap minimum
    
    % Work in log space
    logVals = log(vals);
    
    muLog = mean(logVals, 2, 'omitnan');
    sigLog = std(logVals, 0, 2, 'omitnan');
    
    switch lower(method)
        case 'percentage'
            % For percentage: x / GM => C=0, S=GM
            C = 0;
            S = exp(muLog); % Geometric Mean
            
        case 'zscore'
            % (log(x) - muLog) / sigLog
            C = muLog;
            S = sigLog;
            
        case 'center'
            % log(x) - muLog
            C = muLog;
            S = 1;
            
        case 'scale'
            % log(x) / sigLog
            C = 0;
            S = sigLog;
            
         case 'range'
             % log(x) mapped? Not typical, fallback to zscore log
             C = muLog;
             S = sigLog;
             
        otherwise
             C = muLog;
             S = sigLog;
    end
    
    % Safety
    S(S == 0) = 1; % If S is 0 (const), avoid inf.
    
else
    % Arithmetic Logic
    switch lower(method)
        case 'percentage'
            mu = mean(vals, 2, 'omitnan');
            C = 0;
            S = mu;

        case 'zscore'
            mu = mean(vals, 2, 'omitnan');
            sig = std(vals, 0, 2, 'omitnan');
            C = mu;
            S = sig;

        case 'center'
            mu = mean(vals, 2, 'omitnan');
            C = mu;
            S = 1;

        case 'scale'
            sig = std(vals, 0, 2, 'omitnan');
            C = 0;
            S = sig;

        case 'range'
            mn = min(vals, [], 2, 'omitnan');
            mx = max(vals, [], 2, 'omitnan');
            C = mn;
            S = mx - mn;

        otherwise
            % Default to zscore
            mu = mean(vals, 2, 'omitnan');
            sig = std(vals, 0, 2, 'omitnan');
            C = mu;
            S = sig;
    end
    
    % Safety for division by zero
    S(S == 0) = eps;
end

end
