function tblSynth = mcu_matchQntl(tbl, nBins, varargin)
% MCU_MATCHQNTL  Perform quantile matching to create synthetic units.
%
%   tblSynth = MATCH_QNTL(tbl, nBins, 'flgPool', true) matches units
%   based on firing rate ranks.
%
%   INPUTS:
%       tbl     - (table) Source data with vars: Group, Day, Name, and numeric data.
%       nBins   - (int) Number of quantiles (e.g., 5, 10).
%       varargin:
%           'flgPool' - (bool) If true, pools all units by Group. 
%                       If false (default), bins per Animal.
%           'var'     - (char) Variable to use for sorting/ranking (default: 'fr').
%           'avgType' - (char) 'mean', 'median' (default), or 'geomean'.
%
%   OUTPUTS:
%       tblSynth - (table) Table of synthetic units with columns:
%                  [Group, Name, Bin, Var (BSL), ss_Var (BAC3), ...]
%
%   EXAMPLE:
%       tblSynth = match_qntl(tbl, 5, 'flgPool', false, 'var', 'fr', 'avgType', 'median');
%

%% ========================================================================
%  INPUTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'nBins', @isnumeric);
addParameter(p, 'flgPool', false, @islogical);
addParameter(p, 'var', 'fr', @ischar);
addParameter(p, 'avgType', 'median', @(x) any(validatestring(x, {'mean', 'median', 'geomean'})));
parse(p, tbl, nBins, varargin{:});

flgPool = p.Results.flgPool;
sortVar = p.Results.var; % Rename to avoid conflict
avgType = p.Results.avgType;

% Identify all numeric variables to process
% Exclude metadata columns
metaVars = {'Group', 'Name', 'Day', 'Bin', 'id', 'Timepoint'}; 
allVars = tbl.Properties.VariableNames;

% Find numeric vars
isNum = varfun(@isnumeric, tbl, 'OutputFormat', 'uniform');
numVars = allVars(isNum);

% Remove metadata from numeric list if they were flagged numeric
varsToProc = setdiff(numVars, metaVars);

% Ensure sort variable exists
if ~ismember(sortVar, tbl.Properties.VariableNames)
     error('Sort variable "%s" not found in table.', sortVar);
end


%% ========================================================================
%  QUANTILE MATCHING
%  ========================================================================

tblSynth = table();

if flgPool
    % ---------------------------------------------------------------------
    % OPTION A: POOLED MATCHING (By Group)
    % ---------------------------------------------------------------------
    % Ignores animal identity, treats all units in a group as one population.
    
    grps = unique(tbl.Group);
    
    for iGrp = 1:length(grps)
        grp = char(grps(iGrp));
        
        % Filter by Group
        tGrp = tbl(tbl.Group == grp, :);
        
        % Process BSL/BAC3
        tRow = proc_pair(tGrp, nBins, varsToProc, sortVar, avgType);
        
        if ~isempty(tRow)
            % Add Metadata
            tRow.Group = repmat(categorical({grp}), height(tRow), 1);
            tRow.Name = repmat({'Pooled'}, height(tRow), 1);
            
            % Reorder
            tRow = movevars(tRow, {'Group', 'Name', 'Bin'}, 'Before', 1);
            
            % Append
            tblSynth = [tblSynth; tRow];
        end
    end
    
else
    % ---------------------------------------------------------------------
    % OPTION B: ANIMAL MATCHING (Per Animal)
    % ---------------------------------------------------------------------
    % Preserves animal identity, matches units within each animal.
    
    animals = unique(tbl.Name);
    
    for iAni = 1:length(animals)
        
        name = animals(iAni); 
        
        % Filter by Animal
        tAni = tbl(ismember(tbl.Name, name), :);
        
        % Process BSL/BAC3
        tRow = proc_pair(tAni, nBins, varsToProc, sortVar, avgType);
        
        if ~isempty(tRow)
            % Get Group (First value, assume constant)
            grp = char(tAni.Group(1));
            
            % Add Metadata
            tRow.Group = repmat(categorical({grp}), height(tRow), 1);
            tRow.Name = repmat(name, height(tRow), 1);
            
            % Reorder
            tRow = movevars(tRow, {'Group', 'Name', 'Bin'}, 'Before', 1);
            
            % Append
            tblSynth = [tblSynth; tRow];
        end
    end
end

end


%% ========================================================================
%  HELPER: PROCESS PAIR (BSL & BAC3)
%  ========================================================================
function tWide = proc_pair(tIn, nBins, vars, sortVar, avgType)

    tWide = table();

    % Split
    tBsl = tIn(tIn.Day == 'BSL', :);
    tBac = tIn(tIn.Day == 'BAC3', :);
    
    % Sort by Specified Variable (Descending)
    tBsl = sortrows(tBsl, sortVar, 'descend');
    tBac = sortrows(tBac, sortVar, 'descend');

    % Bin
    tBslBin = bin_table(tBsl, nBins, vars, avgType);
    tBacBin = bin_table(tBac, nBins, vars, avgType);
    
    % Rename Columns (Wide Format)
    % Baseline: No Suffix (Canonical)
    % SteadyState: Prefix 'ss_' (MEA Style)
    
    % tBslBin variables remain as is (e.g. 'fr', 'frBspk')
    
    % tBacBin variables get 'ss_' prefix
    tBacBin.Properties.VariableNames = cellfun(@(x) ['ss_' x], ...
        tBacBin.Properties.VariableNames, 'uni', false);
    
    % Combine Horizontal
    % Note: tBslBin has 'Bin', tBacBin has 'ss_Bin'
    % We assume they align by row index 1..nBins
    tWide = [tBslBin, tBacBin];
    
    % Clean up helper bin column from steady state
    if ismember('ss_Bin', tWide.Properties.VariableNames)
        tWide.ss_Bin = [];
    end
    
end


%% ========================================================================
%  HELPER: BIN TABLE
%  ========================================================================
function tBin = bin_table(tIn, nBins, vars, avgType)

    nUnits = height(tIn);
    
    % Create Bins (1 = Top Value, N = Bottom Value)
    % Using ceil to distribute evenly
    binIdx = ceil((1:nUnits)' / nUnits * nBins);
    
    % Pre-allocate
    tBin = table();
    tBin.Bin = (1:nBins)';
    
    % Select Aggregation Function
    switch avgType
        case 'mean'
            fh = @mean;
        case 'median'
            fh = @median;
        case 'geomean'
            fh = @geomean;
        otherwise
            error('Invalid avgType. Use ''mean'', ''median'', or ''geomean''.');
    end
    
    for iV = 1:length(vars)
        v = vars{iV};
        
        % Robust Aggregation (NaN if empty bin)
        vals = accumarray(binIdx, tIn.(v), [nBins 1], fh, NaN);
        
        tBin.(v) = vals;
    end
    
end
