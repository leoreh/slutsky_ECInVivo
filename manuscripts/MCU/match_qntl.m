function tblSynth = match_qntl(tbl, nBins, varargin)
% MATCH_QNTL  Perform quantile matching to create synthetic units.
%
%   tblSynth = MATCH_QNTL(tbl, nBins, 'flgPool', true) matches units
%   based on firing rate ranks.
%
%   INPUTS:
%       tbl     - (table) Source data with vars: fr, Group, Day, Name.
%       nBins   - (int) Number of quantiles (e.g., 5, 10).
%       varargin:
%           'flgPool' - (bool) If true, pools all units by Group. 
%                       If false (default), bins per Animal.
%           'vars'    - (cell) Variables to average e.g., {'fr', 'pBspk'}.
%
%   OUTPUTS:
%       tblSynth - (table) Table of synthetic units with columns:
%                  [Group, Name, Bin, Var_BSL, Var_BAC3, ...]
%
%   EXAMPLE:
%       tblSynth = match_qntl(tbl, 5, 'flgPool', false);
%

%% ========================================================================
%  INPUTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'nBins', @isnumeric);
addParameter(p, 'flgPool', false, @islogical);
addParameter(p, 'vars', {'fr', 'frBspk', 'frSspk', 'pBspk'}, @iscell);
parse(p, tbl, nBins, varargin{:});

flgPool = p.Results.flgPool;
vars = p.Results.vars;


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
        tRow = proc_pair(tGrp, nBins, vars);
        
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
        tRow = proc_pair(tAni, nBins, vars);
        
        if ~isempty(tRow)
            % Get Group (First value, assume constant)
            grp = char(tAni.Group(1));
            
            % Add Metadata
            tRow.Group = repmat(categorical({grp}), height(tRow), 1);
            tRow.Name = repmat({name}, height(tRow), 1);
            
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
function tWide = proc_pair(tIn, nBins, vars)

    tWide = table();

    % Split
    tBsl = tIn(tIn.Day == 'BSL', :);
    tBac = tIn(tIn.Day == 'BAC3', :);
    
    % Sort by Firing Rate (Descending)
    tBsl = sortrows(tBsl, 'fr', 'descend');
    tBac = sortrows(tBac, 'fr', 'descend');

    % Bin
    tBslBin = bin_table(tBsl, nBins, vars);
    tBacBin = bin_table(tBac, nBins, vars);
    
    % Rename Columns (Wide Format)
    tBslBin.Properties.VariableNames = cellfun(@(x) [x '_BSL'], ...
        tBslBin.Properties.VariableNames, 'uni', false);
    tBacBin.Properties.VariableNames = cellfun(@(x) [x '_BAC3'], ...
        tBacBin.Properties.VariableNames, 'uni', false);
    
    % Combine Horizontal
    % Note: tBslBin and tBacBin both have 'Bin_BSL'/'Bin_BAC3' derived indices
    % We assume they align by row index 1..nBins
    tWide = [tBslBin, tBacBin];
    
    % Fix Bin Column
    tWide.Bin_BAC3 = [];
    tWide = renamevars(tWide, 'Bin_BSL', 'Bin');
    
end


%% ========================================================================
%  HELPER: BIN TABLE
%  ========================================================================
function tBin = bin_table(tIn, nBins, vars)

    nUnits = height(tIn);
    
    % Create Bins (1 = Top FR, N = Bottom FR)
    % Using ceil to distribute evenly
    binIdx = ceil((1:nUnits)' / nUnits * nBins);
    
    % Pre-allocate
    tBin = table();
    tBin.Bin = (1:nBins)';
    
    for iV = 1:length(vars)
        v = vars{iV};
        
        % Robust Mean (NaN if empty bin)
        vals = accumarray(binIdx, tIn.(v), [nBins 1], @mean, NaN);
        
        tBin.(v) = vals;
    end
    
end
