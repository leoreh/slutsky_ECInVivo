function hFig = mcu_frQq(tbl, varargin)
% MCU_FRQQ Performs QQ plots for firing rate distributions.
%
%   HFIG = MCU_FRQQ(TBL) generates a 2x2 tiling of QQ plots to assess
%   Log-Normality, group differences, and stability (BSL vs Recovery) of
%   firing rate distributions.
%
%   This function is useful for determining the appropriate cutoff threshold
%   during preprocessing.
%
%   INPUTS:
%       tbl         - (table) Data table containing 'fr' and grouping vars.
%                     If dataSet='vivo', must contain 'Day'.
%                     If dataSet='mea', must contain 'frSs'.
%       varargin    - (param/value) Optional parameters:
%                     'dataSet' : (char) 'vivo' (default) or 'mea'.
%                     'cutoff_Z': (scalar) Z-score for cutoff line.
%                                 Default: -1.5.
%
%   OUTPUTS:
%       hFig        - (Handle) Figure handle.
%
%   See also: MCU_RCVVEC, MEA_RCVVEC, QQPLOT

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'dataSet', 'vivo', @(x) any(validatestring(x, {'vivo', 'mea'})));
addParameter(p, 'cutoff_Z', -1.5, @isnumeric);

parse(p, tbl, varargin{:});
dataSet = p.Results.dataSet;
cutoff_Z = p.Results.cutoff_Z;


%% ========================================================================
%  EXTRACT DATA
%  ========================================================================

% Ensure Group variable exists
if ~ismember('Group', tbl.Properties.VariableNames)
    error('mcu_frQq:missingVar', 'Table must contain "Group" column.');
end

idxWt = tbl.Group == 'Control';
idxMcu = tbl.Group == 'MCU-KO';

if strcmpi(dataSet, 'vivo')
    
    % Vivo Logic (Long format with Day column)
    if ~ismember('Day', tbl.Properties.VariableNames)
        error('mcu_frQq:missingVar', 'Vivo data requires "Day" column.');
    end
    
    frBsl = tbl.fr(tbl.Day == 'BSL' & idxWt);
    frBac = tbl.fr(tbl.Day == 'BAC3' & idxWt);
    
    frBslMcu = tbl.fr(tbl.Day == 'BSL' & idxMcu);
    frBacMcu = tbl.fr(tbl.Day == 'BAC3' & idxMcu);
    
elseif strcmpi(dataSet, 'mea')
    
    % MEA Logic (Wide format with fr and frSs columns)
    % Assumes 'fr' is Baseline and 'frSs' is Steady State.
    if ~ismember('frSs', tbl.Properties.VariableNames)
        error('mcu_frQq:missingVar', 'MEA data requires "frSs" column.');
    end
    
    frBsl = tbl.fr(idxWt);
    frBac = tbl.frSs(idxWt); % Steady State acts as Recovery
    
    frBslMcu = tbl.fr(idxMcu);
    frBacMcu = tbl.frSs(idxMcu);
    
end


%% ========================================================================
%  PLOTTING
%  ========================================================================

hFig = figure('Position', [300 300 900 900], 'Color', 'w', 'Name', 'Data Distribution Check');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. Log-Normality Check (BSL Control)
% ------------------------------------
nexttile;
hQ = qqplot(log(frBsl)); 
xlabel('Theoretical Normal Quantiles');
ylabel('Sample Log(FR) Quantiles');
title(sprintf('Log-Normality Check (BSL Control)'));
grid on; axis square;   
hQ(1).LineWidth = 1;

% 2. Group Check (Control vs MCU)
% -------------------------------
nexttile;
hQ = qqplot(log(frBsl), log(frBslMcu));
xlabel('Control BSL Log(FR)');
ylabel('MCU BSL Log(FR)');
title('Group: Control vs MCU (BSL)');
grid on; axis square;
refline(1,0);
hQ(1).LineWidth = 1;

% 3. Stability Check (BSL vs Recovery Control)
% --------------------------------------------
nexttile;  hold on;
hQ = qqplot(log(frBsl), log(frBac));
xlabel('BSL Log(FR) Quantiles');
ylabel('Recovery Log(FR) Quantiles'); % Generalized label
title('Stability: BSL vs Recovery (Control)');
grid on; axis square;
refline(1,0); 
hQ(1).LineWidth = 1;

% Calculate and plot cutoff value in log space
logCutoff = prctile(log(frBsl), normcdf(cutoff_Z)*100);
yline(logCutoff, 'r--', 'Cutoff', 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'bottom');

% 4. Stability Check (BSL vs Recovery MCU)
% ----------------------------------------
nexttile;
hQ = qqplot(log(frBslMcu), log(frBacMcu));
xlabel('MCU BSL Log(FR)');
ylabel('MCU Recovery Log(FR)');
title('Stability: BSL vs Recovery (MCU)');
grid on; axis square;
refline(1,0); 
hQ(1).LineWidth = 1;

end
