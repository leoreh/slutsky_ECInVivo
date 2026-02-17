function [hFig, qqData] = mcu_rcvQq(tbl, varargin)
% MCU_RCVQQ Performs QQ plots for distribution analysis.
%
%   HFIG = MCU_RCVQQ(TBL, 'var', 'fr') generates a 2x2 tiling of QQ plots 
%   for the specified variable to assess Log-Normality, group differences, 
%   and stability (BSL vs Recovery/SteadyState).
%
%   This function assists in determining cutoff thresholds and validating
%   distributional assumptions.
%
%   INPUTS:
%       tbl         - (table) Data table.
%                     If dataSet='vivo', must contain 'Day'.
%                     If dataSet='mea', must contain 'ss_<var>' 
%       varargin    - (param/value) Optional parameters:
%                     'var'     : (char) Variable to plot (default: 'fr').
%                     'dataSet' : (char) 'vivo' (default) or 'mea'.
%                     'cutoff_z': (scalar) Z-score for cutoff line (default: -1.5).
%
%   OUTPUTS:
%       hFig        - (Handle) Figure handle.
%       qqData      - (table) Processed plot data. Col 'Data' is [x, y].
%
%   See also: MCU_RCVVEC, MEA_RCVVEC, QQPLOT
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'var', 'fr', @ischar);
addParameter(p, 'dataSet', 'vivo', @(x) any(validatestring(x, {'vivo', 'mea'})));
addParameter(p, 'cutoff_z', -1.5, @isnumeric);

parse(p, tbl, varargin{:});
varName = p.Results.var;
dataSet = p.Results.dataSet;
cutoff_z = p.Results.cutoff_z;


%% ========================================================================
%  EXTRACT DATA
%  ========================================================================

% Ensure Group variable exists
if ~ismember('Group', tbl.Properties.VariableNames)
    error('mcu_rcvQq:missingVar', 'Table must contain "Group" column.');
end

idxWt = tbl.Group == 'Control';
idxMcu = tbl.Group == 'MCU-KO';

if strcmpi(dataSet, 'vivo')
    
    % Vivo Logic (Long format with Day column)
    if ~ismember('Day', tbl.Properties.VariableNames)
        error('mcu_rcvQq:missingVar', 'Vivo data requires "Day" column.');
    end
    
    % Extract Baseline and Recovery for varName
    vBsl = tbl.(varName)(tbl.Day == 'BSL' & idxWt);
    vBac = tbl.(varName)(tbl.Day == 'BAC3' & idxWt);
    
    vBslMcu = tbl.(varName)(tbl.Day == 'BSL' & idxMcu);
    vBacMcu = tbl.(varName)(tbl.Day == 'BAC3' & idxMcu);
    
elseif strcmpi(dataSet, 'mea')
    
    % MEA Logic (Wide format)
    % Baseline is 'varName'
    % SteadyState is 'ss_varName'
    
    if strcmp(varName, 'fr')
        ssVarName = 'ss_fr';
    else
        ssVarName = ['ss_' varName];
    end
    
    if ~ismember(ssVarName, tbl.Properties.VariableNames)
        error('mcu_rcvQq:missingVar', 'MEA data requires "%s" column for variable "%s".', ssVarName, varName);
    end
    
    vBsl = tbl.(varName)(idxWt);
    vBac = tbl.(ssVarName)(idxWt); % Steady State acts as Recovery
    
    vBslMcu = tbl.(varName)(idxMcu);
    vBacMcu = tbl.(ssVarName)(idxMcu);
    
end


%% ========================================================================
%  PLOTTING
%  ========================================================================

figName = sprintf('Distribution Check: %s', varName);
hFig = figure('Position', [300 300 900 900], 'Color', 'w', 'Name', figName);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. Log-Normality Check (BSL Control)
% ------------------------------------
nexttile;
hP = probplot('lognormal', vBsl);
xlabel(sprintf('%s (Hz)', varName));
ylabel('Probability');
title(sprintf('Log-Normality Check (BSL Control)'));
grid on; axis square;   

% Store Plot Data
qq1.Name = "BSL Control (Log-Normal)";
qq1.X_Label = "Data (Hz)";
qq1.Y_Label = "Probability";
% Convert Z-scores to Probabilities
qq1.Data = [hP(1).XData(:), normcdf(hP(1).YData(:))];


% 2. Group Check (Control vs MCU)
% -------------------------------
nexttile;
hQ = qqplot(log(vBsl), log(vBslMcu));
xlabel(sprintf('Control BSL Log(%s)', varName));
ylabel(sprintf('MCU BSL Log(%s)', varName));
title('Group: Control vs MCU (BSL)');
grid on; axis square;
refline(1,0); 
hQ(1).LineWidth = 1;

% Store Plot Data
qq2.Name = "Group: Control vs MCU";
qq2.X_Label = "Control Log(Val)";
qq2.Y_Label = "MCU Log(Val)";
qq2.Data = [hQ(1).XData(:), hQ(1).YData(:)];


% 3. Stability Check (BSL vs Recovery Control)
% --------------------------------------------
nexttile;  hold on;
hQ = qqplot(log(vBsl), log(vBac));
xlabel(sprintf('BSL Log(%s) Quantiles', varName));
ylabel(sprintf('Recovery Log(%s) Quantiles', varName));
title('Stability: BSL vs Recovery (Control)');
grid on; axis square;
refline(1,0); 
hQ(1).LineWidth = 1;

% Store Plot Data
qq3.Name = "Stability: Control";
qq3.X_Label = "BSL Log(Val)";
qq3.Y_Label = "Recovery Log(Val)";
qq3.Data = [hQ(1).XData(:), hQ(1).YData(:)];


% 4. Stability Check (BSL vs Recovery MCU)
% ----------------------------------------
nexttile;
hQ = qqplot(log(vBslMcu), log(vBacMcu));
xlabel(sprintf('MCU BSL Log(%s)', varName));
ylabel(sprintf('MCU Recovery Log(%s)', varName));
title('Stability: BSL vs Recovery (MCU)');
grid on; axis square;
refline(1,0); 
hQ(1).LineWidth = 1;

% Store Plot Data
qq4.Name = "Stability: MCU";
qq4.X_Label = "BSL Log(Val)";
qq4.Y_Label = "Recovery Log(Val)";
qq4.Data = [hQ(1).XData(:), hQ(1).YData(:)];


% =========================================================================
% COMPILE OUTPUT
% =========================================================================
qqData = struct2table([qq1; qq2; qq3; qq4], 'AsArray', true);

% =========================================================================
% ADD CUTOFF LINES
% =========================================================================

% Calculate and plot cutoff value in log space (Theoretical)
% We used Z-score based on the fit, not the empirical percentile
mu = mean(log(vBsl));
sigma = std(log(vBsl));
logCutoff = mu + (cutoff_z * sigma);
cutoff_Hz = exp(logCutoff);

% Add Cutoff to Tile 1 (BSL Normality)
% nexttile(1); hold on;
% xline(cutoff_Hz, 'r--', 'Cutoff', 'LineWidth', 1.5, ...
%     'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');

% Add Cutoff to Tile 2 (Control vs MCU)
nexttile(2); hold on;
xline(logCutoff, 'r--', 'Cutoff', 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');
yline(logCutoff, 'r--', 'Cutoff', 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'bottom');

% Add Cutoff to Tile 3 (BSL vs Recovery Control)
nexttile(3); hold on;
xline(logCutoff, 'r--', 'Cutoff', 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');
    
% Add Cutoff to Tile 4 (MCU Stability)
nexttile(4); hold on;
xline(logCutoff, 'r--', 'Cutoff', 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal');

end
