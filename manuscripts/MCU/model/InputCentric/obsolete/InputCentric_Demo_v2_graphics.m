%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS SCRIPT: Input-Centric Homeostasis Demo
% VERSION 2: Hybrid Architecture
%
% This script handles all visualization for the simulation.
% It is called:
%   1. Once at initialization (creates figure and axes)
%   2. Repeatedly during simulation (updates data)
%
% EXPECTS IN WORKSPACE:
%   - config: Configuration structure
%   - history: History structure with trial data
%   - E, I, thetaE: Current state variables
%   - avgE, avgI, avgInput: Current trial averages
%   - trial: Current trial number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION
if ~exist('hFig', 'var') || ~isvalid(hFig)
    
    % Create figure
    hFig = figure('Name', 'Input-Centric Homeostasis Demo', ...
                  'Color', 'w', 'Position', [50 50 1400 1000]);
    
    % Get configuration (with fallbacks for robustness)
    if ~exist('config', 'var')
        if exist('nE', 'var'), nE = nE; else, nE = 80; end
        if exist('nI', 'var'), nI = nI; else, nI = 20; end
        learning_rule = 'input_centric';
        SetpointMean = 20;
        PerturbRangeE = [150, 250];
        PerturbRangeI = [150, 250];
        PerturbAmpE = 5.0;
        PerturbAmpI = 5.0;
        nTrial = 500;
    else
        nE = config.nE;
        nI = config.nI;
        learning_rule = config.learning_rule;
        SetpointMean = config.E.target;
        PerturbRangeE = config.perturb.E(1:2);
        PerturbRangeI = config.perturb.I(1:2);
        PerturbAmpE = config.perturb.E(3);
        PerturbAmpI = config.perturb.I(3);
        nTrial = config.nTrial;
    end
    
    % Get plotting indices
    if ~exist('idxPlotE', 'var')
        idxPlotE = round(linspace(1, nE, min(nE, 10)));
    end
    if ~exist('idxPlotI', 'var')
        idxPlotI = round(linspace(1, nI, min(nI, 5)));
    end
    nPlotE = length(idxPlotE);
    nPlotI = length(idxPlotI);
    
    % ---------------------------------------------------------------------
    % Color Maps (Gradient families for visual separation)
    % ---------------------------------------------------------------------
    % Excitatory (Greens/Teals)
    cmap_E = [linspace(0.4, 0, nPlotE)', ...
               linspace(0.9, 0.5, nPlotE)', ...
               linspace(0.4, 0.3, nPlotE)'];
    
    % Inhibitory (Reds/Oranges)
    cmap_I = [linspace(1, 0.5, nPlotI)', ...
               linspace(0.8, 0, nPlotI)', ...
               linspace(0.2, 0, nPlotI)'];
    
    % Input (Blues)
    cmap_Input = [linspace(0.4, 0, nPlotE)', ...
                   linspace(0.4, 0, nPlotE)', ...
                   linspace(1, 0.6, nPlotE)'];
    
    % Theta (Purples)
    cmap_Theta = [linspace(0.8, 0.4, nPlotE)', ...
                   linspace(0.4, 0, nPlotE)', ...
                   linspace(0.9, 0.5, nPlotE)'];
    
    % ---------------------------------------------------------------------
    % LAYOUT: 4x3 Grid
    % ---------------------------------------------------------------------
    % Column 1-2: Time series
    % Column 3: Histograms
    % Row 1: Excitatory Rates
    % Row 2: Inhibitory Rates
    % Row 3: Synaptic Input (Regulated Variable)
    % Row 4: Threshold (Theta)
    
    % =====================================================================
    % ROW 1: EXCITATORY RATES
    % =====================================================================
    axRatesE = subplot(4, 3, [1 2]);
    hold(axRatesE, 'on');
    
    % Individual neuron traces
    hLinesE_indiv = gobjects(nPlotE, 1);
    for i = 1:nPlotE
        hLinesE_indiv(i) = plot(axRatesE, NaN, NaN, ...
                               'Color', [cmap_E(i,:) 0.7], ...
                               'LineWidth', 1.0);
    end
    
    % Mean trace
    hLineE = plot(axRatesE, NaN, NaN, 'g', 'LineWidth', 3, ...
                  'DisplayName', 'Mean');
    
    ylabel(axRatesE, 'E Rate (Hz)');
    title(axRatesE, 'Excitatory Population');
    grid(axRatesE, 'on');
    
    % Draw setpoint if output-centric
    if strcmp(learning_rule, 'output_centric')
        yline(axRatesE, SetpointMean, 'k--', 'LineWidth', 2, ...
              'Label', 'Setpoint');
        title(axRatesE, 'Excitatory Population (Regulated Output)');
    end
    
    % Histogram
    axHistE = subplot(4, 3, 3);
    hold(axHistE, 'on');
    hHistObjE = histogram(axHistE, 'FaceColor', 'g', ...
                          'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistE, 'Distribution (E)');
    xlabel(axHistE, 'Hz');
    
    % =====================================================================
    % ROW 2: INHIBITORY RATES
    % =====================================================================
    axRatesI = subplot(4, 3, [4 5]);
    hold(axRatesI, 'on');
    
    % Individual neuron traces
    hLinesI_indiv = gobjects(nPlotI, 1);
    for i = 1:nPlotI
        hLinesI_indiv(i) = plot(axRatesI, NaN, NaN, ...
                               'Color', [cmap_I(i,:) 0.7], ...
                               'LineWidth', 1.0);
    end
    
    % Mean trace
    hLineI = plot(axRatesI, NaN, NaN, 'r', 'LineWidth', 3, ...
                  'DisplayName', 'Mean');
    
    ylabel(axRatesI, 'I Rate (Hz)');
    title(axRatesI, 'Inhibitory Population');
    grid(axRatesI, 'on');
    
    % Histogram
    axHistI = subplot(4, 3, 6);
    hold(axHistI, 'on');
    hHistObjI = histogram(axHistI, 'FaceColor', 'r', ...
                          'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistI, 'Distribution (I)');
    xlabel(axHistI, 'Hz');
    
    % =====================================================================
    % ROW 3: SYNAPTIC INPUT (Regulated Variable)
    % =====================================================================
    axInput = subplot(4, 3, [7 8]);
    hold(axInput, 'on');
    
    % Individual neuron traces
    hLinesInput_indiv = gobjects(nPlotE, 1);
    for i = 1:nPlotE
        hLinesInput_indiv(i) = plot(axInput, NaN, NaN, ...
                                    'Color', [cmap_Input(i,:) 0.7], ...
                                    'LineWidth', 1.0);
    end
    
    % Mean trace
    hLineInput = plot(axInput, NaN, NaN, 'b', 'LineWidth', 3, ...
                      'DisplayName', 'Mean');
    
    ylabel(axInput, 'Input (a.u.)');
    grid(axInput, 'on');
    
    % Draw setpoint if input-centric
    if strcmp(learning_rule, 'input_centric')
        yline(axInput, SetpointMean, 'k--', 'LineWidth', 2, ...
              'Label', 'Setpoint');
        title(axInput, 'Regulated Input (Exc + Ext)', 'FontWeight', 'bold');
    else
        title(axInput, 'Synaptic Input (Exc + Ext)', 'FontWeight', 'bold');
    end
    
    % Histogram
    axHistInput = subplot(4, 3, 9);
    hold(axHistInput, 'on');
    hHistObjInput = histogram(axHistInput, 'FaceColor', 'b', ...
                              'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistInput, 'Distribution (Input)');
    xlabel(axHistInput, 'a.u.');
    
    % =====================================================================
    % ROW 4: THRESHOLD (Theta)
    % =====================================================================
    axTheta = subplot(4, 3, [10 11]);
    hold(axTheta, 'on');
    
    % Individual neuron traces
    hLinesTheta_indiv = gobjects(nPlotE, 1);
    for i = 1:nPlotE
        hLinesTheta_indiv(i) = plot(axTheta, NaN, NaN, ...
                                   'Color', [cmap_Theta(i,:) 0.7], ...
                                   'LineWidth', 1.0);
    end
    
    % Mean trace
    hLineTheta = plot(axTheta, NaN, NaN, 'm', 'LineWidth', 3, ...
                      'DisplayName', 'Mean');
    
    ylabel(axTheta, '\theta_E', 'Interpreter', 'tex');
    xlabel(axTheta, 'Trial');
    title(axTheta, 'Firing Threshold', 'FontWeight', 'bold');
    grid(axTheta, 'on');
    
    % Histogram
    axHistTheta = subplot(4, 3, 12);
    hold(axHistTheta, 'on');
    hHistObjTheta = histogram(axHistTheta, 'FaceColor', 'm', ...
                              'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistTheta, 'Distribution (\theta)', 'Interpreter', 'tex');
    xlabel(axHistTheta, 'a.u.');
    
    % ---------------------------------------------------------------------
    % PERTURBATION INDICATORS
    % ---------------------------------------------------------------------
    % Draw gray bars and labels showing when perturbations occur
    
    % Excitatory perturbations (on E-related plots)
    axes_E = [axRatesE, axInput, axTheta];
    if PerturbAmpE ~= 0
        for ax = axes_E
            yl = ylim(ax);
            line(ax, PerturbRangeE, [yl(1), yl(1)], ...
                 'Color', [0.6, 0.6, 0.6], 'LineWidth', 6);
            text(ax, mean(PerturbRangeE), yl(1), ...
                 sprintf('+%g Ext E', PerturbAmpE), ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', ...
                 'FontWeight', 'bold', 'Color', [0.3 0.6 0.3], ...
                 'FontSize', 8, 'Interpreter', 'none');
        end
    end
    
    % Inhibitory perturbations (on I-related plots)
    axes_I = [axRatesI];
    if PerturbAmpI ~= 0
        for ax = axes_I
            yl = ylim(ax);
            line(ax, PerturbRangeI, [yl(1), yl(1)], ...
                 'Color', [0.6, 0.6, 0.6], 'LineWidth', 6);
            text(ax, mean(PerturbRangeI), yl(1), ...
                 sprintf('+%g Ext I', PerturbAmpI), ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', ...
                 'FontWeight', 'bold', 'Color', [0.8 0.3 0.3], ...
                 'FontSize', 8, 'Interpreter', 'none');
        end
    end
    
end

%% UPDATE LOGIC
if exist('trial', 'var') && trial > 0
    
    x_vec = 1:trial;
    
    % ---------------------------------------------------------------------
    % 1. Excitatory Rates
    % ---------------------------------------------------------------------
    if trial > 1
        for i = 1:nPlotE
            if isvalid(hLinesE_indiv(i))
                set(hLinesE_indiv(i), ...
                    'XData', x_vec, ...
                    'YData', history.allE(1:trial, i));
            end
        end
    end
    set(hLineE, 'XData', x_vec, 'YData', history.meanE(1:trial));
    xlim(axRatesE, [1, max(nTrial, trial+10)]);
    
    % Update histogram
    if exist('avgE', 'var')
        hHistObjE.Data = avgE;
    end
    
    % ---------------------------------------------------------------------
    % 2. Inhibitory Rates
    % ---------------------------------------------------------------------
    if trial > 1
        for i = 1:nPlotI
            if isvalid(hLinesI_indiv(i))
                set(hLinesI_indiv(i), ...
                    'XData', x_vec, ...
                    'YData', history.allI(1:trial, i));
            end
        end
    end
    set(hLineI, 'XData', x_vec, 'YData', history.meanI(1:trial));
    xlim(axRatesI, [1, max(nTrial, trial+10)]);
    
    % Update histogram
    if exist('avgI', 'var')
        hHistObjI.Data = avgI;
    end
    
    % ---------------------------------------------------------------------
    % 3. Synaptic Input
    % ---------------------------------------------------------------------
    if trial > 1
        for i = 1:nPlotE
            if isvalid(hLinesInput_indiv(i))
                set(hLinesInput_indiv(i), ...
                    'XData', x_vec, ...
                    'YData', history.allInput(1:trial, i));
            end
        end
    end
    set(hLineInput, 'XData', x_vec, 'YData', history.meanInput(1:trial));
    xlim(axInput, [1, max(nTrial, trial+10)]);
    
    % Update histogram
    if exist('avgInput', 'var')
        hHistObjInput.Data = avgInput;
    end
    
    % ---------------------------------------------------------------------
    % 4. Threshold (Theta)
    % ---------------------------------------------------------------------
    if trial > 1
        for i = 1:nPlotE
            if isvalid(hLinesTheta_indiv(i))
                set(hLinesTheta_indiv(i), ...
                    'XData', x_vec, ...
                    'YData', history.allTheta(1:trial, i));
            end
        end
    end
    set(hLineTheta, 'XData', x_vec, 'YData', history.meanTheta(1:trial));
    xlim(axTheta, [1, max(nTrial, trial+10)]);
    
    % Update histogram
    if exist('thetaE', 'var')
        hHistObjTheta.Data = thetaE;
    end
    
    % Force graphics update
    drawnow limitrate;
    
end

