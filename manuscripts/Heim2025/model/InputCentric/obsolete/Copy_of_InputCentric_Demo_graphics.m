% Graphics Script for InputCentric_Demo.m
% This script initializes the figure and defines the update logic.

if ~exist('hFig', 'var') || ~isvalid(hFig)
    % --- INITIALIZATION ---
    hFig = figure('Name', 'Input-Centric Homeostasis Demo', 'Color', 'w', 'Position', [50 50 1400 1000]);

    % Requires nE, nI, PerturbRangeE, PerturbRangeI in workspace
    if ~exist('nE','var'), nE=80; end
    if ~exist('nI','var'), nI=20; end
    if ~exist('PerturbAmp','var'), PerturbAmp=4; end

    % Robustly handling Perturbation Ranges
    if ~exist('PerturbRangeE','var')
        if exist('PerturbOnsetTrial','var'), PerturbRangeE=[PerturbOnsetTrial, 300]; else, PerturbRangeE=[150, 300]; end
    end
    if ~exist('PerturbRangeI','var')
        if exist('PerturbOnsetTrial','var'), PerturbRangeI=[PerturbOnsetTrial, 300]; else, PerturbRangeI=[150, 300]; end
    end

    % Robustly handling Perturbation Amplitudes
    if ~exist('PerturbAmpE','var')
        if exist('PerturbAmp','var'), PerturbAmpE = PerturbAmp; else, PerturbAmpE = 0; end
    end
    if ~exist('PerturbAmpI','var'), PerturbAmpI = 0; end

    % Handle Learning Rule Context
    if ~exist('learning_rule','var'), learning_rule='input_centric'; end
    if ~exist('HomeoSetpoint','var')
        if exist('InputSet','var'), HomeoSetpoint=ones(nE,1)*InputSet; else, HomeoSetpoint=zeros(nE,1); end
    end
    SetpointMean = mean(HomeoSetpoint);


    % Ensure nPlot and indices exist (fallback)
    if ~exist('nPlot','var'), nPlot=min(nE, 5); end
    if ~exist('idxPlotE','var'), idxPlotE=round(linspace(1, nE, nPlot)); end
    if ~exist('idxPlotI','var'), idxPlotI=round(linspace(1, nI, nPlot)); end
    nPlotE = length(idxPlotE);
    nPlotI = length(idxPlotI);

    % --- LAYOUT ---
    % Grid 4x3.
    % Col 1-2: Timeseries. Col 3: Histogram.
    % Row 1: E Rates
    % Row 2: I Rates
    % Row 3: Input
    % Row 4: Theta

    % Helper for Colors (Gradient Families with High Variance)
    % Using 'parula' or built-in maps provides better separation than simple linear gradients

    % Excitatory (Greens/Teals)
    % We construct a custom map varying from bright yellow-green to dark teal
    cmap_E = [linspace(0.4, 0, nPlotE)', linspace(0.9, 0.5, nPlotE)', linspace(0.4, 0.3, nPlotE)'];

    % Inhibitory (Reds/Oranges)
    % Variation from yellow-orange to dark red
    cmap_I = [linspace(1, 0.5, nPlotI)', linspace(0.8, 0, nPlotI)', linspace(0.2, 0, nPlotI)'];

    % Input (Blues)
    cmap_Input = [linspace(0.4, 0, nPlotE)', linspace(0.4, 0, nPlotE)', linspace(1, 0.6, nPlotE)'];

    % Theta (Purples)
    cmap_Theta = [linspace(0.8, 0.4, nPlotE)', linspace(0.4, 0, nPlotE)', linspace(0.9, 0.5, nPlotE)'];


    % =========================================================================
    % ROW 1: EXCITATORY RATES
    % =========================================================================
    axRatesE = subplot(4,3, [1 2]); hold(axRatesE, 'on');
    hLinesE_indiv = gobjects(nPlotE, 1);
    for i = 1:nPlotE
        hLinesE_indiv(i) = plot(axRatesE, NaN, NaN, 'Color', [cmap_E(i,:) 0.7], 'LineWidth', 1.0);
    end
    hLineE = plot(axRatesE, NaN, NaN, 'g', 'LineWidth', 3, 'DisplayName', 'Mean'); % Thicker Green
    ylabel(axRatesE, 'E Rate (Hz)'); title(axRatesE, 'Excitatory Populations'); grid(axRatesE, 'on');

    % If Output Centric, draw Setpoint here
    if strcmp(learning_rule, 'output_centric')
        hLineSetRate = line(axRatesE, [1 nTrial], [SetpointMean SetpointMean], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Setpoint');
        title(axRatesE, 'Excitatory Populations (Regulated Output)');
    end

    % Fix: Use NumBins instead of BinWidth. This constrains the bin count (safe for low var)
    % while allowing range to expand dynamically.
    axHistE = subplot(4,3,3); hold(axHistE, 'on');
    hHistObjE = histogram(axHistE, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistE, 'Distribution (E)'); xlabel(axHistE, 'Hz');

    % =========================================================================
    % ROW 2: INHIBITORY RATES
    % =========================================================================
    axRatesI = subplot(4,3, [4 5]); hold(axRatesI, 'on');
    hLinesI_indiv = gobjects(nPlotI, 1);
    for i = 1:nPlotI
        hLinesI_indiv(i) = plot(axRatesI, NaN, NaN, 'Color', [cmap_I(i,:) 0.7], 'LineWidth', 1.0);
    end
    hLineI = plot(axRatesI, NaN, NaN, 'r', 'LineWidth', 3, 'DisplayName', 'Mean'); % Thicker Red
    ylabel(axRatesI, 'I Rate (Hz)'); title(axRatesI, 'Inhibitory Populations'); grid(axRatesI, 'on');

    % Fix: NumBins
    axHistI = subplot(4,3,6); hold(axHistI, 'on');
    hHistObjI = histogram(axHistI, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistI, 'Distribution (I)'); xlabel(axHistI, 'Hz');

    % =========================================================================
    % ROW 3: SYNAPTIC INPUT (Regulated Variable)
    % =========================================================================
    axInput = subplot(4,3, [7 8]); hold(axInput, 'on');
    hLinesInput_indiv = gobjects(nPlotE, 1);
    for i = 1:nPlotE
        hLinesInput_indiv(i) = plot(axInput, NaN, NaN, 'Color', [cmap_Input(i,:) 0.7], 'LineWidth', 1.0);
    end
    hLineInput = plot(axInput, NaN, NaN, 'b', 'LineWidth', 3, 'DisplayName', 'Mean'); % Thicker Blue

    % If Input Centric, draw Setpoint here
    if strcmp(learning_rule, 'input_centric')
        hLineSetInput = line(axInput, [1 nTrial], [SetpointMean SetpointMean], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Setpoint');
        title(axInput, 'Regulated Input (Exc + Ext)', 'FontWeight', 'bold');
    else
        % Generic title if not regulated
        title(axInput, 'Synaptic Input (Exc + Ext)', 'FontWeight', 'bold');
    end

    % Updated Logic: tempInput is W_EE*E + Ext + Noise (Exc Drive only)
    ylabel(axInput, 'Input (a.u.)');
    grid(axInput, 'on');

    % Fix: NumBins
    axHistInput = subplot(4,3,9); hold(axHistInput, 'on');
    hHistObjInput = histogram(axHistInput, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistInput, 'Distribution (Input)'); xlabel(axHistInput, 'a.u.');

    % =========================================================================
    % ROW 4: THETA (THRESHOLD)
    % =========================================================================
    axTheta = subplot(4,3, [10 11]); hold(axTheta, 'on');
    hLinesTheta_indiv = gobjects(nPlotE, 1);
    for i = 1:nPlotE
        hLinesTheta_indiv(i) = plot(axTheta, NaN, NaN, 'Color', [cmap_Theta(i,:) 0.7], 'LineWidth', 1.0);
    end
    hLineTheta = plot(axTheta, NaN, NaN, 'm', 'LineWidth', 3, 'DisplayName', 'Mean'); % Thicker Purple

    % Use TeX interpreter for greek theta
    ylabel(axTheta, '\theta_E', 'Interpreter', 'tex');
    xlabel(axTheta, 'Trial');
    title(axTheta, 'Firing Threshold', 'FontWeight', 'bold'); % Renamed title
    grid(axTheta, 'on');

    % Fix: NumBins
    axHistTheta = subplot(4,3,12); hold(axHistTheta, 'on');
    hHistObjTheta = histogram(axHistTheta, 'FaceColor', 'm', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistTheta, 'Distribution (\theta)', 'Interpreter', 'tex'); xlabel(axHistTheta, 'a.u.');


    % --- STATIC PERTURBATION INDICATOR ---
    % Separate Logic for E and I inputs using RANGES [Onset, Offset]

    % 1. Graphs relating to Excitatory Population or E-Inputs
    axes_E = [axRatesE, axInput, axTheta];
    if PerturbAmpE ~= 0
        for ax = axes_E
            y_pos = 0;
            line(ax, PerturbRangeE, [y_pos, y_pos], 'Color', [0.6, 0.6, 0.6], 'LineWidth', 6);
            text(ax, mean(PerturbRangeE), y_pos, sprintf('+%g Ext E', PerturbAmpE), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontWeight', 'bold', 'Color', [0.3 0.6 0.3], 'FontSize', 8, 'Interpreter', 'none');
        end
    end

    % 2. Graphs relating to Inhibitory Population
    axes_I = [axRatesI];
    if PerturbAmpI ~= 0
        for ax = axes_I
            y_pos = 0;
            line(ax, PerturbRangeI, [y_pos, y_pos], 'Color', [0.6, 0.6, 0.6], 'LineWidth', 6);
            text(ax, mean(PerturbRangeI), y_pos, sprintf('+%g Ext I', PerturbAmpI), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontWeight', 'bold', 'Color', [0.8 0.3 0.3], 'FontSize', 8, 'Interpreter', 'none');
        end
    end

end

% --- UPDATE LOGIC ---
if exist('trial', 'var') && trial > 0
    x_vec = 1:trial;

    % 1. Excitatory Rates
    if trial > 1
        for i = 1:nPlotE
            if isvalid(hLinesE_indiv(i)), set(hLinesE_indiv(i), 'XData', x_vec, 'YData', history.allE(1:trial, idxPlotE(i))); end
        end
    end
    set(hLineE, 'XData', x_vec, 'YData', history.meanE(1:trial));
    xlim(axRatesE, [1, max(nTrial, trial+10)]);

    % 2. Inhibitory Rates
    if trial > 1
        for i = 1:nPlotI
            if isvalid(hLinesI_indiv(i)), set(hLinesI_indiv(i), 'XData', x_vec, 'YData', history.allI(1:trial, idxPlotI(i))); end
        end
    end
    set(hLineI, 'XData', x_vec, 'YData', history.meanI(1:trial));
    xlim(axRatesI, [1, max(nTrial, trial+10)]);

    % 3. Input
    if trial > 1
        for i = 1:nPlotE
            if isvalid(hLinesInput_indiv(i)), set(hLinesInput_indiv(i), 'XData', x_vec, 'YData', history.allInput(1:trial, idxPlotE(i))); end
        end
    end
    set(hLineInput, 'XData', x_vec, 'YData', history.meanInput(1:trial));
    xlim(axInput, [1, max(nTrial, trial+10)]);

    % 4. Theta
    if trial > 1
        for i = 1:nPlotE
            if isvalid(hLinesTheta_indiv(i)), set(hLinesTheta_indiv(i), 'XData', x_vec, 'YData', history.allTheta(1:trial, idxPlotE(i))); end
        end
    end
    set(hLineTheta, 'XData', x_vec, 'YData', history.meanTheta(1:trial));
    xlim(axTheta, [1, max(nTrial, trial+10)]);


    % --- UPDATE HISTOGRAMS ---
    % Use normal data, NumBins handles the scaling safely.
    if exist('avgE', 'var'), hHistObjE.Data = avgE; end
    if exist('avgI_vec', 'var'), hHistObjI.Data = avgI_vec; end
    if exist('avgInput', 'var'), hHistObjInput.Data = avgInput; end
    if exist('thetaE', 'var'), hHistObjTheta.Data = thetaE; end

    drawnow limitrate;
end
