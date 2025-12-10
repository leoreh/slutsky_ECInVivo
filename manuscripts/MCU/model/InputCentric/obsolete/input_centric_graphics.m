% Graphics Script for InputCentric_Demo.m
% -------------------------------------------------------------------------
% This script visualizes the 'pop' structure dynamics.
% It is designed to be called:
%   1. Once at the start (to initialize figures).
%   2. Inside the simulation loop (to update data).
%
% EXPECTS IN WORKSPACE:
%   - hFig (handle)
%   - pop (structure array)
%   - iTrial (current trial index)
%   - nTrial (total trials)
% -------------------------------------------------------------------------

%% INITIALIZATION
if ~exist('hFig', 'var') || ~isvalid(hFig)
    
    % --- CONFIGURATION ---
    hFig = figure('Name', 'Input-Centric Homeostasis Demo', ...
                  'Color', 'w', 'Position', [50 50 1400 1000]);
              
    % Extract Pointers for Readability
    % We assume pop(1)=Exc, pop(2)=Inh based on the model definition
    pE = pop(1); 
    pI = pop(2);
    
    % Plotting Indices (Ensure we only plot a subset to save FPS)
    idxPlotE = pE.idxPlot; 
    idxPlotI = pI.idxPlot;
    nPlotE = length(idxPlotE);
    nPlotI = length(idxPlotI);

    % --- COLOR MAPS (Replicating the exact "Old" look) ---
    % Excitatory (Greens/Teals)
    cmap_E = [linspace(0.4, 0, nPlotE)', linspace(0.9, 0.5, nPlotE)', linspace(0.4, 0.3, nPlotE)'];
    
    % Inhibitory (Reds/Oranges)
    cmap_I = [linspace(1, 0.5, nPlotI)', linspace(0.8, 0, nPlotI)', linspace(0.2, 0, nPlotI)'];
    
    % Input (Blues)
    cmap_Input = [linspace(0.4, 0, nPlotE)', linspace(0.4, 0, nPlotE)', linspace(1, 0.6, nPlotE)'];
    
    % Theta (Purples)
    cmap_Theta = [linspace(0.8, 0.4, nPlotE)', linspace(0.4, 0, nPlotE)', linspace(0.9, 0.5, nPlotE)'];

    % --- LAYOUT SETUP (4x3 Grid) ---
    
    % 1. EXCITATORY RATES (Row 1)
    axRatesE = subplot(4,3, [1 2]); hold(axRatesE, 'on');
    hLinesE_indiv = gobjects(nPlotE, 1);
    for i = 1:nPlotE
        hLinesE_indiv(i) = plot(axRatesE, NaN, NaN, 'Color', [cmap_E(i,:) 0.6], 'LineWidth', 1.0);
    end
    hLineE = plot(axRatesE, NaN, NaN, 'g', 'LineWidth', 3, 'DisplayName', 'Mean'); 
    ylabel(axRatesE, 'Firing Rate (Hz)'); 
    title(axRatesE, ['Excitatory Population (' pE.learningRule ')']); 
    grid(axRatesE, 'on');
    
    if strcmpi(pE.learningRule, 'output_centric')
        yline(axRatesE, pE.target, 'k--', 'LineWidth', 2, 'Label', 'Target');
    end

    % Histogram E (FaceAlpha 0.5 matches old style)
    axHistE = subplot(4,3,3); hold(axHistE, 'on');
    hHistObjE = histogram(axHistE, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistE, 'Distribution (E)'); xlabel(axHistE, 'Hz');


    % 2. INHIBITORY RATES (Row 2)
    axRatesI = subplot(4,3, [4 5]); hold(axRatesI, 'on');
    hLinesI_indiv = gobjects(nPlotI, 1);
    for i = 1:nPlotI
        hLinesI_indiv(i) = plot(axRatesI, NaN, NaN, 'Color', [cmap_I(i,:) 0.6], 'LineWidth', 1.0);
    end
    hLineI = plot(axRatesI, NaN, NaN, 'r', 'LineWidth', 3, 'DisplayName', 'Mean'); 
    ylabel(axRatesI, 'Firing Rate (Hz)'); 
    title(axRatesI, 'Inhibitory Population'); 
    grid(axRatesI, 'on');
    
    % Histogram I
    axHistI = subplot(4,3,6); hold(axHistI, 'on');
    hHistObjI = histogram(axHistI, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistI, 'Distribution (I)'); xlabel(axHistI, 'Hz');


    % 3. SENSED VARIABLE (Row 3)
    axInput = subplot(4,3, [7 8]); hold(axInput, 'on');
    hLinesInput_indiv = gobjects(nPlotE, 1);
    for i = 1:nPlotE
        hLinesInput_indiv(i) = plot(axInput, NaN, NaN, 'Color', [cmap_Input(i,:) 0.6], 'LineWidth', 1.0);
    end
    hLineInput = plot(axInput, NaN, NaN, 'b', 'LineWidth', 3, 'DisplayName', 'Mean'); 
    ylabel(axInput, 'Input Current (a.u.)');
    title(axInput, 'Sensed Excitatory Drive (Exc + Ext + Noise)', 'FontWeight', 'bold');
    grid(axInput, 'on');
    
    if strcmpi(pE.learningRule, 'input_centric')
        yline(axInput, pE.target, 'k--', 'LineWidth', 2, 'Label', 'Target');
    end

    % Histogram Input
    axHistInput = subplot(4,3,9); hold(axHistInput, 'on');
    hHistObjInput = histogram(axHistInput, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistInput, 'Distribution (Input)'); xlabel(axHistInput, 'a.u.');


    % 4. THRESHOLD (Theta) (Row 4)
    axTheta = subplot(4,3, [10 11]); hold(axTheta, 'on');
    hLinesTheta_indiv = gobjects(nPlotE, 1);
    for i = 1:nPlotE
        hLinesTheta_indiv(i) = plot(axTheta, NaN, NaN, 'Color', [cmap_Theta(i,:) 0.6], 'LineWidth', 1.0);
    end
    hLineTheta = plot(axTheta, NaN, NaN, 'm', 'LineWidth', 3, 'DisplayName', 'Mean'); 
    ylabel(axTheta, '\theta (Threshold)', 'Interpreter', 'tex');
    xlabel(axTheta, 'Trial Number');
    title(axTheta, 'Homeostatic Threshold Adaptation', 'FontWeight', 'bold');
    grid(axTheta, 'on');

    % Histogram Theta
    axHistTheta = subplot(4,3,12); hold(axHistTheta, 'on');
    hHistObjTheta = histogram(axHistTheta, 'FaceColor', 'm', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(axHistTheta, 'Distribution (\theta)', 'Interpreter', 'tex'); xlabel(axHistTheta, 'a.u.');


    % --- DRAW PERTURBATION INDICATORS ---
    % Loop through populations to draw lines where perturbations occur
    axes_list = {axRatesE, axRatesI, axInput, axTheta};
    
    for iPop = 1:length(pop)
        pert = pop(iPop).perturb; % [Onset, Offset, Amp]
        if pert(3) ~= 0 
            % Define color: Greenish for E, Reddish for I
            if iPop == 1, col = [0.3 0.6 0.3]; txt='Ext E'; 
            else, col = [0.8 0.3 0.3]; txt='Ext I'; end
            
            for k = 1:length(axes_list)
                ax = axes_list{k};
                yl = ylim(ax);
                % Draw the gray duration bar
                line(ax, [pert(1), pert(2)], [yl(1), yl(1)], ...
                     'Color', [0.7 0.7 0.7], 'LineWidth', 6);
                 
                % Draw the text label
                text(ax, mean([pert(1), pert(2)]), yl(1), ...
                     sprintf('+%0.1f %s', pert(3), txt), ...
                     'Color', col, 'FontWeight', 'bold', ...
                     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
                     'FontSize', 8, 'Interpreter', 'none');
            end
        end
    end
end


%% UPDATE DYNAMICS
if exist('iTrial', 'var') && iTrial > 0
    
    % Shortcuts
    x_vec = 1:iTrial;
    pE = pop(1);
    pI = pop(2);
    
    % 1. Update Lines (Excitatory)
    % ---------------------------------------------------------
    if iTrial > 1
        % Update Individual Traces
        for i = 1:nPlotE
            if isvalid(hLinesE_indiv(i))
                set(hLinesE_indiv(i), 'XData', x_vec, 'YData', pE.hist.allRate(1:iTrial, i)); 
            end
        end
        % Update Mean
        set(hLineE, 'XData', x_vec, 'YData', pE.hist.meanRate(1:iTrial));
        
        % Update Histograms: MUST USE TRIAL AVERAGE (avg_r), NOT INSTANTANEOUS
        % We access the 'avg_r' variable calculated in the main loop if available,
        % otherwise we might have to re-calculate it or use a stored field.
        % Assuming Main Loop defines 'avg_r' (local var) or we use the last history entry.
        
        % Robust Method: Use the history of the current trial
        % (Reshaping the last row of allRate isn't enough as it's a subset).
        % We will use the 'avg_r' variable from the main script workspace.
        if exist('avg_r_E', 'var'), hHistObjE.Data = avg_r_E; end
    end
    
    % 2. Update Lines (Inhibitory)
    % ---------------------------------------------------------
    if iTrial > 1
        for i = 1:nPlotI
            if isvalid(hLinesI_indiv(i))
                set(hLinesI_indiv(i), 'XData', x_vec, 'YData', pI.hist.allRate(1:iTrial, i)); 
            end
        end
        set(hLineI, 'XData', x_vec, 'YData', pI.hist.meanRate(1:iTrial));
        
        if exist('avg_r_I', 'var'), hHistObjI.Data = avg_r_I; end
    end
    
    % 3. Update Lines (Input)
    % ---------------------------------------------------------
    if iTrial > 1 && ~strcmpi(pE.learningRule, 'none')
        % Only update input history if we are tracking it
        set(hLineInput, 'XData', x_vec, 'YData', pE.hist.meanSensed(1:iTrial));
        
        % Input Histogram
        if exist('avg_sensed_E', 'var'), hHistObjInput.Data = avg_sensed_E; end
    end
    
    % 4. Update Lines (Theta)
    % ---------------------------------------------------------
    if iTrial > 1
        for i = 1:nPlotE
             if isvalid(hLinesTheta_indiv(i))
                set(hLinesTheta_indiv(i), 'XData', x_vec, 'YData', pE.hist.allTheta(1:iTrial, i)); 
            end
        end
        set(hLineTheta, 'XData', x_vec, 'YData', pE.hist.meanTheta(1:iTrial));
        
        % Theta changes slowly, so current state is fine
        hHistObjTheta.Data = pE.theta;
    end

    % --- AXIS SCROLLING ---
    % Smooth scrolling like the original
    max_x = max(nTrial, iTrial + 20);
    xlim(axRatesE, [1, max_x]);
    xlim(axRatesI, [1, max_x]);
    xlim(axInput,  [1, max_x]);
    xlim(axTheta,  [1, max_x]);
    
    % --- PERFORMANCE CRITICAL ---
    drawnow limitrate; % Prevents graphics from slowing down simulation
end