% Graphics Script for ifrh.m
% -------------------------------------------------------------------------
% This script visualizes the simulation state.
% EXPECTS IN WORKSPACE:
%   - hFig (handle)
%   - pop (structure array)
%   - iTrial (current trial index)
%   - avg_out, avg_in, theta (state vectors)
% -------------------------------------------------------------------------

%% INITIALIZATION
if ~exist('hFig', 'var') || ~isvalid(hFig)

    % --- CONFIGURATION ---
    hFig = figure('Name', 'IfRH Simulation', ...
        'Color', 'w', 'Position', [50 50 1400 1000]);

    % --- COLOR MAPS ---
    % Define distinctive color maps for different populations/variables
    % E: Greens, I: Reds, Input: Blues, Theta: Purples
    maps = struct();
    maps.E = @(n) [linspace(0.4, 0, n)', linspace(0.9, 0.5, n)', linspace(0.4, 0.3, n)'];
    maps.I = @(n) [linspace(1, 0.5, n)', linspace(0.8, 0, n)', linspace(0.2, 0, n)'];
    maps.In = @(n) [linspace(0.4, 0, n)', linspace(0.4, 0, n)', linspace(1, 0.6, n)'];
    maps.Th = @(n) [linspace(0.8, 0.4, n)', linspace(0.4, 0, n)', linspace(0.9, 0.5, n)'];

    % --- LAYOUT SETUP (4x3 Grid) ---
    % Store handles in a structured way for easy updates
    H = struct();

    % 1. EXCITATORY RATES (Row 1)
    ax = subplot(4,3, [1 2]); hold(ax, 'on');
    H.axRatesE = ax;
    nPlot = length(pop(1).idxPlot);
    H.linesE = gobjects(nPlot, 1);
    cmap = maps.E(nPlot);
    for i = 1:nPlot
        H.linesE(i) = plot(ax, NaN, NaN, 'Color', [cmap(i,:) 0.6], 'LineWidth', 1.0);
    end
    H.meanE = plot(ax, NaN, NaN, 'g', 'LineWidth', 3, 'DisplayName', 'Mean');
    ylabel(ax, 'Firing Rate (Hz)'); title(ax, ['Excitatory (' pop(1).learningRule ')']); grid(ax, 'on');
    if strcmpi(pop(1).learningRule, 'output_centric')
        yline(ax, pop(1).target, 'k--', 'LineWidth', 2);
    end

    % Hist E
    ax = subplot(4,3,3); hold(ax, 'on');
    H.histE = histogram(ax, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(ax, 'Distribution (E)'); xlabel(ax, 'Hz');

    % 2. INHIBITORY RATES (Row 2)
    ax = subplot(4,3, [4 5]); hold(ax, 'on');
    H.axRatesI = ax;
    nPlot = length(pop(2).idxPlot);
    H.linesI = gobjects(nPlot, 1);
    cmap = maps.I(nPlot);
    for i = 1:nPlot
        H.linesI(i) = plot(ax, NaN, NaN, 'Color', [cmap(i,:) 0.6], 'LineWidth', 1.0);
    end
    H.meanI = plot(ax, NaN, NaN, 'r', 'LineWidth', 3, 'DisplayName', 'Mean');
    ylabel(ax, 'Firing Rate (Hz)'); title(ax, 'Inhibitory'); grid(ax, 'on');

    % Hist I
    ax = subplot(4,3,6); hold(ax, 'on');
    H.histI = histogram(ax, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(ax, 'Distribution (I)'); xlabel(ax, 'Hz');

    % 3. SENSED INPUT (Row 3)
    ax = subplot(4,3, [7 8]); hold(ax, 'on');
    H.axIn = ax;
    nPlot = length(pop(1).idxPlot); % Plot Input for E pop
    H.linesIn = gobjects(nPlot, 1);
    cmap = maps.In(nPlot);
    for i = 1:nPlot
        H.linesIn(i) = plot(ax, NaN, NaN, 'Color', [cmap(i,:) 0.6], 'LineWidth', 1.0);
    end
    H.meanIn = plot(ax, NaN, NaN, 'b', 'LineWidth', 3, 'DisplayName', 'Mean');
    ylabel(ax, 'Input Current'); title(ax, 'Sensed Drive (Exc + Ext)'); grid(ax, 'on');
    if strcmpi(pop(1).learningRule, 'input_centric')
        yline(ax, pop(1).target, 'k--', 'LineWidth', 2);
    end

    % Hist Input
    ax = subplot(4,3,9); hold(ax, 'on');
    H.histIn = histogram(ax, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(ax, 'Distribution (Input)'); xlabel(ax, 'a.u.');

    % 4. THRESHOLD (Row 4)
    ax = subplot(4,3, [10 11]); hold(ax, 'on');
    H.axTh = ax;
    nPlot = length(pop(1).idxPlot);
    H.linesTh = gobjects(nPlot, 1);
    cmap = maps.Th(nPlot);
    for i = 1:nPlot
        H.linesTh(i) = plot(ax, NaN, NaN, 'Color', [cmap(i,:) 0.6], 'LineWidth', 1.0);
    end
    H.meanTh = plot(ax, NaN, NaN, 'm', 'LineWidth', 3, 'DisplayName', 'Mean');
    ylabel(ax, '\theta'); xlabel(ax, 'Trial'); title(ax, 'Threshold'); grid(ax, 'on');

    % Hist Theta
    ax = subplot(4,3,12); hold(ax, 'on');
    H.histTh = histogram(ax, 'FaceColor', 'm', 'FaceAlpha', 0.5, 'NumBins', 50);
    title(ax, 'Distribution (\theta)'); xlabel(ax, 'a.u.');

    % --- PERTURBATION BARS ---
    axes_list = {H.axRatesE, H.axRatesI, H.axIn, H.axTh};
    for iPop = 1:length(pop)
        p = pop(iPop).perturb;
        if p(3) ~= 0
            if iPop==1, col=[0.3 0.6 0.3]; txt='Ext E'; else, col=[0.8 0.3 0.3]; txt='Ext I'; end
            for k=1:length(axes_list)
                ax = axes_list{k}; yl = ylim(ax);
                line(ax, [p(1), p(2)], [yl(1) yl(1)], 'Color', [0.7 0.7 0.7], 'LineWidth', 6);
                text(ax, mean(p(1:2)), yl(1), sprintf('+%0.1f %s', p(3), txt), ...
                    'Color',col, 'Vert','bottom', 'Horiz','center', 'FontSize',8, 'FontWeight','bold');
            end
        end
    end
end

%% UPDATE DYNAMICS
if exist('iTrial', 'var') && iTrial > 1 && exist('H', 'var')

    x_vec = 1:iTrial;

    % Update Excitatory (Pop 1)
    p = pop(1);
    for i = 1:length(H.linesE)
        set(H.linesE(i), 'XData', x_vec, 'YData', p.hist.allRate(1:iTrial, i));
    end
    set(H.meanE, 'XData', x_vec, 'YData', p.hist.meanRate(1:iTrial));
    if exist('avg_out','var'), H.histE.Data = avg_out(p.idx); end

    % Update Inhibitory (Pop 2)
    p = pop(2);
    for i = 1:length(H.linesI)
        set(H.linesI(i), 'XData', x_vec, 'YData', p.hist.allRate(1:iTrial, i));
    end
    set(H.meanI, 'XData', x_vec, 'YData', p.hist.meanRate(1:iTrial));
    if exist('avg_out','var'), H.histI.Data = avg_out(p.idx); end

    % Update Input (Pop 1 Sensed)
    p = pop(1);
    if ~strcmp(p.learningRule,'none')
        for i = 1:length(H.linesIn)
            % Note: we don't strictly track individual input histories in the simplified struct
            % but if we did, meaningful update would go here.
            % For now, just mean is available or we re-use rate/theta patterns if we added input hist.
            % The main script stores 'meanSensed'.
        end
        % Just update mean for input as typically individual inputs aren't stored in history struct effectively yet
        set(H.meanIn, 'XData', x_vec, 'YData', p.hist.meanSensed(1:iTrial));
        if exist('avg_in','var'), H.histIn.Data = avg_in(p.idx); end
    end

    % Update Threshold (Pop 1)
    p = pop(1);
    for i = 1:length(H.linesTh)
        set(H.linesTh(i), 'XData', x_vec, 'YData', p.hist.allTheta(1:iTrial, i));
    end
    set(H.meanTh, 'XData', x_vec, 'YData', p.hist.meanTheta(1:iTrial));
    if exist('theta','var'), H.histTh.Data = theta(p.idx); end

    % Scroll X-Axis
    maxx = max(nTrial, iTrial+20);
    xlim(H.axRatesE, [1 maxx]); xlim(H.axRatesI, [1 maxx]);
    xlim(H.axIn, [1 maxx]);     xlim(H.axTh, [1 maxx]);

    drawnow limitrate;
end
