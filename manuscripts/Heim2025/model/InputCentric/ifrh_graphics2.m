% =========================================================================
%  IFRH GRAPHICS ENGINE
% =========================================================================
%  Visualizes simulation state for ifrh.m.
%  Expects workspace variables: hFig, pop, iTrial, avg_out, theta, etc.
% =========================================================================

%% 1. CONFIGURATION
%  ------------------------------------------------------------------------

% Viewport Settings
cfg.scroll_buffer = 100;       % How many trials to look ahead when scrolling
cfg.hist_bins     = 50;        % Number of bins for histograms

% Visual Styling
cfg.fig_pos       = [50 50 1400 1000];
cfg.color.E       = [0.2 0.7 0.2]; % Green
cfg.color.I       = [0.8 0.2 0.2]; % Red
cfg.color.In      = [0.2 0.4 0.8]; % Blue
cfg.color.Th      = [0.8 0.2 0.8]; % Purple
cfg.alpha_trace   = 0.3;           % Transparency for individual traces
cfg.alpha_hist    = 0.5;           % Transparency for histograms

% Color Generators (Gradient for individual lines)
% Generates a matrix of colors fading to black/white based on index
get_cmap = @(base, n) [linspace(base(1), base(1)*0.5, n)', ...
                       linspace(base(2), base(2)*0.5, n)', ...
                       linspace(base(3), base(3)*0.5, n)'];

%% 2. INITIALIZATION (Runs once)
%  ------------------------------------------------------------------------
if ~exist('hFig', 'var') || ~isvalid(hFig)

    % --- Setup Figure ---
    hFig = figure('Name', 'IFRH Simulation Dashboard', ...
                  'Color', 'w', 'Position', cfg.fig_pos);
    H = struct(); % Handle storage

    % --- A. EXCITATORY POPULATION (Row 1) ---
    % Time Series
    H.axRateE = subplot(4, 3, [1 2]); hold(H.axRateE, 'on');
    H.linesE  = init_lines(H.axRateE, pop(1), cfg.color.E, get_cmap, cfg.alpha_trace);
    ylabel(H.axRateE, 'Rate (Hz)'); title(H.axRateE, ['Excitatory (' pop(1).learningRule ')']);
    if strcmpi(pop(1).learningRule, 'output_centric')
        yline(H.axRateE, pop(1).target, 'k--', 'LineWidth', 1.5);
    end

    % Histogram (Pre-calculate edges for speed)
    H.axHistE = subplot(4, 3, 3); hold(H.axHistE, 'on');
    edges = linspace(0, pop(1).max_rate, cfg.hist_bins);
    H.barE = histogram(H.axHistE, 'BinEdges', edges, 'FaceColor', cfg.color.E, ...
                       'FaceAlpha', cfg.alpha_hist, 'EdgeColor', 'none');
    title(H.axHistE, 'Dist (E)'); xlabel(H.axHistE, 'Hz'); view(90, -90); % Rotate

    % --- B. INHIBITORY POPULATION (Row 2) ---
    % Time Series
    H.axRateI = subplot(4, 3, [4 5]); hold(H.axRateI, 'on');
    H.linesI  = init_lines(H.axRateI, pop(2), cfg.color.I, get_cmap, cfg.alpha_trace);
    ylabel(H.axRateI, 'Rate (Hz)'); title(H.axRateI, 'Inhibitory');

    % Histogram
    H.axHistI = subplot(4, 3, 6); hold(H.axHistI, 'on');
    edges = linspace(0, pop(2).max_rate, cfg.hist_bins);
    H.barI = histogram(H.axHistI, 'BinEdges', edges, 'FaceColor', cfg.color.I, ...
                       'FaceAlpha', cfg.alpha_hist, 'EdgeColor', 'none');
    title(H.axHistI, 'Dist (I)'); xlabel(H.axHistI, 'Hz'); view(90, -90);

    % --- C. SENSED INPUT (Row 3) ---
    % Time Series
    H.axInput = subplot(4, 3, [7 8]); hold(H.axInput, 'on');
    H.linesIn = init_lines(H.axInput, pop(1), cfg.color.In, get_cmap, cfg.alpha_trace);
    ylabel(H.axInput, 'Current (a.u.)'); title(H.axInput, 'Sensed Drive');
    if strcmpi(pop(1).learningRule, 'input_centric')
        yline(H.axInput, pop(1).target, 'k--', 'LineWidth', 1.5);
    end

    % Histogram
    H.axHistIn = subplot(4, 3, 9); hold(H.axHistIn, 'on');
    % Note: Input range is hard to predict, using automatic binning for Input only
    H.barIn = histogram(H.axHistIn, 'FaceColor', cfg.color.In, ...
                        'FaceAlpha', cfg.alpha_hist, 'NumBins', cfg.hist_bins);
    title(H.axHistIn, 'Dist (Input)'); view(90, -90);

    % --- D. THRESHOLDS (Row 4) ---
    % Time Series
    H.axTheta = subplot(4, 3, [10 11]); hold(H.axTheta, 'on');
    H.linesTh = init_lines(H.axTheta, pop(1), cfg.color.Th, get_cmap, cfg.alpha_trace);
    ylabel(H.axTheta, '\theta (mV)'); xlabel(H.axTheta, 'Trial'); title(H.axTheta, 'Thresholds');

    % Histogram
    H.axHistTh = subplot(4, 3, 12); hold(H.axHistTh, 'on');
    H.barTh = histogram(H.axHistTh, 'FaceColor', cfg.color.Th, ...
                        'FaceAlpha', cfg.alpha_hist, 'NumBins', cfg.hist_bins);
    title(H.axHistTh, 'Dist (\theta)'); view(90, -90);

    % --- E. COMMON FORMATTING & ANNOTATIONS ---
    all_axes = [H.axRateE, H.axRateI, H.axInput, H.axTheta];
    linkaxes(all_axes, 'x'); % Link X-zooming
    grid(all_axes, 'on');

    % Add Perturbation Bars
    for iPop = 1:length(pop)
        p = pop(iPop).perturb; % [Start, End, Amp]
        if p(3) ~= 0
            % Determine color/text
            if iPop==1, c=[0.2 0.5 0.2]; txt='Ext E'; else, c=[0.8 0.2 0.2]; txt='Ext I'; end
            
            % Draw on all time-series axes
            for ax = all_axes
                yl = ylim(ax);
                line(ax, [p(1) p(2)], [yl(1) yl(1)], 'Color', [0.6 0.6 0.6], 'LineWidth', 5);
                text(ax, mean(p(1:2)), yl(1), sprintf('+%.1f %s', p(3), txt), ...
                    'Color', c, 'VerticalAlignment','bottom', ...
                    'HorizontalAlignment','center', 'FontSize', 8, 'FontWeight','bold');
            end
        end
    end
end

%% 3. UPDATE ROUTINE (Runs every graphics_step)
%  ------------------------------------------------------------------------
if exist('iTrial', 'var') && iTrial > 1 && exist('H', 'var') && isvalid(hFig)

    x_vec = 1:iTrial;

    % --- A. Vectorized Line Updates (Performance Critical) ---
    % 1. Excitatory Rates
    update_lines(H.linesE, x_vec, pop(1).hist.allRate(1:iTrial, :), pop(1).hist.meanRate(1:iTrial));
    if exist('avg_out','var'), H.barE.Data = avg_out(pop(1).idx); end

    % 2. Inhibitory Rates
    update_lines(H.linesI, x_vec, pop(2).hist.allRate(1:iTrial, :), pop(2).hist.meanRate(1:iTrial));
    if exist('avg_out','var'), H.barI.Data = avg_out(pop(2).idx); end

    % 3. Input (Sensed)
    % Note: We typically only track mean sensed input in the main struct
    set(H.linesIn.mean, 'XData', x_vec, 'YData', pop(1).hist.meanSensed(1:iTrial));
    if exist('avg_in','var'), H.barIn.Data = avg_in(pop(1).idx); end

    % 4. Thresholds
    update_lines(H.linesTh, x_vec, pop(1).hist.allTheta(1:iTrial, :), pop(1).hist.meanTheta(1:iTrial));
    if exist('theta','var'), H.barTh.Data = theta(pop(1).idx); end


    % --- B. Smart Axis Scrolling ---
    % Only update xlim if the lines hit the edge of the current view
    current_lim = H.axRateE.XLim;
    if iTrial >= current_lim(2)
        new_max = iTrial + cfg.scroll_buffer;
        set([H.axRateE, H.axRateI, H.axInput, H.axTheta], 'XLim', [1 new_max]);
    end

    drawnow limitrate;
end

%% 4. HELPER FUNCTIONS (Logic encapsulation)
%  ------------------------------------------------------------------------
%  These standard logic blocks are used to initialize the handles struct
%  ------------------------------------------------------------------------

function h = init_lines(ax, pStruct, baseColor, mapFunc, alpha)
    % Initialize plot objects for individual neurons + mean
    nPlot = length(pStruct.idxPlot);
    colors = mapFunc(baseColor, nPlot);
    
    h.indiv = gobjects(nPlot, 1);
    for i = 1:nPlot
        h.indiv(i) = plot(ax, NaN, NaN, 'Color', [colors(i,:) alpha], 'LineWidth', 0.5);
    end
    h.mean = plot(ax, NaN, NaN, 'Color', baseColor, 'LineWidth', 2, 'DisplayName', 'Mean');
end

function update_lines(hStruct, x, y_data, y_mean)
    % Vectorized update of line objects (Avoids loops)
    if ~isempty(hStruct.indiv)
        % Convert data matrix to cell array for vectorized set()
        y_cells = num2cell(y_data, 1);
        set(hStruct.indiv, 'XData', x, {'YData'}, y_cells');
    end
    set(hStruct.mean, 'XData', x, 'YData', y_mean);
end