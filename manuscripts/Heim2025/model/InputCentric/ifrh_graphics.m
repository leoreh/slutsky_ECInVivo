%% ========================================================================
%  IFRH GRAPHICS ENGINE (TABBED VERSION)
%  ========================================================================

% CONFIGURATION
cfg.hist_bins     = 50;
cfg.y_padding     = 0.1;
cfg.fig_pos       = [50 50 1000 900]; % Slightly narrower, vertical layout
cfg.alpha_trace   = 0.4;
cfg.alpha_hist    = 0.5;

% Color Palettes for different populations (Extendable)
pop_colors = {[0.2 0.7 0.2], [0.8 0.2 0.2], [0.2 0.2 0.2]}; % Green, Red, Black

% Color Map Generator
get_cmap = @(base, n) [linspace(min(1, base(1)*1.5), base(1)*0.3, n)', ...
    linspace(min(1, base(2)*1.5), base(2)*0.3, n)', ...
    linspace(min(1, base(3)*1.5), base(3)*0.3, n)'];

%% ========================================================================
%  INITIALIZATION (Runs once)
%  ========================================================================
if ~exist('hFig', 'var') || ~isvalid(hFig)

    hFig = figure('Name', 'IfRH Simulation Dashboard', ...
        'Color', 'w', 'Position', cfg.fig_pos, ...
        'NumberTitle', 'off');

    % CREATE TABS
    tabGroup = uitabgroup(hFig);
    H = struct(); % Handle storage array

    % Loop through populations to create dynamic tabs
    for i = 1:length(pop)

        % 1. Create Tab
        H(i).tab = uitab(tabGroup, 'Title', pop(i).name);

        % 2. Layout for this tab (4 Rows, 3 Columns)
        %    Columns 1-2: Time Series | Column 3: Histogram
        t = tiledlayout(H(i).tab, 4, 3, 'TileSpacing', 'tight', 'Padding', 'compact');

        % Select Base Color
        cBase = pop_colors{min(i, length(pop_colors))};

        % --- PANEL 1: RATE ---
        H(i).axRate = nexttile(t, [1 2]); hold(H(i).axRate, 'on');
        H(i).lRate  = init_lines(H(i).axRate, pop(i), cBase, get_cmap, cfg.alpha_trace);
        ylabel(H(i).axRate, 'Rate (Hz)'); title(H(i).axRate, 'Population Activity');
        % Master Title for the Tab (Learning Rule)
        ruleName = strrep(pop(i).learningRule, '_', ' ');
        title(t, ['Learning Rule: ' upper(ruleName)], 'FontSize', 14, 'FontWeight', 'bold');

        % Target Line (Rate)
        if strcmpi(pop(i).learningRule, 'output_centric')
            yline(H(i).axRate, mean(pop(i).target), 'k--', 'Target');
        end

        H(i).axHistR = nexttile(t); hold(H(i).axHistR, 'on');
        H(i).barR = histogram(H(i).axHistR, 'Orientation', 'horizontal', ...
            'FaceColor', cBase, 'FaceAlpha', cfg.alpha_hist, 'EdgeColor', 'none', 'NumBins', cfg.hist_bins);
        H(i).axHistR.YAxis.Visible = 'off';

        % --- PANEL 2: SENSED INPUT (Drive) ---
        H(i).axIn = nexttile(t, [1 2]); hold(H(i).axIn, 'on');
        H(i).lIn  = init_lines(H(i).axIn, pop(i), [0.2 0.4 0.8], get_cmap, cfg.alpha_trace); % Blue
        ylabel(H(i).axIn, 'Input (a.u.)'); title(H(i).axIn, 'Excitatory Drive');
        % Target Line (Input)
        if strcmpi(pop(i).learningRule, 'input_centric')
            yline(H(i).axIn, mean(pop(i).target), 'k--', 'Target');
        end

        H(i).axHistIn = nexttile(t); hold(H(i).axHistIn, 'on');
        H(i).barIn = histogram(H(i).axHistIn, 'Orientation', 'horizontal', ...
            'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', cfg.alpha_hist, 'EdgeColor', 'none', 'NumBins', cfg.hist_bins);
        H(i).axHistIn.YAxis.Visible = 'off';

        % --- PANEL 3: THRESHOLDS ---
        H(i).axTh = nexttile(t, [1 2]); hold(H(i).axTh, 'on');
        H(i).lTh  = init_lines(H(i).axTh, pop(i), [0.8 0.2 0.8], get_cmap, cfg.alpha_trace); % Purple
        ylabel(H(i).axTh, '\theta (mV)', 'Interpreter', 'tex'); title(H(i).axTh, 'Thresholds');

        H(i).axHistTh = nexttile(t); hold(H(i).axHistTh, 'on');
        H(i).barTh = histogram(H(i).axHistTh, 'Orientation', 'horizontal', ...
            'FaceColor', [0.8 0.2 0.8], 'FaceAlpha', cfg.alpha_hist, 'EdgeColor', 'none', 'NumBins', cfg.hist_bins);
        H(i).axHistTh.YAxis.Visible = 'off';

        % --- PANEL 4: ERROR ---
        if ~strcmpi(pop(i).learningRule, 'none')
            H(i).axErr = nexttile(t, [1 2]); hold(H(i).axErr, 'on');
            % Plot individual lines + mean
            H(i).lErr = init_lines(H(i).axErr, pop(i), [0.5 0.5 0.5], get_cmap, cfg.alpha_trace);
            yline(H(i).axErr, 0, 'k-'); % Zero line
            ylabel(H(i).axErr, 'Error'); title(H(i).axErr, 'Homeostatic Error (Sensed - Target)');
            xlabel(H(i).axErr, 'Trial');

            H(i).axHistErr = nexttile(t); hold(H(i).axHistErr, 'on');
            H(i).barErr = histogram(H(i).axHistErr, 'Orientation', 'horizontal', ...
                'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', cfg.alpha_hist, 'EdgeColor', 'none', 'NumBins', cfg.hist_bins);
            H(i).axHistErr.YAxis.Visible = 'off';

            % Update linkage
            linkaxes([H(i).axRate, H(i).axIn, H(i).axTh, H(i).axErr], 'x');
            linkaxes([H(i).axErr, H(i).axHistErr], 'y');
            grid(H(i).axErr, 'on');
        else
            % Fill space if no error plot
            nexttile(t, [1 2]); axis off;
            nexttile(t); axis off;
            linkaxes([H(i).axRate, H(i).axIn, H(i).axTh], 'x');
        end

        % --- LINKING & FORMATTING ---
        linkaxes([H(i).axRate, H(i).axHistR], 'y');
        linkaxes([H(i).axIn,   H(i).axHistIn], 'y');
        linkaxes([H(i).axTh,   H(i).axHistTh], 'y');

        xlim(H(i).axRate, [1 max(100, length(pop(1).hist.r))]);
        grid(H(i).axRate, 'on'); grid(H(i).axIn, 'on');
        grid(H(i).axTh, 'on');
        linkaxes([H(i).axRate, H(i).axIn, H(i).axTh, H(i).axErr], 'x');

        % Add Perturbation Bar (only on Rate plot)
        draw_perturbation_bar(H(i).axRate, pop(i), cBase);
    end
end

%% ========================================================================
%  UPDATE ROUTINE
%  ========================================================================
if exist('iTrial', 'var') && iTrial > 1 && exist('H', 'var') && isvalid(hFig)

    x_vec = 1:iTrial;

    % Loop through populations and update the graphics objects in H(i)
    for i = 1:length(pop)

        % 1. RATE
        [sub_R, mean_R] = get_plot_data(pop(i).hist.r, iTrial, pop(i).idxPlot);
        update_lines(H(i).lRate, x_vec, sub_R, mean_R);
        H(i).barR.Data = avg_out(pop(i).idx);
        update_ylim(H(i).axRate, sub_R, cfg.y_padding);

        % 2. INPUT
        [sub_In, mean_In] = get_plot_data(pop(i).hist.avg_in, iTrial, pop(i).idxPlot);
        update_lines(H(i).lIn, x_vec, sub_In, mean_In);
        H(i).barIn.Data = avg_in(pop(i).idx);
        update_ylim(H(i).axIn, sub_In, cfg.y_padding);

        % 3. THRESHOLD
        [sub_Th, mean_Th] = get_plot_data(pop(i).hist.theta, iTrial, pop(i).idxPlot);
        update_lines(H(i).lTh, x_vec, sub_Th, mean_Th);
        H(i).barTh.Data = theta(pop(i).idx);
        update_ylim(H(i).axTh, sub_Th, cfg.y_padding);

        % 4. ERROR 
        if ~strcmpi(pop(i).learningRule, 'none')
            [sub_Err, mean_Err] = get_plot_data(pop(i).hist.err, iTrial, pop(i).idxPlot);
            update_lines(H(i).lErr, x_vec, sub_Err, mean_Err);
            H(i).barErr.Data = theta(pop(i).idx);
            update_ylim(H(i).axErr, sub_Err, cfg.y_padding);
        end

        drawnow limitrate;
    end
end


    %% ========================================================================
    %  HELPER FUNCTIONS (Same as before, copy-paste them here)
    %  ========================================================================

    function h = init_lines(ax, pStruct, baseColor, mapFunc, alpha)
    % Initialize plot objects for individual neurons + mean
    % Returns struct: h.indiv (array of lines), h.mean (single line)

    nPlot = length(pStruct.idxPlot);
    colors = mapFunc(baseColor, nPlot);

    h.indiv = gobjects(nPlot, 1);
    for i = 1:nPlot
        % Use thin lines for individuals
        h.indiv(i) = plot(ax, NaN, NaN, 'Color', [colors(i,:) alpha], 'LineWidth', 0.5);
    end
    % Use thick, darker line for mean
    meanColor = baseColor * 0.8;
    h.mean = plot(ax, NaN, NaN, 'Color', meanColor, 'LineWidth', 3, 'DisplayName', 'Mean');
    end

    function update_lines(hStruct, x, y_data, y_mean)
    % Vectorized update of line objects
    % 'y_data' should be (Trials x Neurons)

    if ~isempty(hStruct.indiv)
        % Create cell array of columns for fast set()
        y_cells = num2cell(y_data, 1);
        set(hStruct.indiv, 'XData', x, {'YData'}, y_cells');
    end
    set(hStruct.mean, 'XData', x, 'YData', y_mean);
    end

    function [subset_data, mean_data] = get_plot_data(full_hist, iTrial, idx_plot)
    % Helper to slice history and calculate mean
    % full_hist: (nTrial x N) matrix of all neurons
    % idx_plot:  Indices of neurons to extract for individual plotting

    % Extract valid portion up to current trial
    current_slice = full_hist(1:iTrial, :);

    % Calculate Mean across ALL neurons (dim 2)
    mean_data = mean(current_slice, 2, 'omitnan');

    % Extract Subset for individual traces
    if isempty(idx_plot)
        subset_data = current_slice;
    else
        % Ensure indices are valid
        valid_idx = idx_plot(idx_plot <= size(full_hist, 2));
        subset_data = current_slice(:, valid_idx);
    end
    end

    function update_ylim(ax, data, padding)
    % Dynamically updates Y-limits based on data range
    if isempty(data) || all(isnan(data(:))), return; end

    % Calculate bounds
    ymin = min(data(:), [], 'omitnan');
    ymax = max(data(:), [], 'omitnan');

    if isempty(ymin) || isempty(ymax) || isnan(ymin) || isnan(ymax), return; end
    if ymin == ymax, ymax = ymin + 1e-6; end % Prevent flat limits (error)

    % Apply padding
    range = ymax - ymin;
    if range == 0, range = 1; end
    new_lim = [ymin - range * padding, ymax + range * padding];

    set(ax, 'YLim', new_lim);
    end

    function draw_perturbation_bar(ax, popStruct, color)
    % Adds perturbation indicators to a specific axis
    p = popStruct.perturb; % [Start, End, Amp]
    if p(3) ~= 0

        % Gray bar at bottom
        line(ax, [p(1) p(2)], [0 0], 'Color', [0.7 0.7 0.7], 'LineWidth', 5);

        % Text label
        text(ax, mean(p(1:2)), 0, sprintf('+%.1f', p(3)), ...
            'Color', color, 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
    end
    end
