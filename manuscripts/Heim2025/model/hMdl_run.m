function hMdl_run()
    % Main function to run and compare the two homeostatic models.
    % Author: Your Name/AI Assistant
    % Date: October 26, 2023

    close all; clear; clc;

    rng(42); % For reproducibility

    fprintf('Starting Homeostasis Model Comparison...\n');

    % --- Simulation Parameters ---
    params.Ne = 800; % Number of excitatory neurons
    params.Ni = 200; % Number of inhibitory neurons
    params.N = params.Ne + params.Ni;

    % Time parameters
    params.T = 10 * 1000; % Example: 10 seconds (change as needed)
    params.dt = 1;
    params.sim_steps = params.T / params.dt;

    % Perturbation parameters (as fractions of T)
    params.perturb_fraction = 0.9; % Fraction of excitatory neurons to perturb
    params.perturb_start_time = round(0.3 * params.T); % 30% into the simulation
    params.perturb_end_time   = round(0.8 * params.T); % 80% into the simulation
    params.perturb_strength = -20; % Stronger perturbation

    % Homeostasis parameters
    params.eta = 2e-5; % Homeostasis learning rate (slow!)
    params.f_target = 5; % Target firing rate in Hz for Rule A
    params.homeo_update_interval = round(0.05 * params.T); % e.g., every 5% of T

    % --- Run Simulation for Rule A: Output-Centric ---
    fprintf('Running Simulation for Rule A (Output-Centric)...\n');
    params.rule = 'output';
    results_A = network_simulation(params);

    % --- Run Simulation for Rule B: Input-Centric ---
    fprintf('Running Simulation for Rule B (Input-Centric)...\n');
    params.rule = 'input';
    results_B = network_simulation(params);

    % --- Plotting Results ---
    fprintf('Plotting results...\n');
    plot_comparison(results_A, results_B, params);
    fprintf('Done.\n');
end

function results = network_simulation(params)
    % Core function to run the network simulation for a given rule.

    % --- Unpack Parameters ---
    Ne = params.Ne; Ni = params.Ni; N = params.N;
    dt = params.dt; sim_steps = params.sim_steps;

    % --- Neuron Model Parameters (Izhikevich) ---
    re = rand(Ne, 1);
    ri = rand(Ni, 1);
    a = [0.02 * ones(Ne, 1); 0.02 + 0.08 * ri];
    b = [0.2 * ones(Ne, 1);  0.25 - 0.05 * ri]; % 'b' is the homeostatic variable for E-neurons
    c = [-65 + 15 * re.^2;   -65 * ones(Ni, 1)];
    d = [8 - 6 * re.^2;      2 * ones(Ni, 1)];

    b_initial = b; % Store initial 'b' values

    % --- Synaptic Connections ---
    S = [0.5 * rand(Ne + Ni, Ne), -1 * rand(Ne + Ni, Ni)];
    % No self-connections
    for i = 1:N
        S(i, i) = 0;
    end
    % Sparsify connections (e.g., 10% probability)
    S(rand(N, N) > 0.1) = 0;

    % --- Initialization ---
    v = -65 * ones(N, 1); % Initial membrane potential
    u = b .* v;           % Initial recovery variable
    firings = [];         % Stores [time_step, neuron_id]

    % History tracking for plotting
    num_intervals = floor(sim_steps / params.homeo_update_interval);
    FR_hist = zeros(N, num_intervals);
    I_hist = zeros(Ne, num_intervals); % Only track for E-neurons
    b_hist = zeros(Ne, num_intervals);

    % --- Determine Perturbed Population ---
    e_indices = 1:Ne;
    num_perturbed = floor(params.perturb_fraction * Ne);
    perturbed_indices = e_indices(1:num_perturbed);
    unperturbed_indices = e_indices(num_perturbed+1:end);

    % --- Establish Target Input Current for Rule B ---
    % Run a short pre-simulation to find the baseline average input current
    fprintf('  Establishing baseline for target variables...\n');
    I_baseline_vals = [];
    v_pre = -65 * ones(N, 1); u_pre = b .* v_pre;
    for t = 1:5000 % 5 second pre-run
        I_bg = [5 * randn(Ne, 1); 2 * randn(Ni, 1)]; % Background drive
        fired = find(v_pre >= 30);
        I_syn = S * (v_pre >= 30);
        I_total = I_bg + I_syn;
        
        if mod(t, 10) == 0 % Sample every 10ms
            I_baseline_vals = [I_baseline_vals, I_total(1:Ne)];
        end

        v_pre(fired) = c(fired);
        u_pre(fired) = u_pre(fired) + d(fired);
        v_pre = v_pre + dt * (0.04*v_pre.^2 + 5*v_pre + 140 - u_pre + I_total);
        u_pre = u_pre + dt * (a .* (b .* v_pre - u_pre));
    end
    I_target = mean(I_baseline_vals, 'all');
    params.I_target = I_target; % Store for plotting
    fprintf('  Target Input Current (for Rule B) set to: %.2f\n', I_target);


    % --- Main Simulation Loop ---
    interval_spikes = [];
    interval_I_sum = zeros(Ne, 1);
    interval_step_count = 0;

    for t = 1:sim_steps
        % Background and Perturbation Currents
        I_bg = [5 * randn(Ne, 1); 2 * randn(Ni, 1)]; % Background drive
        I_perturb = zeros(N, 1);
        if t >= params.perturb_start_time && t < params.perturb_end_time
            I_perturb(perturbed_indices) = params.perturb_strength;
        end

        % Synaptic Input
        fired = find(v >= 30);
        if ~isempty(fired)
            firings = [firings; t+0*fired, fired];
            interval_spikes = [interval_spikes; t+0*fired, fired];
        end
        I_syn = S * (v >= 30);
        I_total = I_bg + I_syn + I_perturb;
        
        % Accumulate input current for E-neurons for this interval
        interval_I_sum = interval_I_sum + I_total(1:Ne);
        interval_step_count = interval_step_count + 1;

        % Neuron state update
        v(fired) = c(fired);
        u(fired) = u(fired) + d(fired);
        v = v + dt * (0.04*v.^2 + 5*v + 140 - u + I_total);
        u = u + dt * (a .* (b .* v - u));

        % --- Homeostasis Update ---
        if mod(t, params.homeo_update_interval) == 0
            interval_idx = t / params.homeo_update_interval;
            
            % Calculate current firing rate for all neurons (in Hz)
            f_current = zeros(N, 1);
            for n_idx = 1:N
                n_spikes = sum(interval_spikes(:, 2) == n_idx);
                f_current(n_idx) = n_spikes / (params.homeo_update_interval / 1000);
            end
            FR_hist(:, interval_idx) = f_current;

            % Calculate average input current for E-neurons
            I_current = interval_I_sum / interval_step_count;
            I_hist(:, interval_idx) = I_current;
            
            % Apply the homeostasis rule ONLY to Excitatory neurons
            if strcmp(params.rule, 'output')
                error = FR_hist(1:Ne, interval_idx) - params.f_target;
            elseif strcmp(params.rule, 'input')
                error = I_hist(:, interval_idx) - params.I_target;
            end
            
            % Update homeostatic variable 'b'
            delta_b = params.eta * error;
            b(1:Ne) = b(1:Ne) + delta_b;
            % Clamp 'b' to reasonable values to prevent instability
            b(1:Ne) = max(0.1, min(0.3, b(1:Ne)));

            b_hist(:, interval_idx) = b(1:Ne);
            
            % Reset interval counters
            interval_spikes = [];
            interval_I_sum = zeros(Ne, 1);
            interval_step_count = 0;
        end
    end
    
    % Package results
    results.firings = firings;
    results.FR_hist = FR_hist;
    results.I_hist = I_hist;
    results.b_hist = b_hist;
    results.params = params;
    results.perturbed_indices = perturbed_indices;
    results.unperturbed_indices = unperturbed_indices;
    results.i_indices = (Ne+1):N;
end

function plot_comparison(resA, resB, params)
    % Function to generate comparison plots for the two simulations.
    
    figure('Position', [100, 100, 1600, 1000], 'Color', 'w');
    
    % --- Helper function for smoothing ---
    smooth_data = @(data, sigma) imgaussfilt(data, 1, 'FilterDomain', 'frequency'); % or use sigma=0 for no smoothing
    
    % --- Column 1: Output-Centric Results (Rule A) ---
    ax1 = subplot(4, 2, 1);
    plot(resA.firings(:,1)/1000, resA.firings(:,2), 'k.', 'MarkerSize', 1);
    hold on;
    % Highlight perturbation period
    patch([params.perturb_start_time, params.perturb_end_time, params.perturb_end_time, params.perturb_start_time]/1000, [0, 0, params.N, params.N], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    title('Rule A: Output-Centric (Raster Plot)');
    xlabel('Time (s)'); ylabel('Neuron ID');
    xlim([0, params.T/1000]); ylim([0, params.N]);
    
    ax2 = subplot(4, 2, 3);
    time_axis = (1:size(resA.FR_hist, 2)) * params.homeo_update_interval / 1000;
    plot(time_axis, smooth_data(mean(resA.FR_hist(resA.perturbed_indices,:)), 3), 'r', 'LineWidth', 2);
    hold on;
    plot(time_axis, smooth_data(mean(resA.FR_hist(resA.unperturbed_indices,:)), 3), 'b', 'LineWidth', 2);
    plot(time_axis, smooth_data(mean(resA.FR_hist(resA.i_indices,:)), 3), 'g', 'LineWidth', 2);
    yline(params.f_target, 'k--', 'Label', 'Target Firing Rate');
    patch([params.perturb_start_time, params.perturb_end_time, params.perturb_end_time, params.perturb_start_time]/1000, [0, 0, 30, 30], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    title('Firing Rates (Smoothed)');
    xlabel('Time (s)'); ylabel('Avg Firing Rate (Hz)');
    legend('Perturbed E', 'Unperturbed E', 'Inhibitory');
    xlim([0, params.T/1000]); ylim([0, 30]);

    ax3 = subplot(4, 2, 5);
    plot(time_axis, smooth_data(mean(resA.I_hist(resA.perturbed_indices - min(resA.perturbed_indices) + 1,:)), 3), 'r', 'LineWidth', 2);
    hold on;
    plot(time_axis, smooth_data(mean(resA.I_hist(resA.unperturbed_indices - min(resA.unperturbed_indices) + 1,:)), 3), 'b', 'LineWidth', 2);
    patch([params.perturb_start_time, params.perturb_end_time, params.perturb_end_time, params.perturb_start_time]/1000, [-10, -10, 20, 20], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    title('Synaptic Input Current (E-Neurons)');
    xlabel('Time (s)'); ylabel('Avg Input Current');
    legend('Perturbed E', 'Unperturbed E');
    xlim([0, params.T/1000]); ylim([-10, 20]);
    
    ax4 = subplot(4, 2, 7);
    plot(time_axis, mean(resA.b_hist(resA.perturbed_indices - min(resA.perturbed_indices) + 1,:)), 'r', 'LineWidth', 2);
    hold on;
    plot(time_axis, mean(resA.b_hist(resA.unperturbed_indices - min(resA.unperturbed_indices) + 1,:)), 'b', 'LineWidth', 2);
    patch([params.perturb_start_time, params.perturb_end_time, params.perturb_end_time, params.perturb_start_time]/1000, [0.1, 0.1, 0.3, 0.3], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    title('Homeostatic Variable `b` (E-Neurons)');
    xlabel('Time (s)'); ylabel('Avg `b` value');
    legend('Perturbed E', 'Unperturbed E', 'Location', 'south');
    xlim([0, params.T/1000]); ylim([0.1, 0.3]);
    
    % --- Column 2: Input-Centric Results (Rule B) ---
    ax5 = subplot(4, 2, 2);
    plot(resB.firings(:,1)/1000, resB.firings(:,2), 'k.', 'MarkerSize', 1);
    hold on;
    patch([params.perturb_start_time, params.perturb_end_time, params.perturb_end_time, params.perturb_start_time]/1000, [0, 0, params.N, params.N], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    title('Rule B: Input-Centric (Raster Plot)');
    xlabel('Time (s)'); ylabel('Neuron ID');
    xlim([0, params.T/1000]); ylim([0, params.N]);

    ax6 = subplot(4, 2, 4);
    time_axis = (1:size(resB.FR_hist, 2)) * params.homeo_update_interval / 1000;
    plot(time_axis, smooth_data(mean(resB.FR_hist(resB.perturbed_indices,:)), 3), 'r', 'LineWidth', 2);
    hold on;
    plot(time_axis, smooth_data(mean(resB.FR_hist(resB.unperturbed_indices,:)), 3), 'b', 'LineWidth', 2);
    plot(time_axis, smooth_data(mean(resB.FR_hist(resB.i_indices,:)), 3), 'g', 'LineWidth', 2);
    patch([params.perturb_start_time, params.perturb_end_time, params.perturb_end_time, params.perturb_start_time]/1000, [0, 0, 30, 30], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    title('Firing Rates (Smoothed)');
    xlabel('Time (s)'); ylabel('Avg Firing Rate (Hz)');
    legend('Perturbed E', 'Unperturbed E', 'Inhibitory');
    xlim([0, params.T/1000]); ylim([0, 30]);

    ax7 = subplot(4, 2, 6);
    plot(time_axis, smooth_data(mean(resB.I_hist(resB.perturbed_indices - min(resB.perturbed_indices) + 1,:)), 3), 'r', 'LineWidth', 2);
    hold on;
    plot(time_axis, smooth_data(mean(resB.I_hist(resB.unperturbed_indices - min(resB.unperturbed_indices) + 1,:)), 3), 'b', 'LineWidth', 2);
    yline(resB.params.I_target, 'k--', 'Label', 'Target Input Current');
    patch([params.perturb_start_time, params.perturb_end_time, params.perturb_end_time, params.perturb_start_time]/1000, [-10, -10, 20, 20], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    title('Synaptic Input Current (E-Neurons)');
    xlabel('Time (s)'); ylabel('Avg Input Current');
    legend('Perturbed E', 'Unperturbed E');
    xlim([0, params.T/1000]); ylim([-10, 20]);
    
    ax8 = subplot(4, 2, 8);
    plot(time_axis, mean(resB.b_hist(resB.perturbed_indices - min(resB.perturbed_indices) + 1,:)), 'r', 'LineWidth', 2);
    hold on;
    plot(time_axis, mean(resB.b_hist(resB.unperturbed_indices - min(resB.unperturbed_indices) + 1,:)), 'b', 'LineWidth', 2);
    patch([params.perturb_start_time, params.perturb_end_time, params.perturb_end_time, params.perturb_start_time]/1000, [0.1, 0.1, 0.3, 0.3], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    title('Homeostatic Variable `b` (E-Neurons)');
    xlabel('Time (s)'); ylabel('Avg `b` value');
    legend('Perturbed E', 'Unperturbed E', 'Location', 'south');
    xlim([0, params.T/1000]); ylim([0.1, 0.3]);
    
    sgtitle('Comparison of Homeostatic Rules', 'FontSize', 16, 'FontWeight', 'bold');
end