%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL: Input-Centric Homeostasis via Intrinsic Plasticity
% VERSION 2: Hybrid Architecture (Organized + Fast)
%
% ARCHITECTURE:
% - Configuration stored in structures (organized, easy to modify)
% - Direct variables (E, I, thetaE) used in performance-critical loops (fast)
% - Direct matrices (W_EE, W_EI, etc.) for connectivity (vectorized)
% - History structure for organized data storage
%
% DESIGN PHILOSOPHY:
% This hybrid approach maintains the performance of direct variable access
% while providing the organizational benefits of structured data. Since we
% only simulate 2 populations (E and I), we use simple indexing rather
% than complex generalization.
%
% THEORETICAL BACKGROUND
% -------------------------------------------------------------------------
% Steady-state, dynamic stability, and homeostasis
% Mean-field theory provides self-consistency equations that define the fixed points
% of a system, which are states where the network's output is in equilibrium with its input.
% E_ss = f_E ( W_{EE}E_{ss} - W_{EI}I_{ss} )
%
% However, the existence of a fixed point does not guarantee its stability, as even
% minor perturbations can cause the system to diverge. An Inhibition-Stabilized Network (ISN)
% represents a specific parameter regime that provides rapid, dynamic stability.
% In an ISN, the excitatory sub-network is inherently unstable, with recurrent connections
% strong enough to cause runaway activity (W_{EE}g_E-1 > 0).
%
% This explosive excitation is controlled by powerful and fast inhibitory feedback,
% which is recruited by any increase in excitatory firing and acts to quickly restore equilibrium.
% However, if the network is subjected to a sustained perturbation, such as a continuous
% external current, it will simply settle at a new stable fixed point that balances the new input.
% The ISN mechanism has no intrinsic goal or setpoint; it finds the nearest stable equilibrium
% for a given set of parameters and inputs. In contrast, homeostasis is the distinct,
% slower process that actively restores a target firing rate.
%
% Analytic Derivation:
% The model considers a two-population network of interconnected excitatory (E) and inhibitory (I) neurons.
% 1. Fast network dynamics equilibrium:
%    f_E^{-1}(E_{ss}) = W_{EE}E_{ss} - W_{EI}I_{ss}
%
% 2. Slow homeostatic plasticity equilibrium:
%    tau_f * df_E/dt proportional to (h_{set,E} - f_E^{-1}(E))
%    At equilibrium: f_E^{-1}(E_{ss}) = h_{set,E}
%
% Combining these reveals that the network's synaptic input must equal the neuron's internal setpoint:
% W_{EE}E_{ss} - W_{EI}I_{ss} = h_{set,E}
%
% BIOLOGICAL BASIS
% -------------------------------------------------------------------------
% The input-centric model requires a physical mechanism that can measure synaptic input,
% compare it to an internal setpoint, and generate an error signal to drive plasticity.
% The proposed biological implementation is the mitochondrion.
% - Sensor: Mitochondrial Calcium Uniporter (MCU) sensing high-freq firing (proxy for input cost).
% - Setpoint: Basal mitochondrial calcium [Ca2+]_mito.
% - Effector: Adjustment of intrinsic excitability (thresholds) via ion channel expression.
%
% Resolution of Sign Reversal in High-Gain Inhibition Regimes
% -------------------------------------------------------------------------
% The stabilization of the network was achieved by refining the "sensed input" variable
% to track exclusively excitatory drive (Input_sense = W_EE*E + I_ext), rather than
% net synaptic input (W_EE*E - W_EI*I). Previously, the inclusion of the inhibitory term
% caused a sign reversal in the homeostatic control loop. In the Inhibition-Stabilized
% Network (ISN) regime used here — characterized by strong inhibitory feedback gain —
% any increase in excitatory firing recruits a disproportionately larger inhibitory response.
% Consequently, as the network activity increases, the net synaptic input paradoxically
% decreases. This inverts the error signal, causing the homeostatic rule to maladaptively
% lower thresholds in response to perceived low input, driving the network into a runaway
% high-activity state. Redefining the sensed variable to represent metabolic cost
% (modeled as excitatory current driving mitochondrial calcium accumulation), restores
% the necessary monotonic relationship between network activity and the error signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
rng(44); % Ensures consistent random numbers for reproducibility


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONFIGURATION STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All simulation parameters organized in a clear structure

% -------------------------------------------------------------------------
% Simulation Timing
% -------------------------------------------------------------------------
config.dt = 0.0001;              % Integration step (s) - 0.1 ms
config.trial_dur = 2;             % Duration per trial (s)
config.nTrial = 500;               % Total number of trials
config.tau_measure = 1.0;          % Window for homeostatic averaging (s)

% Derived
config.trial_steps = config.trial_dur / config.dt;

% -------------------------------------------------------------------------
% Simulation Flags
% -------------------------------------------------------------------------
config.flg_graphics = true;        % Enable real-time plotting
config.flg_homeostatic = true;     % Enable plasticity
config.graphics_step = 5;          % Update plot every N trials

% -------------------------------------------------------------------------
% Learning Rule
% -------------------------------------------------------------------------
% Options: 'input_centric' - regulates excitatory synaptic drive
%          'output_centric' - regulates firing rate
config.learning_rule = 'input_centric';

% -------------------------------------------------------------------------
% Population Sizes
% -------------------------------------------------------------------------
config.nE = 80;                    % Excitatory population size
config.nI = 20;                    % Inhibitory population size

% -------------------------------------------------------------------------
% Excitatory Population Parameters
% -------------------------------------------------------------------------
config.E.tau = 0.010;              % Time constant (s) - 10 ms
config.E.gain = 1;                 % Activation function gain
config.E.theta_init = 4.8;         % Initial firing threshold
config.E.max_rate = 100;           % Maximum firing rate (Hz)
config.E.learningRule = 'input_centric';  % Plasticity rule
config.E.alpha = 0.005;            % Learning rate

% Setpoint depends on learning rule (set below after rule definition)

% -------------------------------------------------------------------------
% Inhibitory Population Parameters
% -------------------------------------------------------------------------
config.I.tau = 0.002;              % Time constant (s) - 2 ms (fast for ISN)
config.I.gain = 4;                 % Activation function gain
config.I.theta_init = 15;          % Firing threshold (scalar, not per-neuron)
config.I.max_rate = 250;           % Maximum firing rate (Hz)
config.I.learningRule = 'none';    % No plasticity for I population

% -------------------------------------------------------------------------
% Connectivity Parameters (ISN Regime)
% -------------------------------------------------------------------------
% Scaling factors adapted from Soldado et al. (PNAS, 2022)
config.connectivity.J_EE = 5.00;   % E -> E strength
config.connectivity.J_EI = 1.52;   % I -> E strength
config.connectivity.J_IE = 10.0;   % E -> I strength
config.connectivity.J_II = 2.25;   % I -> I strength
config.connectivity.W_mean = 0.1;  % Mean weight
config.connectivity.W_sigma = 0.04; % Weight heterogeneity

% -------------------------------------------------------------------------
% Perturbation Settings
% -------------------------------------------------------------------------
% Format: [OnsetTrial, OffsetTrial, Amplitude]
config.perturb.E = [150, 250, 5.0];  % External input to E population
config.perturb.I = [150, 250, 5.0];  % External input to I population

% -------------------------------------------------------------------------
% Noise Parameters (Ornstein-Uhlenbeck Process)
% -------------------------------------------------------------------------
config.noise.tau = 0.1;            % Correlation time (s)
config.noise.mu = 0;                % Mean
config.noise.sigma = 4.0;           % Amplitude

% -------------------------------------------------------------------------
% Homeostatic Setpoint Calculation
% -------------------------------------------------------------------------
switch config.learning_rule
    case 'input_centric'
        % GOAL: Maintain excitatory synaptic drive at specific level
        % Setpoint ensures E-cells provide sufficient drive to keep I-cells active
        SafetyMargin = 5.0;
        config.E.target = config.I.theta_init + SafetyMargin;  % ~20 a.u.
        
    case 'output_centric'
        % GOAL: Maintain firing rate at specific level
        config.E.target = 5.0;  % Hz
        
    otherwise
        error('Unknown learning rule: %s', config.learning_rule);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Activation Function (ReLU)
% -------------------------------------------------------------------------
F = @(x, gain, Thr) gain * max(0, x - Thr);

% -------------------------------------------------------------------------
% State Variables (Direct variables for performance)
% -------------------------------------------------------------------------
E = rand(config.nE, 1) * 5;           % Excitatory firing rates
I = rand(config.nI, 1) * 14;           % Inhibitory firing rates
thetaE = ones(config.nE, 1) * config.E.theta_init;  % E thresholds (per-neuron)
thetaI = config.I.theta_init;          % I threshold (scalar)

% Noise state
OUE = zeros(config.nE, 1);
OUI = zeros(config.nI, 1);

% -------------------------------------------------------------------------
% Connectivity Matrices (Direct matrices for vectorized operations)
% -------------------------------------------------------------------------
% Generate random connectivity with ISN scaling
J_EE = config.connectivity.J_EE / sqrt(config.nE);
J_EI = config.connectivity.J_EI / sqrt(config.nI);
J_IE = config.connectivity.J_IE / sqrt(config.nE);
J_II = config.connectivity.J_II / sqrt(config.nI);

W_EE = J_EE * (config.connectivity.W_mean + ...
                config.connectivity.W_sigma * randn(config.nE, config.nE));
W_EI = J_EI * (config.connectivity.W_mean + ...
                config.connectivity.W_sigma * randn(config.nE, config.nI));
W_IE = J_IE * (config.connectivity.W_mean + ...
                config.connectivity.W_sigma * randn(config.nI, config.nE));
W_II = J_II * (config.connectivity.W_mean + ...
                config.connectivity.W_sigma * randn(config.nI, config.nI));

% Remove self-connections
W_EE(1:config.nE+1:end) = 0;
W_II(1:config.nI+1:end) = 0;

% -------------------------------------------------------------------------
% History Structure (Organized data storage)
% -------------------------------------------------------------------------
history.meanE = NaN(config.nTrial, 1);
history.meanI = NaN(config.nTrial, 1);
history.meanInput = NaN(config.nTrial, 1);
history.meanTheta = NaN(config.nTrial, 1);
history.perturbE = zeros(config.nTrial, 1);
history.perturbI = zeros(config.nTrial, 1);
history.regvar = NaN(config.nTrial, 1);  % Regulated variable (input or rate)

% Store full time series for selected neurons (subset to save memory)
nPlotE = min(config.nE, 10);
nPlotI = min(config.nI, 5);
idxPlotE = round(linspace(1, config.nE, nPlotE));
idxPlotI = round(linspace(1, config.nI, nPlotI));

history.allE = NaN(config.nTrial, nPlotE);
history.allI = NaN(config.nTrial, nPlotI);
history.allTheta = NaN(config.nTrial, nPlotE);
history.allInput = NaN(config.nTrial, nPlotE);

% -------------------------------------------------------------------------
% Pre-allocate trial storage (for averaging)
% -------------------------------------------------------------------------
trial_steps = config.trial_steps;
tempE = zeros(config.nE, trial_steps);
tempI = zeros(config.nI, trial_steps);
tempInput = zeros(config.nE, trial_steps);  % Excitatory drive only


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize graphics
if config.flg_graphics
    InputCentric_Demo_v2_graphics;
end

for trial = 1:config.nTrial
    
    % ---------------------------------------------------------------------
    % External Input (Perturbation)
    % ---------------------------------------------------------------------
    ExtInputE = 0;
    ExtInputI = 0;
    
    if trial >= config.perturb.E(1) && trial <= config.perturb.E(2)
        ExtInputE = config.perturb.E(3);
    end
    
    if trial >= config.perturb.I(1) && trial <= config.perturb.I(2)
        ExtInputI = config.perturb.I(3);
    end
    
    history.perturbE(trial) = ExtInputE;
    history.perturbI(trial) = ExtInputI;
    
    % ---------------------------------------------------------------------
    % FAST DYNAMICS: Within-Trial Recurrent Network
    % ---------------------------------------------------------------------
    % Vectorized operations using direct matrices for maximum performance
    
    for t = 1:trial_steps
        
        % Update noise (Ornstein-Uhlenbeck process)
        OUE = OUE + (config.dt/config.noise.tau) * (config.noise.mu - OUE) + ...
                    sqrt(2*config.dt/config.noise.tau) * config.noise.sigma * randn(config.nE, 1);
        OUI = OUI + (config.dt/config.noise.tau) * (config.noise.mu - OUI) + ...
                    sqrt(2*config.dt/config.noise.tau) * config.noise.sigma * randn(config.nI, 1);
        
        % Calculate synaptic drives (vectorized matrix operations)
        TotalDriveE = W_EE*E - W_EI*I + ExtInputE + OUE;
        TotalDriveI = W_IE*E - W_II*I + ExtInputI + OUI;
        
        % Update firing rates (Euler integration)
        dEdt = (-E + F(TotalDriveE, config.E.gain, thetaE)) / config.E.tau;
        dIdt = (-I + F(TotalDriveI, config.I.gain, thetaI)) / config.I.tau;
        
        E = E + dEdt * config.dt;
        I = I + dIdt * config.dt;
        
        % Hard saturation (prevents numerical explosion)
        E(E > config.E.max_rate) = config.E.max_rate;
        I(I > config.I.max_rate) = config.I.max_rate;
        
        % Store state for averaging
        tempE(:, t) = E;
        tempI(:, t) = I;
        
        % Store excitatory drive (for input-centric plasticity)
        % Note: This is W_EE*E + Ext + Noise (excitatory drive only, not net)
        tempInput(:, t) = W_EE*E + ExtInputE + OUE;
    end
    
    % ---------------------------------------------------------------------
    % SLOW DYNAMICS: Homeostatic Plasticity
    % ---------------------------------------------------------------------
    
    % Calculate trial averages over measurement window
    trial_ss = floor(config.tau_measure / config.dt);
    range = (trial_steps - trial_ss + 1):trial_steps;
    avgE = mean(tempE(:, range), 2);
    avgI = mean(tempI(:, range), 2);
    avgInput = mean(tempInput(:, range), 2);
    
    % Update history
    history.meanE(trial) = mean(avgE);
    history.meanI(trial) = mean(avgI);
    history.meanInput(trial) = mean(avgInput);
    history.meanTheta(trial) = mean(thetaE);
    history.allE(trial, :) = avgE(idxPlotE)';
    history.allI(trial, :) = avgI(idxPlotI)';
    history.allTheta(trial, :) = thetaE(idxPlotE)';
    history.allInput(trial, :) = avgInput(idxPlotE)';
    
    % Apply plasticity rule
    if config.flg_homeostatic && ~strcmpi(config.E.learningRule, 'none')
        
        % Determine what variable is being regulated
        switch config.learning_rule
            case 'input_centric'
                MeasuredVal = avgInput;  % Regulate excitatory drive
            case 'output_centric'
                MeasuredVal = avgE;      % Regulate firing rate
        end
        
        % Calculate error signal
        % Positive error = measured value too high -> increase threshold
        err = MeasuredVal - config.E.target;
        
        % Update thresholds
        thetaE = thetaE + config.E.alpha * err;
        
        % Store regulated variable for analysis
        history.regvar(trial) = mean(MeasuredVal);
    end
    
    % ---------------------------------------------------------------------
    % Visualization
    % ---------------------------------------------------------------------
    if config.flg_graphics && (mod(trial, config.graphics_step) == 0 || trial == 1)
        InputCentric_Demo_v2_graphics;
    end
    
end

