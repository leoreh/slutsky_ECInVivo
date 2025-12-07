%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL: Input-Centric Homeostasis via Intrinsic Plasticity
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolution of Sign Reversal in High-Gain Inhibition Regimes
% The stabilization of the network was achieved by refining the "sensed
% input" variable to track exclusively excitatory drive ($Input_{sense} =
% W_{EE}E + I_{ext}$), rather than net synaptic input ($W_{EE}E -
% W_{EI}I$). Previously, the inclusion of the inhibitory term caused a sign
% reversal in the homeostatic control loop. In the Inhibition-Stabilized
% Network (ISN) regime used here — characterized by strong inhibitory
% feedback gain — any increase in excitatory firing recruits a
% disproportionately larger inhibitory response. Consequently, as the
% network activity increases, the net synaptic input paradoxically
% decreases. This inverts the error signal, causing the homeostatic rule to
% maladaptively lower thresholds in response to perceived low input,
% driving the network into a runaway high-activity state. Redefining the
% sensed variable to represent metabolic cost (modeled as excitatory
% current driving mitochondrial calcium accumulation), restores the
% necessary monotonic relationship between network activity and the error
% signal. This modification is biophysically justified, as mitochondrial
% calcium uptake is primarily driven by excitatory (e.g., NMDA-mediated or
% voltage-gated) calcium influx rather than inhibitory chloride currents,
% and it successfully allows the homeostatic mechanism to locate and
% maintain a stable operating point within the ISN regime.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
% close all
rng(44) % Ensures consistent random numbers for reproducibility


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TIMING PARAMETERS
dt = 0.0001;                % 0.1 ms integration step
trial_dur = 2;              % Duration of one trial (seconds)
trial_steps = trial_dur/dt; % Number of simulation steps per trial
nTrial = 500;               % Total number of trials

% FLAGS
flg_graphics = 1;           % Real-time plotting
flg_homeostatic = 1;        % Enable plasticity
graphics_step = 5;          % Update plot every N trials (speed optimization)

% LEARNING RULE 
learning_rule = 'input_centric'; % Options: 'input_centric', 'output_centric'


% PERTURBATION
% Defined as [Onset, Offset] trials for E and I separately
PerturbRangeE = [150, 250];
PerturbRangeI = [150, 250];

PerturbAmpE = 5.0;       % External Input to E population (a.u.)
PerturbAmpI = 5.0;       % External Input to I population (a.u.)

% NOISE PARAMETERS
% Ornstein-Uhlenbeck process
OUtau = 0.1;            % Correlation time
OUmu = 0;               % Mean noise
OUsigma = 4.0;          % Noise amplitude (scaled up for visibility)

% NEURON PARAMTERS
% Activation function (ReLU)
F = @(x,gain,Thr) gain*max(0,x-Thr);

nE = 80 * 1;    % Excitatory population size
nI = 20 * 1;    % Inhibitory population size

% Intrinsic Properties
theta_init = 4.8;       % firing threshold
thetaE = ones(nE,1) * theta_init;
thetaI = 15;
gainE = 1;
gainI = 4;

% Time constants
tauE = 0.010;  % 10 ms
tauI = 0.002;  % 2 ms (Fast inhibition is critical for ISN)

% Saturation prevents numerical explosion
maxfE = 100;    % Hz
maxfI = 250;    % Hz

% Initial State Variables
E = rand(nE,1) * 5;
I = rand(nI,1) * 14;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONNECTIVITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_sigma = 0.04;  % Heterogeneity
W_mean = 0.1;    % Mean weight

% Connectivity (ISN Regime)
J_EE = 5.00  / sqrt(nE);
J_EI = 1.52  / sqrt(nI);
J_IE = 10.0  / sqrt(nE);
J_II = 2.25  / sqrt(nI);

W_EE = J_EE * (W_mean + W_sigma * randn(nE, nE));
W_EI = J_EI * (W_mean + W_sigma * randn(nE, nI));
W_IE = J_IE * (W_mean + W_sigma * randn(nI, nE));
W_II = J_II * (W_mean + W_sigma * randn(nI, nI));

% Remove self-connections
W_EE(1:nE+1:end) = 0;
W_II(1:nI+1:end) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOMEOSTATIC PLASTICITY RULE DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grouping all rule-dependent logic here for clarity.

% Generic Setpoint Vector (will be populated below)
HomeoSetpoint = zeros(nE, 1);
alpha = 0;

switch learning_rule
    case 'input_centric'
        % GOAL: Maintain total synaptic input at specific level
        % Setpoint logic: Ensure E-cells target a drive that keeps I-cells alive
        SafetyMargin = 5.0;
        InputSetVal = thetaI + SafetyMargin;

        HomeoSetpoint = ones(nE,1) * InputSetVal;

        % Learning Rate
        alpha = 0.005;

    case 'output_centric'
        % GOAL: Maintain Firing Rate (Output) at specific level
        TargetRate = 5.0; % Hz

        HomeoSetpoint = ones(nE,1) * TargetRate;

        % Learning Rate 
        alpha = 0.005;

    otherwise
        error('Unknown learning rule');
end

% Measurement Window
tau_measure = 1.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% History arrays
history.meanE = NaN(nTrial, 1);
history.meanI = NaN(nTrial, 1);
history.meanInput = NaN(nTrial, 1);
history.meanTheta = NaN(nTrial, 1);
history.perturbE = zeros(nTrial, 1);
history.perturbI = zeros(nTrial, 1);

history.allE = NaN(nTrial, nE);
history.allI = NaN(nTrial, nI);
history.allTheta = NaN(nTrial, nE);
history.allInput = NaN(nTrial, nE);

history.regvar = NaN(nTrial, 1); % Store Regulated Variable (Input or Rate)


% Pre-allocation of noise vectors
OUE = zeros(nE,1);
OUI = zeros(nI,1);

% Select fixed subset of neurons for plotting
idxPlotE = round(linspace(1, nE, min(nE, 10)));
idxPlotI = round(linspace(1, nI, min(nI, 5)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATION LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flg_graphics
    InputCentric_Demo_graphics % Initialize Plot
end

for trial = 1 : nTrial

    % EXTERNAL INPUT (Perturbation)
    ExtInputE = 0;
    ExtInputI = 0;

    if trial >= PerturbRangeE(1) && trial <= PerturbRangeE(2)
        ExtInputE = PerturbAmpE;
    end

    if trial >= PerturbRangeI(1) && trial <= PerturbRangeI(2)
        ExtInputI = PerturbAmpI;
    end

    % Store effective perturbation (e.g. just E or both)
    history.perturbE(trial) = ExtInputE;
    history.perturbI(trial) = ExtInputI;

    % FAST DYNAMICS (recurrent E/I)
    tempE = zeros(nE, trial_steps);
    tempI = zeros(nI, trial_steps);
    tempInput = zeros(nE, trial_steps);

    for t = 1 : trial_steps

        % Noise
        OUE = OUE + (dt/OUtau) * (OUmu-OUE) + sqrt(2*dt/OUtau) * OUsigma * randn(nE,1);
        OUI = OUI + (dt/OUtau) * (OUmu-OUI) + sqrt(2*dt/OUtau) * OUsigma * randn(nI,1);

        % Synaptic Input
        TotalDriveE = W_EE*E - W_EI*I + ExtInputE + OUE;
        TotalDriveI = W_IE*E - W_II*I + ExtInputI + OUI;

        % Firing Rates
        dEdt = (-E + F(TotalDriveE, gainE, thetaE)) / tauE;
        dIdt = (-I + F(TotalDriveI, gainI, thetaI)) / tauI;

        E = E + dEdt*dt;
        I = I + dIdt*dt;

        % Saturation
        E(E>maxfE) = maxfE;
        I(I>maxfI) = maxfI;

        % Store state
        tempE(:,t) = E;
        tempI(:,t) = I;

        % Store regulated variable (excitatory drive)
        % Even if output-centric, we track Input for analysis
        tempInput(:,t) = W_EE*E + ExtInputE + OUE;
    end

    % SLOW DYNAMICS (Homeostatic Plasticity)

    % Measure Average Activity
    trial_ss = floor(tau_measure/dt);
    range = (trial_steps - trial_ss + 1) : trial_steps;
    avgE = mean(tempE(:,range), 2);
    avgI = mean(tempI(:,range), 2);
    avgInput = mean(tempInput(:,range), 2);

    % Update History
    history.meanE(trial) = mean(avgE);
    history.meanI(trial) = mean(avgI);
    history.meanInput(trial) = mean(avgInput);
    history.meanTheta(trial) = mean(thetaE);
    history.allE(trial, :) = avgE';
    history.allI(trial, :) = avgI';
    history.allTheta(trial, :) = thetaE';
    history.allInput(trial, :) = avgInput';

    % PLASTICITY RULE
    if flg_homeostatic

        % Measure "Sensed" Variable
        switch learning_rule
            case 'input_centric'
                MeasuredVal = avgInput;
            case 'output_centric'
                MeasuredVal = avgE;
        end

        % Calculate Error
        % (Val - Setpoint): Positive err -> Valence is "Too High" -> Increase Theta
        err = MeasuredVal - HomeoSetpoint;

        % Update Thresholds
        thetaE = thetaE + alpha * err;

        % Store for debug/graphics
        history.regvar(trial) = mean(MeasuredVal);
    end


    % VISUALIZATION
    if flg_graphics && (mod(trial, graphics_step) == 0 || trial == 1)
        InputCentric_Demo_graphics;
    end

end
