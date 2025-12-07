%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

%% ========================================================================
%  CONFIGURATION
%  ========================================================================

rng(44);            % Ensures reproducibile random numbers

% TIMING
dt = 0.0001;        % Integration step (s)
trial_dur = 2;      % Duration per trial (s)
sense_dur = 1.0;    % Window for averaging sensed variable each trial (s)
nTrial = 700;

% FLAGS
flg_graphics = true;
flg_setpointHetero = true;
graphics_step = 10;

% NOISE (Ornstein-Uhlenbeck)
OU_tau = 0.1;
OU_mu = 0;
OU_sigma = 4.0;

% ACTIVATION FUNCTION
% ReLU: Gain * max(0, x - theta)
% Func_Act = @(x, g, th) g .* max(0, x - th);

%% ========================================================================
%  POPULATIONS
%  ========================================================================

% EXCITATORY
pop(1).name = 'Exc';
pop(1).synType = 1;              % +1 = Excitatory
pop(1).N = 80;
pop(1).tau = 0.010;
pop(1).gain = 1;
pop(1).theta_init = 4.8;
pop(1).r_init = 5;
pop(1).max_rate = 100;

% Plasticity
pop(1).learningRule = 'output_centric';

% Perturbation
pop(1).perturb = [700, 700, 10]; % [Start, End, Amp]

% INHIBITORY
pop(2).name = 'Inh';
pop(2).synType = -1;             % -1 = Inhibitory
pop(2).N = 20;
pop(2).tau = 0.002;              % Fast inhibition (Critical for ISN)
pop(2).gain = 4;
pop(2).theta_init = 15;
pop(2).r_init = 14;
pop(2).max_rate = 250;

% Plasticity
pop(2).learningRule = 'none';

% Perturbation
pop(2).perturb = [0, 0, 5.0]; % [Start, End, Amp]


%% ========================================================================
%  COMPILATION
%  ========================================================================

nPop = length(pop);
nUnits = sum([pop.N]);

% Pre-allocate State Vectors
r       = zeros(nUnits, 1);       % Firing Rates
theta   = zeros(nUnits, 1);       % Thresholds
tau     = zeros(nUnits, 1);       % Time Constants
gain    = zeros(nUnits, 1);       % Gains
noise   = zeros(nUnits, 1);       % Noise State
ext     = zeros(nUnits, 1);       % External Input Buffer

% Populate State Vectors
cursor = 1;
for iPop = 1:nPop

    % Define Index Map (Mask)
    pop(iPop).idx = cursor : (cursor + pop(iPop).N - 1);

    % Fill Parameters
    r(pop(iPop).idx)            = rand(pop(iPop).N, 1) * pop(iPop).r_init;
    theta(pop(iPop).idx)        = pop(iPop).theta_init;
    tau(pop(iPop).idx)          = pop(iPop).tau;
    gain(pop(iPop).idx)         = pop(iPop).gain;

    % Initialize History (For Graphics)
    pop(iPop).hist.r        = NaN(nTrial, pop(iPop).N);
    pop(iPop).hist.theta    = NaN(nTrial, pop(iPop).N);
    pop(iPop).hist.avg_in   = NaN(nTrial, pop(iPop).N);

    % Select units for plotting
    pop(iPop).idxPlot = round(linspace(1, pop(iPop).N, min(pop(iPop).N, 10)));

    cursor = cursor + pop(iPop).N;
end

% Determine learning rule parameters
for iPop = 1:nPop
    
    pop(iPop).target = nan;
    switch pop(iPop).learningRule

        case 'input_centric'
            
            % homogenous setpoint
            pop(iPop).target = pop(2).theta_init + 5;

            % Lognormal targets for input drive (Median ~ Theta_Inh + 5)
            % Sigma = 0.1 allows variation while keeping drives positive
            if flg_setpointHetero
                pop(iPop).target = exp(log(pop(2).theta_init + 10) + 0.1 * randn(pop(iPop).N, 1));
            end

        case 'output_centric'
            
            % homogenous setpoint
            pop(iPop).target = pop(iPop).r_init;
            
            % Lognormal targets for firing rates (Median = r_init)
            % Sigma = 0.5 creates a realistic heavy-tailed firing rate distribution
            if flg_setpointHetero
                pop(iPop).target = exp(log(pop(iPop).r_init) + 0.5 * randn(pop(iPop).N, 1));
            end
    end

    pop(iPop).alpha = 0.005;    % Rate of learning
end


%% ========================================================================
%  CONNECTIVITY
%  ========================================================================

% Pre-allocate Connectivity
% Connectivity is split to avoid 'if' statements in the loop.
% W_exc handles all positive drives. W_inh handles all negative drives.
W_exc = zeros(nUnits, nUnits);
W_inh = zeros(nUnits, nUnits);

% Generate Connectivity (ISN Regime)
% We generate weights block-by-block, then stitch them into W_exc or W_inh
W_mean = 0.1;
W_sig  = 0.04;

% Hardcoded scaling factors for ISN balance (Soldado-Magraner et al.)
% J_matrix(post, pre)
J_scaling = [5.00, 1.52;    % E<-E, E<-I
    10.0, 2.25];            % I<-E, I<-I

for iPop = 1:nPop
    for jPop = 1:nPop

        % Generate Block Weights (Scaled by sqrt(N))
        J_val = J_scaling(iPop,jPop) / sqrt(pop(jPop).N);
        W_block = J_val * (W_mean + W_sig * randn(pop(iPop).N, pop(jPop).N));
        W_block(W_block < 0) = 0; % Dale's Law enforcement (no negative magnitudes)

        % Remove self-connections if diagonal block
        if iPop == jPop
            W_block(1:pop(iPop).N + 1:end) = 0;
        end

        % Correct Matrix based on Presynaptic Type
        if pop(jPop).synType == 1
            % Excitatory: Goes to W_exc
            W_exc(pop(iPop).idx, pop(jPop).idx) = W_block;
        else
            % Inhibitory: Goes to W_inh
            % Magnitude stored as positive. Physics loop subtracts it.
            W_inh(pop(iPop).idx, pop(jPop).idx) = W_block;
        end
    end
end


%% ========================================================================
%  SIMULATION LOOP
%  ========================================================================

% Derived constants
trial_steps = floor(trial_dur / dt);
sense_steps = floor(sense_dur / dt);
sense_start = trial_steps - sense_steps + 1;

% Initialize the graphics script
if flg_graphics
    ifrh_graphics;
end

for iTrial = 1:nTrial

    % PERTURBATIONS
    % Apply external input to each population
    ext(:) = 0;
    for iPop = 1:nPop
        p = pop(iPop).perturb; % [Start, End, Amp]
        if iTrial >= p(1) && iTrial <= p(2)
            ext(pop(iPop).idx) = p(3);
        end
    end

    %  ====================================================================
    % FAST DYNAMICS (Physics Engine)

    % Initialize noise
    noise_trial = randn(nUnits, trial_steps);

    % Initialize Accumulators for Learning Rules
    sense_out = zeros(nUnits, 1);
    sense_in = zeros(nUnits, 1);

    for iStep = 1:trial_steps

        % Noise Update
        noise = noise + (dt/OU_tau) * (OU_mu - noise) + ...
            sqrt(2*dt/OU_tau) * OU_sigma * noise_trial(:, iStep);

        % Calculate Synaptic Drives
        Drive_Exc = W_exc * r;
        Drive_Inh = W_inh * r;

        % Total Membrane Drive
        % (Exc - Inh + External + Noise)
        Total_Drive = Drive_Exc - Drive_Inh + ext + noise;

        % Update Firing Rates (Euler)
        dr = (-r + (gain .* max(0, Total_Drive - theta))) ./ tau;
        r = r + dr * dt;

        % Data Accumulation (For error calculation)
        % Only accumulate during the measurement window
        if iStep >= sense_start
            % For output_centric: Track Rate
            sense_out = sense_out + r;

            % For Input-Centric: Track Excitatory Drive + External Input
            sense_in = sense_in + (Drive_Exc + ext);
        end
    end

    % Calculate Averages
    avg_out = sense_out / sense_steps;
    avg_in = sense_in / sense_steps;

    %  ====================================================================
    % SLOW DYNAMICS (Homeostatic Plasticity)
    % Iterate through populations to apply specific rules using Index Maps

    for iPop = 1:nPop

        popIdx = pop(iPop).idx;

        % Store History (For Graphics)
        pop(iPop).hist.r(iTrial, :)       = avg_out(popIdx);
        pop(iPop).hist.theta(iTrial, :)   = theta(popIdx);
        pop(iPop).hist.avg_in(iTrial, :)  = avg_in(popIdx);

        % Calculate Error at steady state 
        switch pop(iPop).learningRule

            case 'input_centric'
                err = avg_in(popIdx) - pop(iPop).target;

            case 'output_centric'
                err = avg_out(popIdx) - pop(iPop).target;

            case 'none'
                err = 0;
        end

        % Apply Learning Rule; Update Global Thresholds
        % If measured > target -> Error > 0 -> Theta increases -> Rate drops
        theta(popIdx) = theta(popIdx) + pop(iPop).alpha * err;
    end

    % VISUALIZATION
    if flg_graphics && (mod(iTrial, graphics_step) == 0 || iTrial == 1)
        ifrh_graphics;
    end
end





%% ========================================================================
%  NOTE: INPUT-CENTRIC SENSING IN ISN REGIMES
%  ========================================================================
%  1. The ISN Paradox (Why Net Input Fails):
%     In an Inhibition-Stabilized Network (ISN), recurrent excitation is
%     strong enough to be unstable, requiring rapid, powerful inhibitory
%     feedback to balance it. Paradoxically, this means that as the
%     excitatory firing rate (E) increases, it recruits such strong
%     inhibition (I) that the *net* synaptic input (W_ee*E - W_ei*I)
%     actually *decreases*.
%
%  2. Control Loop Instability:
%     If the homeostat were to sense Net Input, this paradoxical relationship
%     would invert the error signal. A high firing rate would be read as
%     "low input," causing the homeostat to maladaptively lower thresholds,
%     driving the network into a runaway high-activity state.
%
%  3. The Solution (Excitatory/Metabolic Sensing):
%     To ensure a monotonic relationship between activity and the error signal,
%     this model tracks only Excitatory Drive + External Input (W_ee*E + I_ext).
%     This variable represents the "Metabolic Cost" or "Required Input"
%     (Inverse Transfer Function, f^-1(E)).
%
%  4. Biological Basis:
%     This maps to mitochondrial calcium accumulation. Mitochondria act as
%     metabolic integrators, accumulating calcium primarily through
%     excitatory influx (energy demand) rather than inhibitory currents.
%     This allows the cell to regulate its "budget" (intrinsic excitability)
%     against a stable metabolic setpoint, ensuring stability even in
%     strongly coupled ISN regimes.
%  ========================================================================


%% ========================================================================
%  NOTE: TRIAL STRUCTURE & TIMESCALE SEPARATION
%  ========================================================================
%  This simulation employs a trial-based structure to explicitly separate
%  the fast timescales of neural dynamics (milliseconds) from the slow
%  timescales of homeostatic plasticity.
%
%  1. Separation of Timescales: Each trial runs long enough for the fast
%     neural variables (rates, potentials) to settle into a stable
%     attractor (Fixed Point). Homeostatic adjustments to thresholds are
%     only applied *between* trials. This prevents the slow homeostat from
%     fighting the necessary transient dynamics of the network.
%
%  2. Steady-State Averaging: The error signal is calculated by averaging
%     activity over the final 1.0s of the trial. Mathematically, this
%     moving average acts as a low-pass filter, similar to a differential
%     equation for calcium accumulation. By ignoring the initial transient
%     and smoothing spiking noise, we ensure the homeostatic rule responds
%     only to the global equilibrium state. This aligns with the proposed
%     biological mechanism of MCU-mediated Ca2+ influx, which operates on
%     timescales far slower than individual action potentials.
%  ========================================================================