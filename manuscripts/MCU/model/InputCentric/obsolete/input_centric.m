%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL: Input-Centric Homeostasis 
% 
% A flexible, population-based implementation of homeostatic plasticity 
% in an Inhibition-Stabilized Network (ISN).
%
% ARCHITECTURE:
% - Data is stored in a 'pop' structure array.
% - Connectivity is a cell array of weight matrices.
% - Equations are vectorized to handle N populations dynamically.
%
% DYNAMICS:
% 1. Fast:  tau * dr/dt = -r + F(W*r + I_ext)
% 2. Slow:  theta += alpha * (SensedMetric - Setpoint)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; rng(44); 

%% GLOBAL SIMULATION SETTINGS
% -------------------------------------------------------------------------
dt = 0.0001;                 % Integration step (s)
trial_dur = 2;               % Duration per trial (s)
trial_steps = trial_dur/dt;
nTrial = 500;
tau_measure = 1.0;           % Window for homeostatic averaging (s)

% Flags
flg_graphics = true;
flg_homeostatic = true;
graphics_step = 5;

% Noise (Ornstein-Uhlenbeck)
OU.tau = 0.1; 
OU.mu = 0; 
OU.sigma = 4.0; 

% Activation Function (ReLU)
Func_Act = @(x, gain, th) gain * max(0, x - th);


%% NEURON POPULATIONS 
% -------------------------------------------------------------------------

% Options for learning rule:
% - 'input_centric'
% - 'output_centric'
% - 'none'

% EXCITATORY (E)
pop(1).name = 'Exc';
pop(1).synType = +1;            % +1 Excitatory
pop(1).color = [0.2, 0.7, 0.2]; 
pop(1).N = 80;
pop(1).tau = 0.010;
pop(1).gain = 1;
pop(1).theta_init = 4.8;
pop(1).max_rate = 100;

% Plasticity Settings
pop(1).learningRule = 'input_centric'; 
pop(1).alpha = 0.005;           % Learning rate
pop(1).target = 20;             % Target (Hz for output, a.u. for input)            

% Perturbation Settings [OnsetTrial, OffsetTrial, Amplitude]
pop(1).perturb = [550, 550, 5.0]; 


% INHIBITORY (I) 
pop(2).name = 'Inh';
pop(2).synType = -1;            % -1 Inhibitory
pop(2).color = [0.8, 0.2, 0.2]; 
pop(2).N = 20;
pop(2).tau = 0.002;             % Fast inhibition for ISN stability
pop(2).gain = 4;
pop(2).theta_init = 15;
pop(2).max_rate = 250;

% Plasticity Settings 
pop(2).learningRule = 'none';   
pop(2).target = 10.0;
pop(2).alpha = 0.005;

% Perturbation Settings [OnsetTrial, OffsetTrial, Amplitude]
pop(2).perturb = [550, 550, 5.0];


%% STATE INITIALIZATION
% -------------------------------------------------------------------------

nPop = length(pop);

for iPop = 1:nPop
    % Current State
    pop(iPop).r = rand(pop(iPop).N, 1) * 5;     % Firing Rate
    pop(iPop).theta = ones(pop(iPop).N, 1) * pop(iPop).theta_init;
    pop(iPop).noise = zeros(pop(iPop).N, 1);
    
    % History Pre-allocation (Structure of arrays)
    pop(iPop).hist.meanRate = NaN(nTrial, 1);
    pop(iPop).hist.meanTheta = NaN(nTrial, 1);
    pop(iPop).hist.meanSensed = NaN(nTrial, 1); % The variable being regulated
    pop(iPop).hist.allRate = NaN(nTrial, min(pop(iPop).N, 10)); % Save subset
    pop(iPop).hist.allTheta = NaN(nTrial, min(pop(iPop).N, 10));
    
    % Plotting Indices (Save subset to save RAM)
    pop(iPop).idxPlot = round(linspace(1, pop(iPop).N, min(pop(iPop).N, 10)));
end

%% CONNECTIVITY (ISN REGIME)
% -------------------------------------------------------------------------
% We use a Cell Matrix W{post, pre}.
% W{i,j} is the matrix connecting Population j -> Population i.
% Note: We strictly use positive weights. 
% Signs are applied during the dynamics loop (Exc+, Inh-).

W = cell(nPop, nPop);

% ISN Scaling Factors. Adapted from (Soldado, PNAS, 2022)
J_EE = 5.00 / sqrt(pop(1).N);
J_EI = 1.52 / sqrt(pop(2).N);
J_IE = 10.0 / sqrt(pop(1).N);
J_II = 2.25 / sqrt(pop(2).N);

W_mean = 0.1; 
W_sig  = 0.04;

% Generate Random Connectivity
W{1,1} = J_EE * (W_mean + W_sig * randn(pop(1).N, pop(1).N)); % E <- E
W{1,2} = J_EI * (W_mean + W_sig * randn(pop(1).N, pop(2).N)); % E <- I
W{2,1} = J_IE * (W_mean + W_sig * randn(pop(2).N, pop(1).N)); % I <- E
W{2,2} = J_II * (W_mean + W_sig * randn(pop(2).N, pop(2).N)); % I <- I

% Remove self-connections
W{1,1}(1:pop(1).N+1:end) = 0; 
W{2,2}(1:pop(2).N+1:end) = 0;





%% MAIN SIMULATION LOOP
% -------------------------------------------------------------------------

% Initialize GUI
if flg_graphics
    input_centric_graphics; 
end

for iTrial = 1:nTrial
    
    % PERTURBATIONS
    for iPop = 1:nPop
        p = pop(iPop).perturb;
        if iTrial >= p(1) && iTrial <= p(2)
            pop(iPop).extInput = p(3);
        else
            pop(iPop).extInput = 0;
        end
    end
    
    % FAST DYNAMICS (Within-Trial)
    % ---------------------------------------------------------------------
    % We need temporary storage for the trial to calculate averages
    trialData = repmat(struct('r', [], 'sensed', []), nPop, 1);
    for iPop = 1:nPop
        trialData(iPop).r = zeros(pop(iPop).N, trial_steps);
        trialData(iPop).sensed = zeros(pop(iPop).N, trial_steps);
    end
    
    for iStep = 1:trial_steps
        
        % Calculate Drives for all populations
        % We do this first so updates are synchronous
        for iPop = 1:nPop
            
            % Update Noise (OU Process)
            pop(iPop).noise = pop(iPop).noise + (dt/OU.tau)*(OU.mu - pop(iPop).noise) + ...
                            sqrt(2*dt/OU.tau)*OU.sigma*randn(pop(iPop).N, 1);
            
            % Synaptic Summation
            exc_drive = 0;
            inh_drive = 0;

            % Loop through sources
            for jPop = 1:nPop
                % Calculate raw input magnitude (Weights are positive)
                syn_input = W{iPop,jPop} * pop(jPop).r;

                % Sort into drives based on Synaptic Type (+1 or -1)
                if pop(jPop).synType == 1
                    exc_drive = exc_drive + syn_input;
                elseif pop(jPop).synType == -1
                    inh_drive = inh_drive + syn_input;
                end
            end

            % Total Membrane Drive (Exc - Inh + Ext + Noise)
            pop(iPop).drive = exc_drive - inh_drive + pop(iPop).extInput + pop(iPop).noise;

            % Store "Sensed" Variable (For Plasticity)
            % Input-centric specifically senses Excitatory drive + Ext
            pop(iPop).sensed_now = exc_drive + pop(iPop).extInput + pop(iPop).noise;
        end
        
        % Update Rates (Euler Integration)
        for iPop = 1:nPop
            dr = (-pop(iPop).r + Func_Act(pop(iPop).drive, pop(iPop).gain, pop(iPop).theta)) / pop(iPop).tau;
            pop(iPop).r = pop(iPop).r + dr * dt;
            
            % Hard Saturation
            pop(iPop).r(pop(iPop).r > pop(iPop).max_rate) = pop(iPop).max_rate;
            
            % Store Data
            trialData(iPop).r(:, iStep) = pop(iPop).r;
            trialData(iPop).sensed(:, iStep) = pop(iPop).sensed_now;
        end
    end
    
    
    % SLOW DYNAMICS (Homeostatic Plasticity)
    % ---------------------------------------------------------------------
    ss_steps = floor(tau_measure/dt);
    range = (trial_steps - ss_steps + 1) : trial_steps;
    
    for iPop = 1:nPop
        % Calculate trial averages
        avg_r = mean(trialData(iPop).r(:, range), 2);
        avg_sensed = mean(trialData(iPop).sensed(:, range), 2);
        
        % Make available to graphics script
        if iPop == 1
            avg_r_E = avg_r;
            avg_sensed_E = avg_sensed;
        elseif iPop == 2
            avg_r_I = avg_r;
        end
        
        % Record History
        pop(iPop).hist.meanRate(iTrial) = mean(avg_r);
        pop(iPop).hist.meanTheta(iTrial) = mean(pop(iPop).theta);
        pop(iPop).hist.allRate(iTrial, :) = avg_r(pop(iPop).idxPlot)';
        pop(iPop).hist.allTheta(iTrial, :) = pop(iPop).theta(pop(iPop).idxPlot)';
        
        % Plasticity Logic
        if flg_homeostatic && ~strcmpi(pop(iPop).learningRule, 'none')
            
            switch pop(iPop).learningRule
                case 'input_centric'
                    % Regulate Synaptic Drive
                    measured = avg_sensed; 
                    pop(iPop).hist.meanSensed(iTrial) = mean(measured);
                    
                case 'output_centric'
                    % Regulate Firing Rate
                    measured = avg_r;
                    pop(iPop).hist.meanSensed(iTrial) = mean(measured);
            end
            
            % Error Signal (Too high -> Increase Theta)
            err = measured - pop(iPop).target;
            
            % Update Thresholds
            pop(iPop).theta = pop(iPop).theta + pop(iPop).alpha * err;
        end
    end
    
    % VISUALIZATION
    if flg_graphics && (mod(iTrial, graphics_step) == 0 || iTrial == 1)
        input_centric_graphics;
    end
end