function M = cyto2mito(C, varargin)
% MITO2CYTO Calculates mitochondrial calcium accumulation from cytosolic calcium.
%
%   M = MITO2CYTO(C, KD, N) calculates the mitochondrial matrix calcium (M)
%   given cytosolic calcium concentration (C), dissociation constant (KD),
%   and Hill coefficient (N). It first computes the influx (J) using the
%   Hill equation and then integrates it using either a linear filter or
%   Michaelis-Menten efflux kinetics.
%
%   INPUTS:
%       C           - (vec) Cytosolic Calcium Concentration.
%       varargin    - (param/value) Optional parameters:
%                     'Kd'      : (num) Dissociation constant for MCU.
%                     'n'       : (num) Hill coefficient for MCU.
%                     'dt'      : (num) Time step {0.001} (s).
%                     'tauM'    : (num) Matrix decay time constant {20} (s).
%                     'Vmax'    : (num) Max Efflux Rate {0.05} (Conc/s).
%                     'Km'      : (num) Efflux Affinity {1.0} (Conc).
%                     'flgVmax' : (log) Use Michaelis-Menten Efflux {true}.
%                                 If false, uses linear decay (tauM).
%
%   OUTPUTS:
%       M           - (vec) Mitochondrial Matrix Calcium Concentration.
%
%   See also: SPK2CA, SPK2CA_SIM

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'C', @isnumeric);
addParameter(p, 'Kd', 10, @isnumeric);
addParameter(p, 'n', 4, @isnumeric);
addParameter(p, 'dt', 0.001, @isnumeric);
addParameter(p, 'tauM', 20, @isnumeric);
addParameter(p, 'Vmax', 0.0001, @isnumeric);
addParameter(p, 'Km', 0.00005, @isnumeric);
addParameter(p, 'flgVmax', true, @islogical);

parse(p, C, varargin{:});
par = p.Results;

%% ========================================================================
%  COMPUTE INFLUX (J)
%  ========================================================================

% Hill Equation
% J = C^n / (C^n + Kd^n)
Cn = C .^ par.n;
Kn = par.Kd ^ par.n;
J  = Cn ./ (Cn + Kn);

%% ========================================================================
%  COMPUTE MATRIX CALCIUM (M)
%  ========================================================================

% Derived Constants
alphaM    = exp(-par.dt / par.tauM);
gainIn    = (1 - alphaM);      % Scales J to be physically compatible
Vmax_step = par.Vmax * par.dt; % Max efflux per time step

if ~par.flgVmax
    % --- Linear Filter Mode ---
    % M(t) = alpha * M(t-1) + (1-alpha) * J(t)
    % This is the standard 'Leaky Integrator'
    M = filter(gainIn, [1, -alphaM], J);

else
    % --- Michaelis-Menten Efflux Mode ---
    % Influx is linear (driven by J), but Efflux saturates at Vmax

    nBins = length(J);
    M = zeros(size(J), 'like', J);
    M(1) = 0; % Initialize State

    for t = 2:nBins
        % Influx
        In = J(t) * gainIn;

        % Efflux (Saturating)
        prevM = M(t-1);
        Out = (Vmax_step * prevM) / (par.Km + prevM);

        % Update State
        M(t) = prevM + In - Out;

        % Biological Boundary
        if M(t) < 0, M(t) = 0; end
    end
end

end     % EOF
