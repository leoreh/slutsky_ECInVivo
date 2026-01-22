function rippSpks = ripp_spks(spkTimes, rippTimes, ctrlTimes, varargin)
% RIPP_SPKS Analyzes spiking statistics (rates/gain) during ripples.
%
% SUMMARY:
%   Calculates scalar metrics: Firing Rates, Gain, Modulation, and P-values.
%   Does NOT generate PETH maps (see ripp_spkMaps).
%
% INPUT:
%   rippTimes   - [N x 2] start/end times (s).
%   spkTimes    - {Nunits x 1} spike times (s).
%   ctrlTimes   - [N x 2] control event times (s).
%
% OUTPUT:
%   rippSpks    - Structure with fields: .frRipp, .frRand, .frZ, .pVal, etc.

% =========================================================================
%  ARGUMENTS
% =========================================================================
p = inputParser;
addRequired(p, 'spkTimes', @iscell);
addRequired(p, 'rippTimes', @isnumeric);
addRequired(p, 'ctrlTimes', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
parse(p, spkTimes, rippTimes, ctrlTimes, varargin{:});

basepath = p.Results.basepath;
flgSave = p.Results.flgSave;

% =========================================================================
%  CALCULATION
% =========================================================================
[~, basename] = fileparts(basepath);
savefile = fullfile(basepath, [basename, '.rippSpks.mat']);

nUnits = length(spkTimes);

% Initialize Output
rippSpks = struct();
rippSpks.pVal = nan(nUnits, 1);
rippSpks.sigMod = false(nUnits, 1);
rippSpks.frZ = nan(nUnits, 1);
rippSpks.frMod = nan(nUnits, 1);
rippSpks.frRipp = nan(nUnits, 1);
rippSpks.frRand = nan(nUnits, 1);

% Calculate Instantaneous Rates (Events x Units)
% 'binsize', Inf -> Returns one rate per event
rippSpks.rippRates = times2rate(spkTimes, 'winCalc', rippTimes, 'binsize', Inf);
rippSpks.ctrlRates = times2rate(spkTimes, 'winCalc', ctrlTimes, 'binsize', Inf);

% Per Unit Statistics
for iUnit = 1:nUnits
    
    rRate = rippSpks.rippRates(iUnit, :);
    cRate = rippSpks.ctrlRates(iUnit, :);
    
    if any(isnan(rRate)) || any(isnan(cRate))
        continue;
    end
    
    % Means & Std
    rMu = mean(rRate, 'omitnan');
    cMu = mean(cRate, 'omitnan');
    cSigma = std(cRate, 'omitnan');
    
    % Store Descriptive Stats
    rippSpks.frRipp(iUnit) = rMu;
    rippSpks.frRand(iUnit) = cMu;
    
    % Metrics
    if cSigma > 0
        rippSpks.frZ(iUnit) = (rMu - cMu) / cSigma;
    else
        rippSpks.frZ(iUnit) = NaN; % Avoid infinite gain
    end
    
    % Modulation Index (-1 to 1)
    if (rMu + cMu) > 0
        rippSpks.frMod(iUnit) = (rMu - cMu) / (rMu + cMu);
    end

    % Statistical Significance (Wilcoxon Sign Rank)
    % Paired test comparing distribution of ripple rates vs control rates
    % Note: Requires equal number of events. If N_ripp != N_ctrl, 
    % we use ranksum (unpaired) 
    if length(rRate) == length(cRate)
        rippSpks.pVal(iUnit) = signrank(rRate, cRate);
    else
        rippSpks.pVal(iUnit) = ranksum(rRate, cRate);
    end
    
    rippSpks.sigMod(iUnit) = rippSpks.pVal(iUnit) < 0.05;
end

% =========================================================================
%  SAVE
% =========================================================================
if flgSave
    save(savefile, 'rippSpks', '-v7.3');
    fprintf('Saved spike stats: %s\n', savefile);
end

end