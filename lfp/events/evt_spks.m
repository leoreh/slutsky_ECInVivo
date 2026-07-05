function evtSpks = evt_spks(spkTimes, evtTimes, ctrlTimes, peakTime, varargin)
% EVT_SPKS Analyzes spiking rate modulation during events.
%
%   evtSpks = EVT_SPKS(spkTimes, evtTimes, ctrlTimes, peakTime, varargin)
%
%   SUMMARY:
%       Calculates scalar modulation metrics comparing Event vs Control periods.
%       1. Instantaneous Firing Rates (FR) per event (Fixed Window).
%       2. Mean FR differences (EventMu vs ControlMu).
%       3. Z-Scored Gain.
%       4. Statistical Significance (Wilcoxon Sign-Rank / Rank-Sum).
%       5. Rank Order (Mean and Variance).
%       6. Population Center of Mass (CoM) per event.
%
%   INPUTS:
%       spkTimes    - (Cell) {N_units x 1} Spike times [s].
%       evtTimes   - (Mat)  [N x 2] Event start/end times [s].
%       ctrlTimes   - (Mat)  [N x 2] Control start/end times [s].
%       peakTime    - (Vec)  [N x 1] Peak times of events [s].
%       varargin    - Parameter/Value pairs:
%           'basepath' - (Char) Save location. (Default: pwd).
%           'flgSave'  - (Log)  Save output? (Default: true).
%           'unitType' - (Cat)  [N_units x 1] Categorical array of unit types.
%
%   OUTPUTS:
%       evtSpks    - (Struct) Stats structure with [N_units x 1] fields.
%
%   HISTORY:
%       Updated: 26 Jan 2026
%

% =========================================================================
%  ARGUMENTS
% =========================================================================
p = inputParser;
addRequired(p, 'spkTimes', @iscell);
addRequired(p, 'evtTimes', @isnumeric);
addRequired(p, 'ctrlTimes', @isnumeric);
addRequired(p, 'peakTime', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'unitType', [], @(x) iscategorical(x) || iscell(x) || isempty(x));
addParameter(p, 'winFxd', 0.020, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, spkTimes, evtTimes, ctrlTimes, peakTime, varargin{:});

basepath  = p.Results.basepath;
flgSave   = p.Results.flgSave;
peakTime  = p.Results.peakTime;
unitType  = p.Results.unitType;
winFxd    = p.Results.winFxd;

[~, basename] = fileparts(basepath);
savefile = fullfile(basepath, [basename, '.evtSpks.mat']);

% Fixed window for asym / com (winFxd param; default 0.020 s matches events)

% =========================================================================
%  FR METRICS
% =========================================================================

nUnits = length(spkTimes);

% Get Spike COUNTS (Full Event Duration)
evtCounts = times2rate(spkTimes, 'winCalc', evtTimes, 'binsize', Inf, 'c2r', false);
ctrlCounts = times2rate(spkTimes, 'winCalc', ctrlTimes, 'binsize', Inf, 'c2r', false);

% Calculate RATES (Fixed Window)
winEvt = [peakTime - winFxd, peakTime + winFxd];
midCtrl = mean(ctrlTimes, 2);
winCtrl = [midCtrl - winFxd, midCtrl + winFxd];
evtRates = times2rate(spkTimes, 'winCalc', winEvt, 'binsize', Inf, 'c2r', true);
ctrlRates = times2rate(spkTimes, 'winCalc', winCtrl, 'binsize', Inf, 'c2r', true);

% Mean Rates
frEvt = mean(evtRates, 2, 'omitnan');
frCtrl = mean(ctrlRates, 2, 'omitnan');
sdCtrl = std(ctrlRates, [], 2, 'omitnan');

% Spike Count Stats
cEvt = mean(evtCounts, 2, 'omitnan');
pFire = mean(evtCounts > 0, 2, 'omitnan');

% Modulation Metrics
frZ = (frEvt - frCtrl) ./ sdCtrl;
frZ(sdCtrl < eps) = NaN; % Avoid infs

frSum = frEvt + frCtrl;
frMod = (frEvt - frCtrl) ./ frSum;
frMod(frSum < eps) = NaN;

% Conditional Rate / Count (intensity when active).
% Essentially equal to frEvt ./ pFire.
% Prepare matrix with NaNs where rate is 0 (inactive)
activeMat = evtRates;
activeMat(activeMat == 0) = NaN;
frActive = mean(activeMat, 2, 'omitnan');

activeMat = evtCounts;
activeMat(activeMat == 0) = NaN;
cActive = mean(activeMat, 2, 'omitnan');

% =========================================================================
%  STATISTICAL TEST
% =========================================================================

pVal = nan(nUnits, 1);
flgPair = size(evtRates, 2) == size(ctrlRates, 2);

% Pre-calculate nanflag to speed up loop skip
goodIdx = ~all(isnan(evtRates), 2) & sum(evtCounts, 2) > 0;

for iUnit = 1:nUnits
    if ~goodIdx(iUnit), continue; end

    if flgPair
        pVal(iUnit) = signrank(evtRates(iUnit, :), ctrlRates(iUnit, :));
    else
        pVal(iUnit) = ranksum(evtRates(iUnit, :), ctrlRates(iUnit, :));
    end
end
h0 = pVal < 0.05;

% =========================================================================
%  CENTER OF MASS (PER UNIT)
% =========================================================================
% Average time of spikes relative to event peak

com = nan(nUnits, 1);

peakTime = peakTime(:);
nEvt = length(peakTime);

for iUnit = 1:nUnits
    uSpks = spkTimes{iUnit};
    if isempty(uSpks), continue; end

    % Find index of nearest peak for each spike
    % We use interp1 with 'nearest' to map spike times to the "index" of the peak
    idx = interp1(peakTime, 1:nEvt, uSpks, 'nearest', 'extrap');
    nearestPeaks = peakTime(idx);
    relTime = uSpks - nearestPeaks;

    % Filter by Window
    isValid = abs(relTime) <= winFxd;

    % Calculate Mean
    if any(isValid)
        com(iUnit) = mean(relTime(isValid)) * 1000;     % (ms)
    end
end

% =========================================================================
%  ASYMMETRY (PER UNIT)
% =========================================================================

% Asymmetry (Per Unit)
% Define Intervals (Fixed Window)
timesMid = peakTime;
timesPre  = [timesMid - winFxd, timesMid];
timesPost = [timesMid, timesMid + winFxd];

% Durations
durPre  = timesPre(:,2) - timesPre(:,1);
durPost = timesPost(:,2) - timesPost(:,1);

% Spike Counts (all units)
cPre  = times2rate(spkTimes, 'winCalc', timesPre, 'binsize', Inf, 'c2r', false);
cPost = times2rate(spkTimes, 'winCalc', timesPost, 'binsize', Inf, 'c2r', false);

% Rates
frPre  = sum(cPre, 2) ./ sum(durPre);
frPost = sum(cPost, 2) ./ sum(durPost);

% Asymmetry Index (Unit)
% (Pre - Post) / (Pre + Post)
frSum = frPre + frPost;
asym = (frPre - frPost) ./ frSum;
asym(frSum < eps) = NaN;

% =========================================================================
%  PER UNITTYPE
% =========================================================================

% Initialize Output Containers
rankMean   = nan(nUnits, 1);
rankVar    = nan(nUnits, 1);

% Prepare Iteration (Global vs Types)
if isempty(unitType)
    iterNames = {'Global'};
    iterMasks = {true(nUnits, 1)};
else
    uTypes = unique(unitType);
    uTypes(isundefined(uTypes)) = [];
    uTypes(uTypes == 'Other') = [];

    iterNames = cellstr(uTypes);
    iterMasks = cell(length(uTypes), 1);
    for iType = 1:length(uTypes)
        iterMasks{iType} = unitType == uTypes(iType);
    end
end

% Initialize Events Structure
evtSpks.events = struct();

% Loop
for iIter = 1:length(iterNames)
    currName = iterNames{iIter};
    currMask = iterMasks{iIter};

    if sum(currMask) == 0, continue; end

    % --- Rank ---
    subSpks = spkTimes(currMask);
    [subMean, subVar] = evt_rankOrder(subSpks, evtTimes);

    % Fill Global Arrays
    rankMean(currMask)   = subMean;
    rankVar(currMask)    = subVar;

    % --- Population Stats (Per Event) ---

    % Counts for this subset
    subCounts = evtCounts(currMask, :);

    % Fraction Participation
    % (Active Units in Group / Total Units in Group)
    nActive = sum(subCounts > 0, 1); % [1 x nEvt]
    nTotal  = sum(currMask);
    currFrac = (nActive ./ nTotal)';

    % Asymmetry (Per Event)
    subPre  = cPre(currMask, :);
    subPost = cPost(currMask, :);

    ratePre  = sum(subPre, 1) ./ durPre';
    ratePost = sum(subPost, 1) ./ durPost';

    frSum = ratePre + ratePost;
    currAsym = (ratePre - ratePost) ./ frSum;
    currAsym(frSum < eps) = NaN;

    % --- Center of Mass (Per Event) ---
    % Calculate the center of mass of spikes relative to the event peak.
    % This is done by aggregating all spikes from the current unit group.
    currCom = nan(size(evtTimes, 1), 1);

    % Collect all spikes from the current unit selection
    grpSpks = vertcat(subSpks{:});

    % Find nearest event peak for each spike
    peakIdx = interp1(peakTime, 1:length(peakTime), grpSpks, 'nearest', 'extrap');

    % Calculate relative time
    tRel = grpSpks - peakTime(peakIdx);

    % Filter spikes within fixed window
    inWin = abs(tRel) <= winFxd;

    if any(inWin)
        valIdx = peakIdx(inWin);
        valRel = tRel(inWin);

        % Calculate Mean CoM per Event (in ms)
        currCom = accumarray(valIdx, valRel, [length(peakTime), 1], @mean, NaN) * 1000;
    end

    % --- Store Results ---
    fn = matlab.lang.makeValidName(currName);
    evtSpks.events.(fn).frac = currFrac;
    evtSpks.events.(fn).asym = currAsym';
    evtSpks.events.(fn).com  = currCom;
end



% =========================================================================
%  OUTPUT & SAVE
% =========================================================================

% Pack results
evtSpks.cEvt     = cEvt;
evtSpks.frEvt    = frEvt;
evtSpks.frCtrl    = frCtrl;
evtSpks.frZ       = frZ;
evtSpks.frMod     = frMod;
evtSpks.pFire     = pFire;
evtSpks.frActive  = frActive;
evtSpks.cActive   = cActive;
evtSpks.pVal      = pVal;
evtSpks.h0        = h0;
evtSpks.rankMean  = rankMean;
evtSpks.rankVar   = rankVar;
evtSpks.com       = com;

evtSpks.frPre    = frPre;
evtSpks.frPost   = frPost;
evtSpks.asym = asym;

% Save
if flgSave
    save(savefile, 'evtSpks', '-v7.3');
end

end     % EOF