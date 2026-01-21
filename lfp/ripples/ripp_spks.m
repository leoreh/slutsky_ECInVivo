function rippSpks = ripp_spks(rippTimes, spkTimes, peakTimes, ctrlTimes, varargin)
% RIPP_SPKS Analyzes spiking activity relative to ripple events.
%
% SUMMARY:
% Analyzes single-unit (SU) and multi-unit (MU) activity during ripples.
% Calculates rates, modulation statistics, and generates PETH maps.
%
% INPUT:
%   rippTimes       [N x 2] start/end times (s).
%   spkTimes        {Nunits x 1} spike times (s).
%   peakTimes       [N x 1] ripple peak times (s).
%   ctrlTimes       [N x 2] control event times (s).
%   varargin        'basepath', 'flgPlot', 'flgSave'.
%
% OUTPUT:
%   rippSpks        Structure of results.
%   * Saves 'basename.rippSpks.mat'.
%
% DEPENDENCIES: SubtractIntervals, Restrict, times2rate.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'rippTimes', @isnumeric);
addRequired(p, 'spkTimes', @iscell);
addRequired(p, 'peakTimes', @isnumeric);
addRequired(p, 'ctrlTimes', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'mapDur', [-0.075 0.075], @isnumeric);
parse(p, rippTimes, spkTimes, peakTimes, ctrlTimes, varargin{:});

rippTimes = p.Results.rippTimes;
spkTimes = p.Results.spkTimes;
peakTimes = p.Results.peakTimes;
ctrlTimes = p.Results.ctrlTimes;
basepath = p.Results.basepath;
flgSave = p.Results.flgSave;
mapDur = p.Results.mapDur;

%% ========================================================================
%  INITIALIZATION
%  ========================================================================
cd(basepath);
[~, basename] = fileparts(basepath);
rippSpksFile = fullfile(basepath, [basename, '.rippSpks.mat']);

nUnits = length(spkTimes);
nRipples = size(rippTimes, 1);

% Map Params
nBinsMap = ceil(1250 * diff(mapDur));

% Output Struct
rippSpks = struct();
rippSpks.pVal = nan(nUnits, 1);
rippSpks.sigMod = false(nUnits, 1);
rippSpks.frZ = nan(nUnits, 1);
rippSpks.frModulation = nan(nUnits, 1);
rippSpks.frPrct = nan(nUnits, 1);
rippSpks.frRipp = nan(nUnits, 1);
rippSpks.frRand = nan(nUnits, 1);
rippSpks.rippMap = nan(nUnits, nRipples, nBinsMap);
rippSpks.ctrlMap = nan(nUnits, nRipples, nBinsMap);

% Firing Rates
rippSpks.rippRates = times2rate(spkTimes, 'winCalc', rippTimes, 'binsize', Inf);
rippSpks.ctrlRates = times2rate(spkTimes, 'winCalc', ctrlTimes, 'binsize', Inf);

% Per Unit Stats
for iUnit = 1:nUnits
    
    % Stats
    rRate = rippSpks.rippRates(iUnit, :);
    cRate = rippSpks.ctrlRates(iUnit, :);
    
    if any(isnan(rRate)) || any(isnan(cRate))
        continue
    end

    cSigma = std(cRate, 0, 'omitnan');   
    cMu = mean(cRate, 'omitnan');
    rMu = mean(rRate, 'omitnan');
    
    % Indices
    rippSpks.frRipp(iUnit) = rMu;
    rippSpks.frRand(iUnit) = cMu;
    rippSpks.frZ(iUnit) = (rMu - cMu) / cSigma;
    rippSpks.frPrct(iUnit) = (rMu - cMu) / cMu * 100;
    rippSpks.frMod(iUnit) = (rMu - cMu) / (rMu + cMu);

    % Rate in ripples significantly higher than random
    rippSpks.pVal(iUnit) = signrank(rRate, cRate);

    % Maps
    rippSpks.rippMap(iUnit, :, :) = sync_spksMap(spkTimes{iUnit}, ...
        peakTimes, mapDur, nBinsMap);
    rippSpks.ctrlMap(iUnit, :, :) = sync_spksMap(spkTimes{iUnit}, ...
        diff(ctrlTimes, 1, 2) / 2, mapDur, nBinsMap);
end

% Handle case of 1 unit (mu)
rippSpks.rippMap = squeeze(rippSpks.rippMap);
rippSpks.ctrlMap = squeeze(rippSpks.ctrlMap);


%% ========================================================================
%  SAVE & PLOT
%  ========================================================================
if flgSave
    save(rippSpksFile, 'rippSpks', '-v7.3');
    fprintf('Saved: %s\n', rippSpksFile);
end


end     % EOF


%% ========================================================================
%  HELPER: SYNC MAP
%  ========================================================================
function syncMap = sync_spksMap(spikeTimes, eventTimes, mapDur, nBinsMap)

nEvents = length(eventTimes);
syncMap = zeros(nEvents, nBinsMap);

if isempty(spikeTimes) || isempty(eventTimes), return; end
spikeTimes = sort(spikeTimes(:));
eventTimes = eventTimes(:);

winStart = eventTimes + mapDur(1);
winEnd = eventTimes + mapDur(2);

searchEdges = [-inf; spikeTimes; inf];

startIdx = discretize(winStart, searchEdges);
[~, endBin] = histc(winEnd, searchEdges);
endIdx = endBin - 1;

validWin = find(endIdx >= startIdx);
if isempty(validWin), return; end

spkInEvt = endIdx(validWin) - startIdx(validWin) + 1;
evtIds = repelem(validWin, spkInEvt);

totalSpk = sum(spkInEvt);
spkIndices = zeros(totalSpk, 1);
c = 1;
for k = 1:length(validWin)
    e = validWin(k);
    n = spkInEvt(k);
    spkIndices(c:c+n-1) = startIdx(e) : endIdx(e);
    c = c+n;
end

spkVals = spikeTimes(spkIndices);
relTimes = spkVals - eventTimes(evtIds);
edges = linspace(mapDur(1), mapDur(2), nBinsMap+1);
bins = discretize(relTimes, edges);

validBin = ~isnan(bins);
finalEvtIds = evtIds(validBin);
finalBins = bins(validBin);

syncMap = accumarray([finalEvtIds, finalBins], 1, [nEvents, nBinsMap]);


end
