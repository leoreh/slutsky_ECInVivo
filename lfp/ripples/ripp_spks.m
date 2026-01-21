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
%   varargin        'muTimes' (cell array), 'basepath', 'flgPlot', 'flgSave'.
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
addParameter(p, 'muTimes', [], @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'mapDur', [-0.075 0.075], @isnumeric);
parse(p, rippTimes, spkTimes, peakTimes, ctrlTimes, varargin{:});

rippTimes = p.Results.rippTimes;
spkTimes = p.Results.spkTimes;
peakTimes = p.Results.peakTimes;
ctrlTimes = p.Results.ctrlTimes;
muTimes = p.Results.muTimes;
basepath = p.Results.basepath;
flgPlot = p.Results.flgPlot;
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
fsLfp = 1250;
nBinsMap = floor(fsLfp * diff(mapDur) / 2) * 2 + 1;

% Output Struct
rippSpks = struct();
rippSpks.info.mapDur = mapDur;
rippSpks.info.nBinsMap = nBinsMap;
rippSpks.su = struct('pVal', nan(nUnits, 1), 'sigMod', false(nUnits, 1), ...
    'frGain', nan(nUnits, 1), 'frModulation', nan(nUnits, 1), ...
    'frPrct', nan(nUnits, 1), 'frRipp', nan(nUnits, 1), 'frRand', nan(nUnits, 1), ...
    'rippMap', [], 'ctrlMap', []);
rippSpks.mu = struct('rippMap', [], 'ctrlMap', []);


%% ========================================================================
%  SU ANALYSIS
%  ========================================================================
fprintf('Analyzing %d Single Units...\n', nUnits);

% Firing Rates
rippRates = times2rate(spkTimes, 'winCalc', rippTimes, 'binsize', Inf);
ctrlRates = times2rate(spkTimes, 'winCalc', ctrlTimes, 'binsize', Inf);

rippSpks.su.rippRates = rippRates;
rippSpks.su.ctrlRates = ctrlRates;

% Maps & Stats
rippSpks.su.rippMap = nan(nUnits, nRipples, nBinsMap);
rippSpks.su.ctrlMap = nan(nUnits, nRipples, nBinsMap);

for iUnit = 1:nUnits
    % Stats
    rRate = rippRates(iUnit, :);
    cRate = ctrlRates(iUnit, :);

    if any(~isnan(rRate)) && any(~isnan(cRate))
        [pVal, h] = signrank(rRate, cRate);
        rippSpks.su.pVal(iUnit) = pVal;
        rippSpks.su.sigMod(iUnit) = h;

        mC = mean(cRate, 'omitnan');
        sdC = std(cRate, 0, 'omitnan');
        mR = mean(rRate, 'omitnan');

        rippSpks.su.frGain(iUnit) = (mR - mC) / sdC;
        rippSpks.su.frModulation(iUnit) = (mR - mC) / (mR + mC);
        rippSpks.su.frPrct(iUnit) = (mR - mC) / mC * 100;
        rippSpks.su.frRipp(iUnit) = mR;
        rippSpks.su.frRand(iUnit) = mC;
    end

    % Maps
    rippSpks.su.rippMap(iRipp, :, :) = sync_spksMap(spkTimes{iRipp}, ...
        peakTimes, mapDur, nBinsMap);
    rippSpks.su.ctrlMap(iRipp, :, :) = sync_spksMap(spkTimes{iRipp}, ...
        diff(ctrlTimes, 1, 2) / 2, mapDur, nBinsMap);
end


%% ========================================================================
%  MU ANALYSIS
%  ========================================================================
if ~isempty(muTimes)
    fprintf('Analyzing Multi-Unit Activity...\n');
    rippSpks.mu.rippMap = sync_spksMap(muTimes, peakTimes, mapDur, nBinsMap);
    rippSpks.mu.ctrlMap = sync_spksMap(muTimes, diff(ctrlTimes, 1, 2) / 2, mapDur, nBinsMap);
else
    rippSpks.mu = [];
end

if ~isempty(rippSpks.mu)
    rippRatesMU = times2rate(muTimes, 'winCalc', rippTimes, 'binsize', Inf);
    ctrlRatesMU = times2rate(muTimes, 'winCalc', ctrlTimes, 'binsize', Inf);
    rippSpks.mu.muGain = (rippRatesMU - mean(ctrlRatesMU, 'omitnan')) ./ std(ctrlRatesMU, 'omitnan');
end

%% ========================================================================
%  SAVE & PLOT
%  ========================================================================
if flgSave
    save(rippSpksFile, 'rippSpks', '-v7.3');
    fprintf('Saved: %s\n', rippSpksFile);
end

if flgPlot
    ripp_plotSpks(rippSpks)
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
