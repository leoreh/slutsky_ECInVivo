function maps = ripp_spkPeth(spkTimes, peakTimes, ctrlTimes, varargin)
% RIPP_SPKMAPS Generates spike PETH maps for ripples and control events.
%
% SUMMARY:
%   Calculates 3D maps (Unit x Event x Time) for ripples and controls.
%   Uses fast vectorized binning (discretize/accumarray) for performance.
%
% INPUT:
%   spkTimes    - {Nunits x 1} Cell array of spike times [s].
%   peakTimes   - [N x 1] Ripple peak times [s] (Map centers).
%   ctrlTimes   - [N x 2] Control start/end [s].
%   varargin    - 'basepath', 'mapDur' (default [-0.075 0.075]), 'flgSave'.
%
% OUTPUT:
%   maps        - Structure containing:
%       .ripp     (Nunits x Nripples x Nbins)
%       .ctrl     (Nunits x Ncontrols x Nbins)
%       .bins     (Time vector for x-axis)
%
% DEPENDENCIES: None.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'spkTimes', @iscell);
addRequired(p, 'peakTimes', @isnumeric);
addRequired(p, 'ctrlTimes', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'mapDur', [-0.075 0.075], @isnumeric);
addParameter(p, 'flgSave', false, @islogical);
parse(p, spkTimes, peakTimes, ctrlTimes, varargin{:});

basepath = p.Results.basepath;
mapDur = p.Results.mapDur;
flgSave = p.Results.flgSave;

%% ========================================================================
%  SETUP
%  ========================================================================
[~, basename] = fileparts(basepath);
savefile = fullfile(basepath, [basename, '.rippSpkMaps.mat']);

nUnits = length(spkTimes);
nRipples = length(peakTimes);
nControls = size(ctrlTimes, 1);

% Define Time Bins
% Define Time Bins (Aligned to samples matching ripp_maps)
fs = 1250;
winSamps = round(mapDur(1)*fs) : round(mapDur(2)*fs);
nBins = length(winSamps);
% Edges centered on these samples with 1/fs width
edges = [winSamps - 0.5, winSamps(end) + 0.5] / fs;
timeBins = winSamps / fs;

% Control Centers (Midpoint of control interval)
ctrlCenters = mean(ctrlTimes, 2);

% Initialize Output
maps = struct();
maps.ripp = zeros(nUnits, nRipples, nBins);
maps.ctrl = zeros(nUnits, nControls, nBins);
maps.tstamps = timeBins;

%% ========================================================================
%  GENERATE MAPS
%  ========================================================================
fprintf('Generating Spike Maps for %d units...\n', nUnits);

for iUnit = 1:nUnits
    unitSpks = spkTimes{iUnit};

    if isempty(unitSpks)
        continue;
    end

    % Ripple Map
    maps.ripp(iUnit, :, :) = sync_spksMap(unitSpks, peakTimes, mapDur, edges);

    % Control Map
    maps.ctrl(iUnit, :, :) = sync_spksMap(unitSpks, ctrlCenters, mapDur, edges);
end

%% ========================================================================
%  SAVE
%  ========================================================================
if flgSave
    save(savefile, 'maps', '-v7.3');
    fprintf('Saved spike maps: %s\n', savefile);
end

end

% -------------------------------------------------------------------------
% Helper: Fast PETH Generator (Vectorized)
% -------------------------------------------------------------------------
function syncMap = sync_spksMap(spikeTimes, eventTimes, mapDur, edges)

nBinsMap = length(edges) - 1;
nEvents = length(eventTimes);
syncMap = zeros(nEvents, nBinsMap);

if isempty(spikeTimes) || isempty(eventTimes)
    return;
end

spikeTimes = sort(spikeTimes(:));
eventTimes = eventTimes(:);

% Define search windows per event
winStart = eventTimes + mapDur(1);
winEnd = eventTimes + mapDur(2);

% Use discretize to find indices in the sorted spike vector
% This acts as a vectorized binary search
searchEdges = [-inf; spikeTimes; inf];

% Find first spike AFTER window start
startIdx = discretize(winStart, searchEdges);

% Find last spike BEFORE window end
% histc/discretize returns bin index.
% We use histc here for specific edge behavior matching original code logic
[~, endBin] = histc(winEnd, searchEdges);
endIdx = endBin - 1;

% Identify events that actually contain spikes
validWin = find(endIdx >= startIdx);

if isempty(validWin)
    return;
end

% Vectorize the extraction of relative times
% 1. How many spikes per valid event?
spkInEvt = endIdx(validWin) - startIdx(validWin) + 1;

% 2. Create an event ID vector matching the spikes
evtIds = repelem(validWin, spkInEvt);

% 3. Extract the spikes
% We create a list of indices to grab from spikeTimes
totalSpk = sum(spkInEvt);
spkIndices = zeros(totalSpk, 1);

% We need to fill spkIndices.
% A cumsum approach allows us to do this without a loop if desired,
% but a tight loop over valid events is often fast enough here.
% However, to match your requested speed, we loop only over valid windows.
curr = 1;
for k = 1:length(validWin)
    idx = validWin(k);
    n = spkInEvt(k);
    spkIndices(curr : curr+n-1) = startIdx(idx) : endIdx(idx);
    curr = curr + n;
end

spkVals = spikeTimes(spkIndices);

% 4. Calculate Relative Times and Bin
relTimes = spkVals - eventTimes(evtIds);
bins = discretize(relTimes, edges);

% 5. Accumulate into Map
% Remove any NaNs (spikes exactly on edge cases)
validBin = ~isnan(bins);
finalEvtIds = evtIds(validBin);
finalBins = bins(validBin);

if ~isempty(finalEvtIds)
    syncMap = accumarray([finalEvtIds, finalBins], 1, [nEvents, nBinsMap]);
end
end