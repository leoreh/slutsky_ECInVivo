function hFig = ripp_gui2(basepath, varargin)
% RIPP_GUI2  Interactive SWR viewer: raw LFP + RS firing-order raster.
%
%   hFig = RIPP_GUI2(basepath, varargin)
%
%   SUMMARY:
%       Loads SWR detection results and spike data from disk, then
%       displays a 2-panel interactive view of each event:
%           1. Raw LFP (ripple channel)
%           2. Spike raster (RS units only), sorted on every navigation
%              step by order of first firing in the display window.
%              Earliest-firing unit is at the bottom; latest at the top.
%              Units that do not fire appear above all firing units.
%              Burst spikes are overlaid in red.
%
%   CONTROLS:
%       [<] [>] Buttons or Left/Right Arrow Keys to navigate events.
%       'a' / 'd' Keys also navigate.
%       Edit box to jump directly to a numbered event.
%
%   INPUTS:
%       basepath  - (Char) Session directory. Must contain:
%                       <basename>.session.mat
%                       <basename>.ripp.mat
%                       <basename>.spikes.mat
%                       <basename>.units.mat
%                       <basename>.brst.mat
%                       <basename>.lfp  (binary)
%       varargin  - Parameter/Value pairs:
%           'winPlot'  - (Num) Window half-width around peak [s]. {0.5}
%
%   OUTPUTS:
%       hFig      - (Handle) Figure handle.
%
%   DEPENDENCIES:
%       binary_load, plot_raster.
%
%   HISTORY:
%       Created: 24 Mar 2026

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p,  'basepath', @ischar);
addParameter(p, 'winPlot',  0.5, @isnumeric);

parse(p, basepath, varargin{:});
winPlot = p.Results.winPlot;


%% ========================================================================
%  DATA LOADING
%  ========================================================================

[~, basename] = fileparts(basepath);

% --- Session Metadata ---
s = load(fullfile(basepath, [basename, '.session.mat']), 'session');
session = s.session;

fs     = session.extracellular.srLfp;        % LFP sampling rate [Hz]
nChans = session.extracellular.nChannels;    % Total channels in binary

% Ripple channel (used in binary_load)
if isfield(session, 'channelTags') && isfield(session.channelTags, 'Ripple')
    rippCh = session.channelTags.Ripple;
else
    rippCh = 1;
    warning('[RIPP_GUI2]: Ripple channel not found in session. Defaulting to ch 1.');
end

% Voltage conversion factor (Intan vs TDT)
if round(session.extracellular.sr) == 24414
    bit2uv = 1;       % TDT / Tucker-Davis
else
    bit2uv = 0.195;   % Intan
end

% Load Session & Data
vars = {'ripp', 'spikes', 'brst', 'units'};
v = basepaths2vars('basepaths', {basepath}, 'vars', vars);
ripp = v.ripp;
spikes = v.spikes;
units = v.units;
brst = v.brst;


%% ========================================================================
%  PRE-PROCESSING
%  ========================================================================

% --- RS Unit Spike Times ---
% spikes.times is a cell array of spike times in seconds, indexed by unit.
% brst.spktimes is structured identically (burst spikes only).
idxRS     = units.type == 'RS';
spktimes  = spikes.times(:);       % ensure column cell array
spktimes  = spktimes(idxRS);
brstTimes = brst.spktimes(:);      % ensure column cell array
brstTimes = brstTimes(idxRS);

% --- Raw LFP (ripple channel, full recording) ---
% Load from time 0 so timestamps align with absolute ripp.times.
fname = fullfile(basepath, [basename, '.lfp']);
lfp   = double(binary_load(fname, 'duration', Inf, 'fs', fs, ...
    'nCh', nChans, 'start', 0, 'ch', rippCh, 'downsample', 1, ...
    'bit2uv', bit2uv));

if size(lfp, 2) > 1
    lfp = mean(lfp, 2);     % average if multiple channels returned
end
lfp = lfp(:);               % ensure column vector

% Absolute time axis (matches ripp.times which are stored in absolute seconds)
timestamps = (0 : length(lfp) - 1)' / fs;

% Global y-limits for the LFP panel (robust to outliers)
yLimRaw = prctile(lfp, [0.1, 99.9]);


%% ========================================================================
%  GUI SETUP
%  ========================================================================

hFig = figure('Name', [basename, ' - Ripple GUI2'], 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], ...
    'KeyPressFcn', @onKeyPress);

% ---- Pack all GUI state into UserData ----
guiData            = struct();
guiData.rippTimes  = ripp.times;       % [N x 2] start / end times (s)
guiData.peakTimes  = ripp.peakTime;    % [N x 1] peak times (s)
guiData.lfp        = lfp;
guiData.timestamps = timestamps;
guiData.fs         = fs;
guiData.winPlot    = winPlot;
guiData.yLimRaw    = yLimRaw;
guiData.spktimes   = spktimes;         % RS unit spike times {nRS x 1}
guiData.brstTimes  = brstTimes;        % RS burst spike times {nRS x 1}
guiData.currIdx    = 1;
guiData.nEvents    = size(ripp.times, 1);
guiData.basename   = basename;

% ---- Layout: top = navigation bar, bottom = data panels ----
guiData.grid              = uigridlayout(hFig, [2, 1]);
guiData.grid.RowHeight    = {40, '1x'};
guiData.grid.ColumnWidth  = {'1x'};

% Navigation Panel (top)
hPanNav            = uipanel(guiData.grid);
hPanNav.BorderType = 'none';

hFlow = uiflowcontainer('v0', 'Parent', hPanNav, ...
    'FlowDirection', 'LeftToRight', 'Margin', 5);

% Prev button
hBtnPrev = uicontrol('Parent', hFlow, 'Style', 'pushbutton', 'String', '<', ...
    'Callback', @(s, e) navStep(-1));
set(hBtnPrev, 'WidthLimits', [30 30]);

% Event index edit box
guiData.hEditIdx = uicontrol('Parent', hFlow, 'Style', 'edit', ...
    'String', '1', 'Callback', @onEditIdx, 'BackgroundColor', 'w');
set(guiData.hEditIdx, 'WidthLimits', [50 60]);

% Total event count label
guiData.hTxtTotal = uicontrol('Parent', hFlow, 'Style', 'text', ...
    'String', sprintf('/ %d', guiData.nEvents), 'HorizontalAlignment', 'left');
set(guiData.hTxtTotal, 'WidthLimits', [60 80]);

% Next button
hBtnNext = uicontrol('Parent', hFlow, 'Style', 'pushbutton', 'String', '>', ...
    'Callback', @(s, e) navStep(1));
set(hBtnNext, 'WidthLimits', [30 30]);

% Spacer
hSpacer = uicontrol('Parent', hFlow, 'Style', 'text', 'String', '');
set(hSpacer, 'WidthLimits', [20 20]);

% Peak time info label
guiData.hTxtTime = uicontrol('Parent', hFlow, 'Style', 'text', ...
    'String', 't = 0.000 s', 'HorizontalAlignment', 'left');

% Plot Panel (bottom): 2 tiles stacked vertically
hPanPlot            = uipanel(guiData.grid);
hPanPlot.BorderType = 'none';
guiData.tiled       = tiledlayout(hPanPlot, 2, 1, ...
    'TileSpacing', 'tight', 'Padding', 'compact');

guiData.hAxRaw = nexttile(guiData.tiled);   % Raw LFP
guiData.hAxSpk = nexttile(guiData.tiled);   % Raster

linkaxes([guiData.hAxRaw, guiData.hAxSpk], 'x');

% ---- Store State & Draw ----
hFig.UserData = guiData;

if guiData.nEvents > 0
    updatePlot(hFig);
else
    title(guiData.hAxRaw, 'No Ripples Detected');
end


%% ====================================================================
%  CALLBACKS
%  ====================================================================

    % -----------------------------------------------------------------
    % Key Press: arrow keys or a/d to navigate
    % -----------------------------------------------------------------
    function onKeyPress(~, event)
        switch event.Key
            case {'rightarrow', 'd'},  navStep(1);
            case {'leftarrow',  'a'},  navStep(-1);
        end
    end

    % -----------------------------------------------------------------
    % Navigate: advance or retreat by one event
    % -----------------------------------------------------------------
    function navStep(step)
        data   = hFig.UserData;
        newIdx = data.currIdx + step;
        if newIdx >= 1 && newIdx <= data.nEvents
            data.currIdx  = newIdx;
            hFig.UserData = data;
            updatePlot(hFig);
        end
    end

    % -----------------------------------------------------------------
    % Edit Box: jump to a specific event by index
    % -----------------------------------------------------------------
    function onEditIdx(src, ~)
        data = hFig.UserData;
        val  = str2double(src.String);
        if isnan(val) || val < 1 || val > data.nEvents || val ~= round(val)
            src.String = num2str(data.currIdx);   % restore valid value
            return;
        end
        data.currIdx  = val;
        hFig.UserData = data;
        updatePlot(hFig);
    end

    % -----------------------------------------------------------------
    % Update Plot: redraw both panels for the current event index
    % -----------------------------------------------------------------
    function updatePlot(figH)
        data = figH.UserData;
        idx  = data.currIdx;

        % Safety clamp
        idx          = max(1, min(idx, data.nEvents));
        data.currIdx = idx;

        % Sync controls
        data.hEditIdx.String  = num2str(idx);
        data.hTxtTotal.String = sprintf('/ %d', data.nEvents);

        % Early exit if empty session
        if data.nEvents == 0
            cla(data.hAxRaw);  cla(data.hAxSpk);
            title(data.hAxRaw, 'No Ripples');
            return;
        end

        % Event markers
        pkT = data.peakTimes(idx);
        stT = data.rippTimes(idx, 1);
        enT = data.rippTimes(idx, 2);
        data.hTxtTime.String = sprintf('Peak t = %.3f s', pkT);

        % Display window boundaries
        winStart = pkT - data.winPlot;
        winEnd   = pkT + data.winPlot;

        % Sample index range (1-indexed)
        sampleStart = max(1, round(winStart * data.fs));
        sampleEnd   = min(length(data.lfp), round(winEnd * data.fs));
        rangeIdx    = sampleStart : sampleEnd;
        tVec        = data.timestamps(rangeIdx);

        if isempty(tVec), return; end

        % --- Panel 1: Raw LFP ---
        cla(data.hAxRaw);
        plot(data.hAxRaw, tVec, data.lfp(rangeIdx), 'k');
        xline(data.hAxRaw, [stT, pkT, enT], '--b');
        ylabel(data.hAxRaw, 'LFP (\muV)');
        grid(data.hAxRaw, 'on');
        xlim(data.hAxRaw, [winStart, winEnd]);
        ylim(data.hAxRaw, data.yLimRaw);

        % --- Panel 2: Raster (firing-order sorted) ---
        cla(data.hAxSpk);

        % Prune spike times to the display window (efficiency)
        spkWin  = cellfun(@(s) s(s >= winStart & s <= winEnd), ...
            data.spktimes,  'UniformOutput', false);
        brstWin = cellfun(@(s) s(s >= winStart & s <= winEnd), ...
            data.brstTimes, 'UniformOutput', false);

        % Reorder units: earliest first spike at bottom (index 1),
        % latest at top, non-firing units above all others
        sortOrder  = firingOrderSort(spkWin);
        spkSorted  = spkWin(sortOrder);
        brstSorted = brstWin(sortOrder);

        % All spikes — black
        plot_raster(spkSorted,  'hAx', data.hAxSpk, ...
            'xLim', [winStart, winEnd], 'flgLbls', false, 'clr', [0 0 0]);

        % Burst spikes overlaid — red (skip call if none fire in window)
        if any(~cellfun(@isempty, brstSorted))
            plot_raster(brstSorted, 'hAx', data.hAxSpk, ...
                'xLim', [winStart, winEnd], 'flgLbls', false, 'clr', [1 0 0]);
        end

        xline(data.hAxSpk, [stT, pkT, enT], '--b');
        ylabel(data.hAxSpk, 'Units (RS)');
        xlabel(data.hAxSpk, 'Time (s)');

        figH.UserData = data;
    end

    % -----------------------------------------------------------------
    % Firing Order Sort: returns a permutation vector such that position 1
    % corresponds to the unit with the earliest first spike in spkWin.
    % Units with no spikes in the window are placed at the end (top).
    % -----------------------------------------------------------------
    function sortOrder = firingOrderSort(spkWin)
        nUnits   = length(spkWin);
        firstSpk = inf(nUnits, 1);
        for k = 1 : nUnits
            if ~isempty(spkWin{k})
                firstSpk(k) = min(spkWin{k});
            end
        end
        [~, sortOrder] = sort(firstSpk, 'ascend');
    end

end     % EOF
