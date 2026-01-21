function hFig = ripp_gui(rippTimes, peakTimes, rippSig, spktimes, fs, thr, varargin)
% RIPP_GUI Interactive visualization of detected ripples.
%
% SUMMARY:
%   Displays a 5-panel view of each ripple event:
%   1. Raw LFP
%   2. Filtered LFP (100-300Hz)
%   3. Detection Signal (z-scored)
%   4. EMG
%   5. Spike Raster
%
%   Navigation: Left/Right Arrows or 'a'/'d' keys.
%   Controls: interactive Event Index box, Total Count, Time Info.
%
% INPUTS:
%   rippTimes   - (N x 2) Ripple start/end times [s].
%   peakTimes   - (N x 1) Ripple peak times [s].
%   lfp         - (Vec) Raw LFP signal.
%   spktimes    - (Cell) Spike times for raster.
%   fs          - (Num) LFP sampling rate [Hz].
%   emg         - (Vec) EMG signal.
%   thr         - (Vec) Thresholds used [start, peak, cont...].
%   varargin    - Optional: 'basepath', 'detectAlt' (default 3), 'winPlot' (default 0.5s).
%
% OUTPUT:
%   hFig        - Figure handle.
%
% DEPENDENCIES: filterLFP, plot_raster, setMatlabGraphics.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'rippTimes', @isnumeric);
addRequired(p, 'peakTimes', @isnumeric);
addRequired(p, 'rippSig', @isstruct);
addRequired(p, 'spktimes', @iscell);
addRequired(p, 'fs', @isnumeric);
addRequired(p, 'thr', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'detectAlt', 3, @isnumeric);
addParameter(p, 'winPlot', 0.5, @isnumeric); % +/- window in seconds

parse(p, rippTimes, peakTimes, rippSig, spktimes, fs, thr, varargin{:});

basepath = p.Results.basepath;
detectAlt = p.Results.detectAlt;
winPlot = p.Results.winPlot;

%% ========================================================================
%  PRE-PROCESSING
%  ========================================================================

fprintf('Initializing Ripple GUI... Pre-processing signals...\n');

% Extract from rippSig
lfp = rippSig.lfp;
sigFilt = rippSig.filt;
sigDetect = rippSig.z;
emg = rippSig.emg;

% Time Vector
timestamps = (0:length(lfp)-1)' / fs;


%% ========================================================================
%  GUI SETUP
%  ========================================================================

[~, basename] = fileparts(basepath);
hFig = figure('Name', [basename ' - Ripple GUI'], 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], ...
    'KeyPressFcn', @onKeyPress);

% Store Data
guiData = struct();
guiData.rippTimes = rippTimes;
guiData.peakTimes = peakTimes;
guiData.lfp = lfp;
guiData.sigFilt = sigFilt;
guiData.sigDetect = sigDetect;
guiData.emg = emg;
guiData.spktimes = spktimes;
guiData.timestamps = timestamps;
guiData.fs = fs;
guiData.thr = thr;
guiData.winPlot = winPlot;
guiData.currIdx = 1;
guiData.nEvents = size(rippTimes, 1);
guiData.basename = basename;

% Layout
% Main Grid: Top Panel (Controls) + Bottom Panel (Plots)
guiData.grid = uigridlayout(hFig, [2, 1]);
guiData.grid.RowHeight = {40, '1x'};
guiData.grid.ColumnWidth = {'1x'};

hPanNav = uipanel(guiData.grid);
hPanNav.BorderType = 'none';

% Controls in hPanNav
% We want: [Prev] [EditBox] [/ Total] [Next] ... [Time Info]
hFlow = uiflowcontainer('v0', 'Parent', hPanNav, 'FlowDirection', 'LeftToRight', 'Margin', 5);

% Prev Button
% Prev Button
hBtnPrev = uicontrol('Parent', hFlow, 'Style', 'pushbutton', 'String', '<', ...
    'Callback', @(s,e) navStep(-1));
set(hBtnPrev, 'WidthLimits', [30 30]);

% Edit Box (Event Index)
% Edit Box (Event Index)
guiData.hEditIdx = uicontrol('Parent', hFlow, 'Style', 'edit', ...
    'String', '1', 'Callback', @onEditIdx, ...
    'BackgroundColor', 'w');
set(guiData.hEditIdx, 'WidthLimits', [50 60]);

% Total Count Text
guiData.hTxtTotal = uicontrol('Parent', hFlow, 'Style', 'text', ...
    'String', sprintf('/ %d', guiData.nEvents), ...
    'HorizontalAlignment', 'left');
set(guiData.hTxtTotal, 'WidthLimits', [60 80]);

% Next Button
% Next Button
hBtnNext = uicontrol('Parent', hFlow, 'Style', 'pushbutton', 'String', '>', ...
    'Callback', @(s,e) navStep(1));
set(hBtnNext, 'WidthLimits', [30 30]);

% Spacer
% Spacer
hSpacer = uicontrol('Parent', hFlow, 'Style', 'text', 'String', '');
set(hSpacer, 'WidthLimits', [20 20]);

% Time Info Text
guiData.hTxtTime = uicontrol('Parent', hFlow, 'Style', 'text', ...
    'String', 't = 0.000 s', 'HorizontalAlignment', 'left');

% Bottom Panel (Plots)
hPanPlot = uipanel(guiData.grid);
hPanPlot.BorderType = 'none';
guiData.tiled = tiledlayout(hPanPlot, 5, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

% Axes Handles (Initialize)
guiData.hAxRaw = nexttile(guiData.tiled);
guiData.hAxFilt = nexttile(guiData.tiled);
guiData.hAxDet = nexttile(guiData.tiled);
guiData.hAxEmg = nexttile(guiData.tiled);
guiData.hAxSpk = nexttile(guiData.tiled);

% Link Axes
linkaxes([guiData.hAxRaw, guiData.hAxFilt, guiData.hAxDet, guiData.hAxEmg, guiData.hAxSpk], 'x');

% Store Data
hFig.UserData = guiData;

% Initial Plot
if guiData.nEvents > 0
    updatePlot(hFig);
else
    title(guiData.hAxRaw, 'No Ripples Detected');
end


%% ====================================================================
%  CALLBACKS
%  ====================================================================

    function onKeyPress(src, event)
        data = src.UserData;
        switch event.Key
            case {'rightarrow', 'd'}
                navStep(1);
            case {'leftarrow', 'a'}
                navStep(-1);
        end
    end

    function navStep(step)
        data = hFig.UserData; % Reload
        newIdx = data.currIdx + step;
        if newIdx >= 1 && newIdx <= data.nEvents
            data.currIdx = newIdx;
            hFig.UserData = data;
            updatePlot(hFig);
        end
    end

    function onEditIdx(src, ~)
        data = hFig.UserData; % Reload
        val = str2double(src.String);
        if isnan(val) || val < 1 || val > data.nEvents || val ~= round(val)
            % Restore previous val
            src.String = num2str(data.currIdx);
            return;
        end
        data.currIdx = val;
        hFig.UserData = data;
        updatePlot(hFig);
    end

    function updatePlot(figH)
        data = figH.UserData;
        idx = data.currIdx;

        % Check bounds (safety for external updates)
        if idx > data.nEvents, idx = data.nEvents; end
        if idx < 1, idx = 1; end
        data.currIdx = idx;

        % Update Controls
        data.hEditIdx.String = num2str(idx);
        data.hTxtTotal.String = sprintf('/ %d', data.nEvents);

        % Event Info
        if data.nEvents > 0
            pkT = data.peakTimes(idx);
            stT = data.rippTimes(idx, 1);
            enT = data.rippTimes(idx, 2);
            data.hTxtTime.String = sprintf('Peak t = %.3f s', pkT);
        else
            pkT = 0; stT=0; enT=0;
            data.hTxtTime.String = 'No Data';
        end

        % If no events, just clear and return
        if data.nEvents == 0
            cla(data.hAxRaw); cla(data.hAxFilt); cla(data.hAxDet); cla(data.hAxEmg); cla(data.hAxSpk);
            title(data.hAxRaw, 'No Ripples');
            return;
        end

        winStart = pkT - data.winPlot;
        winEnd = pkT + data.winPlot;

        % Indices for window
        % Using simple logical indexing or find for speed (assuming vector timestamps)
        % For speed, calculate index range directly
        sampleStart = max(1, round(winStart * data.fs));
        sampleEnd = min(length(data.lfp), round(winEnd * data.fs));
        rangeIdx = sampleStart:sampleEnd;
        tVec = data.timestamps(rangeIdx);

        if isempty(tVec), return; end

        % --- 1. RAW LFP ---
        cla(data.hAxRaw);
        plot(data.hAxRaw, tVec, data.lfp(rangeIdx), 'k');
        xline(data.hAxRaw, [stT, pkT, enT], '--b');
        title(data.hAxRaw, ''); % Title info moved to top panel
        ylabel(data.hAxRaw, 'Raw');
        grid(data.hAxRaw, 'on'); xlim(data.hAxRaw, [winStart, winEnd]);

        % --- 2. FILTERED LFP ---
        cla(data.hAxFilt);
        plot(data.hAxFilt, tVec, data.sigFilt(rangeIdx), 'k');
        xline(data.hAxFilt, [stT, pkT, enT], '--b');
        ylabel(data.hAxFilt, 'Filt');
        grid(data.hAxFilt, 'on'); xlim(data.hAxFilt, [winStart, winEnd]);

        % --- 3. DETECTION SIGNAL ---
        cla(data.hAxDet);
        plot(data.hAxDet, tVec, data.sigDetect(rangeIdx), 'k', 'LineWidth', 1.5);
        xline(data.hAxDet, [stT, pkT, enT], '--b');
        yline(data.hAxDet, data.thr(1), '--g', 'Thr1 (Start)');
        yline(data.hAxDet, data.thr(2), '--r', 'Thr2 (Peak)');
        yline(data.hAxDet, data.thr(3), '--m', 'Thr3 (Cont)');
        ylabel(data.hAxDet, 'Z-Score');
        grid(data.hAxDet, 'on'); xlim(data.hAxDet, [winStart, winEnd]);

        % --- 4. EMG ---
        cla(data.hAxEmg);
        plot(data.hAxEmg, tVec, data.emg(rangeIdx), 'k');
        xline(data.hAxEmg, [stT, pkT, enT], '--b');
        ylabel(data.hAxEmg, 'EMG');
        grid(data.hAxEmg, 'on'); xlim(data.hAxEmg, [winStart, winEnd]);

        % --- 5. RASTER ---
        cla(data.hAxSpk);
        % Prepare spikes in window
        % plot_raster takes cell array. We pass all and let it handle/zoom or filter?
        % Passing everything is slow if we have many spikes.
        % Prune spikes to window for efficiency
        spkWin = cell(size(data.spktimes));
        for k = 1:length(data.spktimes)
            s = data.spktimes{k};
            spkWin{k} = s(s >= winStart & s <= winEnd);
        end

        plot_raster(spkWin, 'hAx', data.hAxSpk, 'xLim', [winStart, winEnd], 'flgLbls', false);
        xline(data.hAxSpk, [stT, pkT, enT], '--b');
        ylabel(data.hAxSpk, 'Units');
        xlabel(data.hAxSpk, 'Time (s)');

    end

end
