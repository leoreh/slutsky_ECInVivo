function mancur_stg(varargin)

% manually curate significant synapses. see calc_stgs.m
%
% INPUT:
%   stg         struct containing stg data
%   basepath    string specifying the file directory {pwd}
%
% DEPENDENCIES
%   plot_ccg
%
% 10 may 24 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'stg', []);
addParameter(p, 'basepath', pwd, @ischar);

parse(p, varargin{:});
stg = p.Results.stg;
basepath = p.Results.basepath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize global data structure to contain everything
global g

% files
[~, basename] = fileparts(basepath);
g.stgfile = fullfile(basepath, [basename, '.stg.mat']);

% load
vars = ["swv_metrics"; "units"; "session"];
v = basepaths2vars('basepaths', {basepath}, 'vars', vars);
if isfield(v, 'swv_metrics'), [v.swv] = v.swv_metrics; v = rmfield(v, 'swv_metrics'); end
if isempty(stg)
    load(g.stgfile)
end

% organize
g.fs = v.session.extracellular.sr;
g.wv = v.swv.wv;
g.wv_std = v.swv.wv_std;
g.stg = stg;
g.unitType = [1 2] * v.units.clean;
g.unitType(g.unitType == 0) = 3;

% current synapse params
g.nSyn = size(stg.synIdx, 1);
g.currentIdx = 1;
g.uIdx = g.stg.synIdx(g.currentIdx, [1, 2]);
g.synType = g.stg.synIdx(g.currentIdx, 3);

% initialize user input. 1 for accept, -1 for reject, 0 for undecided
g.accepted = zeros(g.nSyn, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% graphic params
g.clr = 'brk';
g.wv_xVal = [1 : size(g.wv, 2)] / g.fs * 1000;

% figure
setMatlabGraphics(true)
fh = figure('Name', 'STG Reviewer', 'NumberTitle', 'off',...
    'Position', [100, 100, 900, 500], 'CloseRequestFcn', @closeRequest);
set(fh, 'WindowState', 'maximized');
set(fh, 'DefaultAxesFontSize', 16);
g.th = tiledlayout(2, 4);
g.th.TileSpacing = 'tight';


% manually create tilelayout
tlayout = [2, 4];
pWidth = 0.15;
pHeight = 0.4;
leftMarg = 0.05;
topMarg = 0.5;
vertSpac = 0.04;
horzSpac = 0.03;

% create axes
cnt = 1;
for irow = 1 : tlayout(1)
    for icol = 1 : tlayout(2)
        g.axh(cnt) = axes('Parent', fh, 'Position', ...
            [leftMarg + (icol - 1) * (pWidth + horzSpac), ...
            topMarg - (irow - 1) * (pHeight + vertSpac), ...
            pWidth, pHeight]);
        axis(g.axh(cnt), 'square'); % Make axes square
        cnt = cnt + 1;
    end
end

% create 2x2 array of buttons
btnTxt = ["Accept"; "Reject"; "Next"; "Prev"; "Save"; "Nothing"]
btnWidth = 0.05;
btnHeight = 0.05;
leftMarg = 0.8;
topMarg = 0.3;
vertSpac = 0.05;
horzSpac = 0.03;
cnt = 1;
for irow = 1 : 3
    for icol = 1 : 2
        btnPos = [leftMarg + (icol - 1) * (btnWidth + horzSpac), ...
            topMarg - (irow - 1) * (btnHeight + vertSpac), ...
            btnWidth, btnHeight];

        uicontrol('Parent', fh, 'Style', 'pushbutton', 'String', btnTxt{cnt},...
            'Units', 'normalized', 'Position', btnPos,...
            'Callback', @(src, event) set_cbk(src, event, cnt));

        cnt = cnt + 1;
    end
end

% create drop down menu of synaptic pairs
synStr = arrayfun(@(x) sprintf('%d-%d %d', g.stg.synIdx(x, 1), g.stg.synIdx(x, 2), g.stg.synIdx(x, 3)),...
    1:size(g.stg.synIdx, 1), 'UniformOutput', false);
uicontrol('Parent', fh, 'Style', 'popupmenu', 'String', synStr, ...
    'Units', 'normalized', 'Position', [0.8, 0.4, btnWidth * 2 + horzSpac, 0.05], ...
    'Callback', @dropdown_cbk);

% Load and display the first synapse
updateDisplay();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% callbacks for buttons
    function set_cbk(src, event, cnt)
        switch cnt
            case 1      % accept
                g.accepted(g.currentIdx) = 1;
                nextSynapse();

            case 2      % reject
                g.accepted(g.currentIdx) = -1;
                nextSynapse();

            case 3      % next
                nextSynapse();

            case 4      % previous
                if g.currentIdx > 1
                    g.currentIdx = g.currentIdx - 1;
                    g.uIdx = g.stg.synIdx(g.currentIdx, [1, 2]);
                    g.synType = g.stg.synIdx(g.currentIdx, 3);
                end
                updateDisplay();

            case 5      % save
                saveStg()
        end
    end

% handle dropdown selection
    function dropdown_cbk(src, ~)
        g.currentIdx = src.Value;
        g.uIdx = g.stg.synIdx(g.currentIdx, [1, 2]);
        g.synType = g.stg.synIdx(g.currentIdx, 3);
        updateDisplay();
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gui operation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Move to the next synapse if possible
    function nextSynapse()
        if g.currentIdx < g.nSyn
            g.currentIdx = g.currentIdx + 1;
            g.uIdx = g.stg.synIdx(g.currentIdx, [1, 2]);
            g.synType = g.stg.synIdx(g.currentIdx, 3);

            updateDisplay();
        else
            closeRequest();
        end
    end

% save the manual curation
    function saveStg()
        stg = g.stg;
        stg.accepted = g.accepted;
        save(g.stgfile, 'stg')
    end

% called when the user tries to close the figure
    function closeRequest(src, event)
        selection = questdlg('Do you want to save changes before closing?', ...
            'Close Request Function', 'Yes', 'No', 'Cancel', 'Yes');

        switch selection
            case 'Yes'
                saveStg();
                delete(gcf);

            case 'No'
                delete(gcf);

            case 'Cancel'
                return
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update the display with the current synapse
    function updateDisplay()

        % extract data from struct
        dccc = g.stg.cc.dccc(:, g.uIdx(1), g.uIdx(2));
        pred = g.stg.info.pred(:, g.uIdx(1), g.uIdx(2));
        eBins = g.stg.eBins(:, g.uIdx(1), g.uIdx(2));
        iBins = g.stg.iBins(:, g.uIdx(1), g.uIdx(2));
        roi_binIdx = g.stg.info.roi_binIdx;

        % update title
        title(g.th, sprintf('%s\nSynapse %d of %d',...
            g.stg.info.basename, g.currentIdx, g.nSyn),...
            'interpreter', 'none', 'FontSize', 20);

        % acg of presynaptic unit
        axes(g.axh(1)); cla; hold on
        plot_ccg(g.stg.cc.cc150(:, g.uIdx(1), g.uIdx(1)), g.stg.cc.cc150bins,...
            'clr', g.clr(g.unitType(g.uIdx(1))),...
            'pred', [], 'sigbins1', [], 'sigbins2', []);
        title(gca, sprintf('Unit #%d (Presynaptic)', g.uIdx(1)))

        % acg of postsynaptic unit
        axes(g.axh(4)); cla; hold on
        plot_ccg(g.stg.cc.cc150(:, g.uIdx(2), g.uIdx(2)), g.stg.cc.cc150bins,...
            'clr', g.clr(g.unitType(g.uIdx(2))),...
            'pred', [], 'sigbins1', [], 'sigbins2', []);
        title(gca, sprintf('Unit #%d (Presynaptic)', g.uIdx(2)))

        % ccg 50
        axes(g.axh(2)); cla; hold on
        plot_ccg(g.stg.cc.cc50(:, g.uIdx(1), g.uIdx(2)), g.stg.cc.cc50bins, 'clr', 'k',...
            'pred', [], 'sigbins1', roi_binIdx(eBins), 'sigbins2', roi_binIdx(iBins));

        % ccg 150
        axes(g.axh(6)); cla; hold on
        plot_ccg(g.stg.cc.cc150(:, g.uIdx(1), g.uIdx(2)), g.stg.cc.cc150bins, 'clr', 'k',...
            'pred', [], 'sigbins1', [], 'sigbins2', []);

        % ccg 20
        axes(g.axh(7)); cla; hold on
        plot_ccg(g.stg.cc.cc20(:, g.uIdx(1), g.uIdx(2)), g.stg.cc.cc20bins, 'clr', 'k',...
            'pred', [], 'sigbins1', [], 'sigbins2', []);
        xticks([stg.cc.cc20bins(1) : 2 : stg.cc.cc20bins(end)])
        xtickangle(0)

        % deconvoluted cc
        axes(g.axh(3)); cla; hold on
        plot_ccg(dccc, g.stg.cc.cc50bins, 'clr', 'k', 'pred', pred,...
            'sigbins1', roi_binIdx(eBins), 'sigbins2', roi_binIdx(iBins))
        if g.synType == -1
            title(gca, sprintf('iSTG = %.4f', g.stg.iStg(g.uIdx(1), g.uIdx(2))))
        elseif g.synType == 1
            title(gca, sprintf('eSTG = %.4f', g.stg.eStg(g.uIdx(1), g.uIdx(2))))
        end

        % plot waveform of presynaptic unit
        axes(g.axh(5)); cla; hold on
        plot_wv(g.uIdx(1))

        % plot waveform of postsynaptic unit
        axes(g.axh(8)); cla; hold on
        plot_wv(g.uIdx(2))

    end

% plot waveform
    function plot_wv(uIdx)
        wv = g.wv(uIdx, :);
        wv_std = g.wv_std(uIdx, :);
        clr = g.clr(g.unitType(uIdx));
        plot(g.wv_xVal, wv, clr, 'LineWidth', 2)
        patch([g.wv_xVal, flip(g.wv_xVal)], [wv + wv_std, flip(wv - wv_std)],...
            clr, 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
        xlabel('Time [ms]')
        ylabel('Voltage [mV]')
    end


end