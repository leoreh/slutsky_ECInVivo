function ssEmg = mcu_states_EmgMouse(mname, graphics)


% loads the sessions of a mouse. classfies states by EMG (as_emg).
% organizes the output in a struct and plots

if nargin == 1
    graphics = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state bouts by emg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths = mcu_basepaths(mname);
nfiles = length(basepaths);
mpath = fileparts(basepaths{1});

% select session
for ifile = 1 : nfiles
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)


    ssEmg = as_emg('basepath', basepath);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ssEmg structs
vars = ["sleep_statesEmg"];
basepaths = mcu_basepaths(mname);
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
ssEmg = catfields([v(:).ssEmg], 'addim', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics

    % organize data
    boutTimes = squeeze(ssEmg.bouts.times);

    % fig params
    cfg = as_loadConfig('flgEmg', true);
    nstates = 2;
    clr = cfg.colors;
    snames = cfg.names;

    % -------------------------------------------------------------------------
    % plot hypnogram for all sessions

    % open figure
    fh = figure;
    setMatlabGraphics(true)
    set(fh, 'WindowState','maximized');
    tlayout = [nfiles, 1];
    th = tiledlayout(tlayout(1), tlayout(2));
    th.TileSpacing = 'tight';
    th.Padding = 'none';
    title(th, mname, 'interpreter', 'none')

    % hypnogram
    for ifile = 1 : nfiles
        axh = nexttile(th, ifile, [1, 1]); hold on; cla
        plot_hypnogram('boutTimes', boutTimes(:, ifile), 'clr', clr, 'axh', axh,...
            'sstates', [1 : nstates])
        axis tight
        xval = 3600 : 3600 : length(v(ifile).ssEmg.labels);
        xticks(xval)
        xticklabels(string(xval / 3600))
        xlabel('Time (hr)')
    end

    figpath = fullfile(mpath, 'graphics', 'sleepState');
    mkdir(figpath)
    figname = fullfile(figpath, [mname, '_hypnogramEmg.png']);
    saveas(fh, figname)


    % -------------------------------------------------------------------------
    % plot bout states across sessions

    % organize data
    boutLen = squeeze(ssEmg.bouts.boutLen);
    nbouts = squeeze(ssEmg.bouts.nbouts);
    prctDur = squeeze(ssEmg.bouts.prctDur);

    % open figure
    fh = figure;
    setMatlabGraphics(true)
    set(fh, 'WindowState','maximized');
    tlayout = [3, 2];
    th = tiledlayout(tlayout(1), tlayout(2));
    th.TileSpacing = 'tight';
    th.Padding = 'none';
    title(th, mname, 'interpreter', 'none')

    % bout length
    tilebias = 0;
    for istate = 1 : 2

        axh = nexttile(th, tilebias + istate, [1, 1]); hold on; cla
        boutMat = cell2padmat(boutLen(istate, :), 2);
        plot_boxMean('axh', axh, 'dataMat', boutMat, 'plotType', 'bar',...
            'clr', clr{istate})
        axis tight
        ylabel('Bout Length (s)')
        xlabel('Session No.')
        title(snames{istate})
    end

    % number of bouts
    tilebias = 2;
    for istate = 1 : 2

        axh = nexttile(th, tilebias + istate, [1, 1]); hold on; cla
        plot_boxMean('axh', axh, 'dataMat', nbouts, 'plotType', 'bar',...
            'clr', clr{istate})
        axis tight
        ylabel('No. Bouts')
        xlabel('Session No.')
        title(snames{istate})
    end

    % state duration
    axh = nexttile(th, 5, [1, 1]); hold on; cla
    bh = bar(axh, prctDur' / 60 / 60, 'stacked');
    bh(1).FaceColor = clr{1};
    bh(2).FaceColor = clr{2};
    axis tight
    ylabel('State Duration (h)')
    xlabel('Session No.')

    % state duration (%)
    axh = nexttile(th, 6, [1, 1]); hold on; cla
    prctDur = prctDur ./ sum(prctDur) * 100;
    bh = bar(axh, prctDur', 'stacked');
    bh(1).FaceColor = clr{1};
    bh(2).FaceColor = clr{2};
    axis tight
    ylabel('State Duration (%)')
    xlabel('Session No.')

    figpath = fullfile(mpath, 'graphics', 'sleepState');
    mkdir(figpath)
    figname = fullfile(figpath, [mname, '_boutTimesEmg.png']);
    saveas(fh, figname)

end

end

% EOF