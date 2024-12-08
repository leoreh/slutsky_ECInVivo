
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare baseline and mk801
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% groups
queryStr = ["bsl"; "mk801"; "wsh"];
ngrp = length(queryStr);

clear drft2 drft
for igrp = 1 : ngrp

    % load data
    basepaths = mk801_sessions(queryStr{igrp});
    nfiles = length(basepaths);

    % go over files
    clear drft1
    for ifile = 1 : nfiles

        % file params
        basepath = basepaths{ifile};
        cd(basepath)
        [~, basename] = fileparts(basepath);

        drft1(ifile) = drift_file('basepath', basepath, 'graphics', false);

    end

    drft2(igrp) = catfields(drft1, 'addim');

end
drft = catfields(drft2, 'addim');

% dimensions of dgrp are (example m_corr):
% data x unit x state x file x grp

% -------------------------------------------------------------------------
% graphics

ulabel = ["RS"; "FS"];
slabel = ["Full Recordings"; "AW"; "NREM"];

setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [2, length(slabel)];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
set(fh, 'DefaultAxesFontSize', 16);

yLimit = [0, round(max(drft.drate, [], 'all') * 100) / 100];
tbias = 1;
for iunit = 1 : 2
    for istate = 1 : length(slabel)
        axh = nexttile(th, tbias, [1, 1]); cla; hold on
        dataMat = squeeze(drft.drate(iunit, istate, :, :));
        plot_boxMean('dataMat', dataMat, 'plotType', 'bar',...
            'allPnts', false, 'axh', axh)
        plot(axh, [1 : ngrp], dataMat)

        ylim(yLimit)
        ylabel(sprintf('%s Drift Rate [1 / h]', ulabel(iunit)))
        title(axh, slabel(istate))
        tbias = tbias + 1;
    end
end

% data for prism
istate = 1;
iunit = 1;
tmp = squeeze(drft.drate(iunit, istate, :, :));
prismData = reshape(tmp, 2, 10)


igrp = 2;
tmp = squeeze(drft.drate(:, [2, 3], :, igrp));
tmp = permute(tmp, [1, 3, 2]);
prismData = reshape(tmp, 2, 10)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% investigate drift during baseline as a function of PV size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% ogranize and load data

% files
basepaths = {...
    'F:\Data\lh122\lh122_230109_095412';...
    'F:\Data\lh123\lh123_221219_094508';...
    'F:\Data\lh126\lh126_230111_091208';...
    'F:\Data\lh129\lh129_230123_095540';...
    'F:\Data\lh130\lh130_230322_084541';...
    'F:\Data\lh96\lh96_220120_090157';...
    'F:\Data\lh107\lh107_220518_091200';...
    'F:\Data\lh142\lh142_231005_091832';...
    'F:\Data\lh100\lh100_220413_111004';...
    'F:\Data\lh122\lh122_221223_092656';...
    'F:\Data\lh111\lh111_220823_094417';...
    'F:\Data\lh87\lh87_210523_100607';
    };

mnames = cellfun(@(x) regexp(x, 'lh\d+', 'match'), basepaths, 'uni', false);
mnames = cellfun(@(x) x{1}, mnames, 'uni', false);

% general params
varsFile = ["fr"; "units"];
varsName = ["fr"; "units"];

% load data
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"]);
nfiles = length(basepaths);


% -------------------------------------------------------------------------
% analyze

% analysis params
winsize = 3600;
thrLin = 0;
thrFr = 0.005;
thrWin = 10;

% unit params
minUnits = 2;       % minimum no. of units for drift calculation

% go over files
clear drft3 
for ifile = 1 : nfiles

    % file params
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % time stamps
    tstamps = v(ifile).fr.tstamps;
    
    clear drft2
    for sunit = 1 : 2       % rs / fs

        % select units
        unitIdx = find(v(ifile).units.clean(sunit, :));
        fr_mat = v(ifile).fr.strd(unitIdx, :);
        nunits = length(unitIdx);

        % calc drift for increasing number of units
        clear drft1
        cnt = 1;
        for iunit = minUnits : nunits

            drft1(cnt) = drift_calc(fr_mat, tstamps, 'graphics', false,...
                'winsize', winsize, 'thrLin', thrLin, 'thrFr', thrFr,...
                'thrWin', thrWin, 'limUnit', iunit);

            cnt = cnt + 1;
        end

        drft2(sunit) = catfields(drft1, 'addim');
    end

    drft3(ifile) = catfields(drft2, 'addim');
end

drft = catfields(drft3, 'addim');

drate = squeeze(drft.drate);
unitVec = 2 : size(drate, 1) + minUnits - 1;

% -------------------------------------------------------------------------
% graphics: drift vs. pv size

ulabel = ["RS"; "FS"];

setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [2, 2];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
set(fh, 'DefaultAxesFontSize', 16);

tbias = 0;
for sunit = 1 : 2
    tbias = tbias + sunit;
    
    dataMat = squeeze(drate(:, sunit, :));
    
    axh = nexttile(th, tbias, [1, 1]); cla; hold on
    plot(unitVec, dataMat)
    xlim([unitVec(1), unitVec(end)])
    xlabel('No. Units in PV')
    ylabel(sprintf('%s Drift Rate [1 / h]', ulabel(sunit)))
    legend(mnames)

    axh = nexttile(th, tbias + 1, [1, 1]); cla; hold on
    plot_stdShade('dataMat', dataMat, 'xVal', unitVec, 'axh', axh,...
        'clr', [0 0 0], 'alpha', 0.5)
    xlim([unitVec(1), unitVec(end)])
    xlabel('No. Units in PV')
    ylabel(sprintf('%s Drift Rate [1 / h]', ulabel(sunit)))
end


% compare rs vs fs of selected files and nunits
% 123   2
% 126   3
% 130   5
% 87    12

% 10 units      
cntIdx = find(unitVec == 10);
dataMat = squeeze(drate(cntIdx, :, [2, 3, 5]));
fh = figure;
plot_boxMean('dataMat', dataMat', 'plotType', 'bar',...
    'allPnts', false)

pVal = stat_compare1D(dataMat', 'flgRep', false, 'flgParam', true)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare drift in AW vs NREM during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% ogranize and load data

% files
basepaths = {...
    'F:\Data\lh122\lh122_230109_095412';...
    'F:\Data\lh123\lh123_221219_094508';...
    'F:\Data\lh126\lh126_230111_091208';...
    'F:\Data\lh129\lh129_230123_095540';...
    'F:\Data\lh130\lh130_230322_084541';...
    'F:\Data\lh96\lh96_220120_090157';...
    'F:\Data\lh107\lh107_220518_091200';...
    'F:\Data\lh142\lh142_231005_091832';...
    'F:\Data\lh100\lh100_220413_111004';...
    'F:\Data\lh122\lh122_221223_092656';...
    'F:\Data\lh111\lh111_220823_094417';...
    'F:\Data\lh87\lh87_210523_100607';
    };

mnames = cellfun(@(x) regexp(x, 'lh\d+', 'match'), basepaths, 'uni', false);
mnames = cellfun(@(x) x{1}, mnames, 'uni', false);

% general params
varsFile = ["fr"; "units"];
varsName = ["fr"; "units"];

% load data
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"]);
nfiles = length(basepaths);

% -------------------------------------------------------------------------
% analyze drift 

clear drft1 drft
for ifile = 1 : nfiles

    % file params
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    drft1(ifile) = drift_file('basepath', basepath, 'graphics', false);

end
drft = catfields(drft1, 'addim');
iunit = 1;

% -------------------------------------------------------------------------
% graphics: aw vs nrem

cfg = as_loadConfig;
sclr = [{[0, 0, 0]}; cfg.colors([1, 4])];
snames = [{'FULL REC'}; cfg.names([1, 4])];

setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [1, 2];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
set(fh, 'DefaultAxesFontSize', 16);

axh = nexttile(th, 1, [1, 1]); cla; hold on
xval = [1 : size(drft.m_corr, 1)];
for istate = 1 : size(drft.drate, 2);
    lin_coef = squeeze(drft.lin_coef(:, iunit, istate, :));
    lin_mat = [lin_coef(2, :) + lin_coef(1, :) .* xval'];
    plot_stdShade('dataMat', lin_mat, 'xVal', xval, 'clr', sclr{istate},...
        'alpha', 0.2, 'axh', axh)
end
axis tight
ylim([0 1])
ylabel('Drift Rate [1 / hr]')
xlabel(sprintf('\x394 Time [h]'))
legend(snames)

axh = nexttile(th, 2, [1, 1]); cla; hold on
dataMat = squeeze(drft.drate(iunit, :, :));
plot_boxMean('axh', axh, 'dataMat', dataMat', 'plotType', 'bar',...
    'allPnts', false, 'clr', [0.8, 0.8, 0.8])
set(get(gca, 'children'), 'HandleVisibility', 'off');
plot([1 : 3], dataMat)
legend(mnames)
xlim([0.5, 3.5])
xticks([1 : 3])
xticklabels(snames)

