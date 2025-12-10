
% groups
queryStr = ["wt_bsl"; "wt_bac2"; "wt_bac3"; "wt_wsh"];
ngrp = length(queryStr);

clear drft2 drft
for igrp = 1 : ngrp

    % load data
    basepaths = mcu_basepaths(queryStr{igrp});
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
    for istate = 1
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


prismData = reshape(tmp, 1, 10)


igrp = 2;
tmp = squeeze(drft.drate(:, [2, 3], :, igrp));
tmp = permute(tmp, [1, 3, 2]);
prismData = reshape(tmp, 2, 10)


