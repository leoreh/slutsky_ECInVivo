

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate psd per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go over each mouse and analyze all experiment days
grps = [mcu_sessions('wt'), mcu_sessions('mcu')];

% go over baseline for wt vs. mcu. during experiment, only psdEmg should be
% calculated
grps = {'mcu_bsl'; 'wt_bsl'};
flg_bsl = true;

% state params
ftarget = [1 : 0.5 : 100];

% iterate
for igrp = 1 : length(grps)

    queryStr = grps{igrp};
    basepaths = mcu_sessions(queryStr);
    nfiles = length(basepaths);
    mpath = fileparts(basepaths{1});

    for ifile = 1 : nfiles

        % files
        basepath = basepaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)
        sigfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

        % load lfp signal
        sig = load(sigfile, 'eeg');
        sig = sig.eeg;
        load(sigfile, 'fs');

        % calc psd
        if flg_bsl
            psd = psd_states('basepath', basepath, 'sstates', [1, 4, 5],...
                'sig', sig, 'fs', fs, 'saveVar', true,...
                'graphics', false, 'forceA', true, 'ftarget', ftarget,...
                'flgEmg', false);
        end

        % calc psd accordign to emg state separation
        psdEmg = psd_states('basepath', basepath, 'sstates', [1, 2],...
            'sig', sig, 'fs', fs, 'saveVar', true,...
            'graphics', false, 'forceA', true, 'ftarget', ftarget,...
            'flgEmg', true);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate fooof per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go over each mouse and analyze all experiment days
grps = [mcu_sessions('wt'), mcu_sessions('mcu')];
flg_emg = true;

% go over baseline for wt vs. mcu. during experiment, only psdEmg should be
% calculated
grps = {'wt_bsl'; 'mcu_bsl'};
flg_emg = false;

% vars
if flg_emg
    vars = {'psdEmg'};
    saveVar = 'psdEmg_1of';
else
    vars = {'psd'};
    saveVar = 'psd_1of';
end

% iterate
for igrp = 1 : length(grps)

    queryStr = grps{igrp};
    basepaths = mcu_sessions(queryStr);
    nfiles = length(basepaths);
    mpath = fileparts(basepaths{1});
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);

    for ifile = 1 : nfiles

        % files
        basepath = basepaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)
       
        % psd params
        psd = v(ifile).psd;
        freqs = psd.info.faxis;
        sstates = psd.info.sstates;
        nstates = length(sstates);
        clr = psd.info.clr;
        
        % sample and plot several bouts
        bout_idx = randperm(50, 0);
        if ~isempty(bout_idx)
            psd_1of = psd_fooof(psd.bouts.psd{istate}(bout_idx, :), freqs,...
                'flg_plot', true, 'saveVar', false);
        end
        
        % use psd averaged across bouts
        psd_1of = psd_fooof(num2cell(psd.psd, 2), freqs,...
            'flg_plot', true, 'saveVar', false);

        % calculate fooof
        psd_1of = psd_fooof(psd.bouts.psd, freqs,...
            'flg_plot', false, 'saveVar', saveVar);

        % plot
        fh = figure;
        set(fh, 'WindowState', 'maximized');
        tlayout = [1, nstates];
        th = tiledlayout(tlayout(1), tlayout(2));
        th.TileSpacing = 'tight';
        th.Padding = 'none';
        title(th, 'FOOOF', 'interpreter', 'none', 'FontSize', 20)
        set(fh, 'DefaultAxesFontSize', 16);
        
        % limit to boutlength
        blen = cellfun(@(x) diff(x, 1, 2), psd.bouts.times, 'uni', false);
        thr_blen = 0;
        for istate = 1 : nstates
            blen_idx = blen{istate} > thr_blen;
            axh = nexttile(th, istate, [1, 1]); cla; hold on
            psd_cell{1} = squeeze(psd_1of.psd_orig(istate, blen_idx, :))';
            psd_cell{2} = squeeze(psd_1of.psd_fit(istate, blen_idx, :))';
            mcu_psdPlot([], psd_cell, [], 'axh', axh, 'ptype', 'freq',...
                'flg_log', false, 'faxis', freqs, 'clr', clr{istate})
            ylim([1 5])
        end

    end
end


% FOOOF ~ Group + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grps = {'wt_bsl'; 'mcu_bsl'};

clear grppaths
for igrp = 1 : length(grps)
    grppaths{igrp} = string(mcu_sessions(grps{igrp})');
end


% organize for lme
frml = 'FOOOF ~ Group * State + (1|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
% istate = 1;
% band_tbl = lme_tbl(lme_tbl.State == categorical(istate), :);
lme = fitlme(lme_tbl, lme_cfg.frml);

% plot
mcu_lmePlot(lme_tbl, lme);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare baseline psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab data
[psd_data, psd_cfg] = mcu_psdOrg(2);
nstates = length(psd_cfg.sstates);

% open figure
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [2, nstates];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'PSD Analysis', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% psd per state
for istate = 1 : nstates
    axh = nexttile(th, istate, [1, 1]); cla; hold on

    psd_state = cell(1, 2);
    for igrp = 1 : length(psd_data)
        psd_state{igrp} = squeeze(psd_data{igrp}(istate, :, :));
    end

    % analyze psd difference between groups for this state
    [stats, ~] = mcu_psdCBPT(psd_state, psd_cfg);

    mcu_psdPlot(stats, psd_state, psd_cfg, 'axh', axh, 'ptype', 'freq',...
        'clr', psd_cfg.clr{istate})

    title(axh, psd_cfg.snames{istate})
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare psd across high and low EMG states, across time for WT mice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab data
[psd_data, psd_cfg] = mcu_psdOrg(1);
[stats, cfg] = mcu_psdCBPT(psd_data, psd_cfg);

% plot

% open figure
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [2, ndays];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'PSD Analysis', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% extract params
ndays = size(psd_data{1}, 3);
grps = psd_cfg.snames;
clr = psd_cfg.clr;
faxis = psd_cfg.faxis;
txt_yax = psd_cfg.txt_yax;

axh = nexttile(th, 1, [1, ndays]); cla; hold on
mcu_psdPlot(stats, psd_data, psd_cfg, 'axh', axh, 'ptype', 'heat', 'flg_log', true)

% plot psd across days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% psd per state
for iday = 1 : ndays
    axh = nexttile(th, iday + ndays, [1, 1]); cla; hold on

    psd_cfg.expDsgn = 2;
    psd_cfg.grps = psd_cfg.snames;
    psd_day = {squeeze(psd_data{1}(:, :, iday))', squeeze(psd_data{2}(:, :, iday))'};
    [stats, ~] = mcu_psdCBPT(psd_day, psd_cfg);
    mcu_psdPlot(stats, psd_day, psd_cfg, 'axh', axh, 'ptype', 'freq',...
        'flg_log', true)

    title(axh, sprintf('day %d', iday))
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for each group
grps = {'wt', 'mcu'};
clear grppaths
for igrp = 1 : length(grps)

    mnames = mcu_sessions(grps(igrp));

    for imouse = 1 : length(mnames)
        basepaths = mcu_sessions(mnames{imouse});
        grppaths{igrp}(imouse, :) = string(basepaths)';
    end
end

% Band ~ Group * Session + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize for lme
frml = 'Band ~ Group * Day + (1|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, true);

% run lme
istate = 2;
band_tbl = lme_tbl(lme_tbl.State == categorical(istate), :);
lme = fitlme(band_tbl, lme_cfg.frml);

% plot
mcu_lmePlot(band_tbl, lme);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bout Lengths vs. Band Power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for each group
grps = {'wt_bsl', 'mcu_bsl'};
clear grppaths
for igrp = 1 : length(grps)
    grppaths{igrp} = string(mcu_sessions(grps{igrp})');
end


% Band ~ Group * Session + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize for lme
frml = 'Band ~ Group * State * BoutLength + (1|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
lme = fitlme(lme_tbl, lme_cfg.frml);

% plot
mcu_lmePlot(lme_tbl, lme);


% correlation per state and group
for istate = 1 : nstates
    for igrp = 1 : ngroups

        grp_idx = lme_tbl.Group == categorical(igrp - 1) &...
            lme_tbl.State == categorical(istate);
        y_data = lme_tbl.Band(grp_idx);
        x_data = lme_tbl.BoutLength(grp_idx);

        % calculate correlation
        [r, p] = corr(x_data, y_data);
        r_mat(istate, igrp) = r;
        p_mat(istate, igrp) = p;
    end
end

fh = figure;
colors = {'b', 'r'};
for istate = 1 : nstates
    subplot(1, nstates, istate)
    hold on

    band_tbl = lme_tbl(lme_tbl.State == categorical(istate), :);

    % correlation
    for igrp = 1 : ngroups

        grp_idx = band_tbl.Group == categorical(igrp - 1);
        y_data = band_tbl.Band(grp_idx);
        x_data = log10(band_tbl.BoutLength(grp_idx));

        scatter(x_data, y_data, 15, 'filled', 'MarkerFaceAlpha', 0.2)

        % add regression line
        p = polyfit(x_data, y_data, 1);
        xfit = linspace(min(x_data), max(x_data), 100);
        yfit = polyval(p, xfit);
        plot(xfit, yfit, 'LineWidth', 2)
    end

    % aesthetics
    xlabel('Bout Length')
    ylabel('Band Power')
    title(sprintf('State %d', istate))
    box off
    axis tight
end











