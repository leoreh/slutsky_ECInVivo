

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate psd per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go over each mouse and analyze all experiment days
grps = [mcu_sessions('wt'), mcu_sessions('mcu')];

% go over baseline for wt vs. mcu. during experiment, only psdEmg should be
% calculated
grps = {'mcu_bsl'; 'wt_bsl'};
flg_bsl = true;
flg_bsl = false;

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
                'graphics', true, 'forceA', true, 'ftarget', ftarget,...
                'flgEmg', false);
        end

        % calc psd accordign to emg state separation
        psdEmg = psd_states('basepath', basepath, 'sstates', [1, 2],...
            'sig', sig, 'fs', fs, 'saveVar', true,...
            'graphics', true, 'forceA', true, 'ftarget', ftarget,...
            'flgEmg', true);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate fooof per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go over baseline for wt vs. mcu. during experiment
grps = {'wt_bsl'; 'mcu_bsl'};
flg_emg = false;

% go over each mouse and analyze all experiment days
grps = [mcu_sessions('wt'), mcu_sessions('mcu')];
flg_emg = true;

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
        

        % -----------------------------------------------------------------
        % create n PSDs per state by averaging k bouts (random sampling)
        k = 100;        % Number of bouts to sample per average
        npsd = 20;      % Number of averaged PSDs per state

        % Initialize
        nStates = numel(psd.bouts.psd);
        psd_avg = cell(1, nStates);
        psd_avg = cellfun(@(x) zeros(npsd, length(freqs)), psd_avg, 'uni', false);

        % Loop over each state
        for istate = 1 : nStates
            [nBouts, ~] = size(psd.bouts.psd{istate});

            for ipsd = 1 : npsd
                % Randomly sample k bouts with replacement
                sample_idx = randi(nBouts, k, 1);
                sampled_psds = psd.bouts.psd{istate}(sample_idx, :);

                % Compute the mean across sampled bouts
                psd_avg{istate}(ipsd, :) = mean(sampled_psds, 1);
            end
        end
     
        % -----------------------------------------------------------------
        % calculate fooof
        psd_1of = psd_fooof(psd_avg, freqs,...
            'flg_plot', false, 'saveVar', saveVar);

        % sample and plot several bouts
        % bout_idx = randperm(npsd, 5);
        % if ~isempty(bout_idx)
        %     psd_1of = psd_fooof(psd_avg{istate}(bout_idx, :), freqs,...
        %         'flg_plot', true, 'saveVar', false);
        % end
    end
end

fh = figure;
plot(psd.info.faxis, psd.psd(2, :))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOOOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grps = {'wt_bsl'; 'mcu_bsl'};

clear grppaths
for igrp = 1 : length(grps)
    grppaths{igrp} = string(mcu_sessions(grps{igrp})');
end
grppaths{2}(end) = [];

% organize for lme
frml = 'FOOOF ~ Group * State + (1|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
% istate = 1;
% band_tbl = lme_tbl(lme_tbl.State == categorical(istate), :);
lme = fitlme(lme_tbl, lme_cfg.frml);

% plot
fh = mcu_lmePlot(lme_tbl, lme, 'ptype', 'line');


% Baclofen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for each group
grps = {'wt', 'mcu'};
idx_rmDays = [];              % remove bac on and off
clear grppaths
for igrp = 1 : length(grps)

    mnames = mcu_sessions(grps(igrp));  
    for imouse = 1 : length(mnames)
        basepaths = mcu_sessions(mnames{imouse});
        basepaths(idx_rmDays) = [];
        grppaths{igrp}(imouse, :) = string(basepaths)';
    end
end


% FOOOF ~ Group * Day + (1|Mouse)
% organize for lme
frml = 'FOOOF ~ Group * Day + (Day|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, true, 'powsum', 1);

% select unit
istate = categorical({'High EMG'});
plot_tbl = lme_tbl(lme_tbl.State == istate, :);
[plot_tbl, norm_cfg] = mcu_lmeNorm(plot_tbl, 'Day', 'BSL');

% run lme
lme = fitlme(plot_tbl, lme_cfg.frml, 'FitMethod', 'REML');

% plot
fh = mcu_lmePlot(plot_tbl, lme, 'ptype', 'line');

% save
grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})

[prism_data] = fh2prism(fh);
exlTbl = lme2exl(lme, true);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot fooof results across time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for each group
grps = {'wt', 'mcu'};
idx_rmDays = [2, 6];              % remove bac on and off
clear grppaths
for igrp = 1 : length(grps)

    mnames = mcu_sessions(grps(igrp));  
    for imouse = 1 : length(mnames)
        basepaths = mcu_sessions(mnames{imouse});
        basepaths(idx_rmDays) = [];
        grppaths{igrp}(imouse, :) = string(basepaths)';
    end
end
ndays = length(basepaths);
if ndays == 7
    str_day = {'BSL'; 'BAC On'; 'BAC1'; 'BAC2'; 'BAC3'; 'BAC Off'; 'WASH'};
elseif ndays == 5
    str_day = {'BSL'; 'BAC1'; 'BAC2'; 'BAC3'; 'WASH'};
elseif ndays == 1
    str_day = {'BSL'};
else
    str_day = split(num2str(1 : ndays));
end

% vars
flg_emg = true;

if flg_emg
    vars = {'psdEmg_1of'};
else
    vars = {'psd_1of'};
end

% Load and organize vars for specific mouse
igrp = 1;
mice_idx = [1 : 5];
mice_idx = 1;
clear psd_1of_mouse
for imouse = 1 : length(mice_idx)
    basepaths = grppaths{igrp}(mice_idx(imouse), :);
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
    psd_1of_mouse(imouse) = catfields([v(:).psd_1of], 'addim', true);
end
psd_1of = catfields([psd_1of_mouse(:)], 2, true);
npnts = size(psd_1of.freqs, 2);
freqs = squeeze(psd_1of_mouse(imouse).freqs(1, :, 1));

% select params
fooof_param = 'residuals';

% prepare colors 
rgb_variants = rgbVariants([0.9, 0.1, 1], ndays);

% open figure
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [1, 2];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'FOOOF', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

for istate = 1 : 2
    % plot param
    axh = nexttile(th, istate, [1, 1]); cla; hold on
    data_mat = squeeze(psd_1of.(fooof_param)(istate, :, :, :));
    for iday = 1 : ndays
        data_plot = squeeze(data_mat(:, :, iday));
        ph(iday) = plot_stdShade('dataMat', data_plot', 'xVal', freqs,...
            'axh', axh, 'clr', rgb_variants(iday, :), 'alpha', 0.5);
    end

    % formatting
    % if ~contains(fooof_param, 'residuals')
    %     set(axh, 'YScale', 'log')
    % end
    set(axh, 'XScale', 'log')
    xlabel('Frequency [Hz]')
    ylabel('Power (log10)')
    legend(ph, str_day)
end

prism_data = fh2prism(fh);


xx = psd_fooof(psd_avg, freqs,...
    'flg_plot', false, 'saveVar', false);
fh = psd_fooofPlot(xx, 'istate', 2);









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CBPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% compare baseline psd wt vs mcu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab data
[psd_data, psd_cfg] = mcu_psdOrg(2, false);
nstates = length(psd_cfg.sstates);

% open figure
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [1, nstates];
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
    [stats, ~, tbl_stats] = mcu_psdCBPT(psd_state, psd_cfg);
    fpath = 'C:\Users\Leore\Downloads';
    writetable(tbl_stats, fullfile(fpath, 'cbpt_results.xlsx'), 'Sheet', psd_cfg.snames{istate});

    mcu_psdPlot(stats, psd_state, psd_cfg, 'axh', axh, 'ptype', 'freq',...
        'clr', psd_cfg.clr{istate})

    title(axh, psd_cfg.snames{istate})
end

% save
grph_save('fh', fh, 'fname', 'PSD_BSL_EMG', 'frmt', {'ai', 'jpg', 'fig'})
prism_data = fh2prism(fh);



% compare high vs low emg, for wt and mcu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab data - baseline
[psd_data, psd_cfg] = mcu_psdOrg(2, true);

% grab data - bac1
[psd_data, psd_cfg] = mcu_psdOrg(3, true);

grpnames = {'WT', 'MCU-KO'};

% open figure
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [1, 2];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'PSD Analysis', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

for igrp = 1 : 2

% psd per state
axh = nexttile(th, igrp, [1, 1]); cla; hold on

psd_grp = cell(1, 2);
for istate = 1 : 2
    psd_grp{istate} = squeeze(psd_data{igrp}(istate, :, :));
end

% analyze psd difference between groups for this state
[stats, ~, tbl_stats] = mcu_psdCBPT(psd_grp, psd_cfg);
fpath = 'C:\Users\Leore\Downloads';
writetable(tbl_stats, fullfile(fpath, 'cbpt_results.xlsx'), 'Sheet', grpnames{igrp});

mcu_psdPlot(stats, psd_grp, psd_cfg, 'axh', axh, 'ptype', 'freq',...
    'clr', psd_cfg.clr{istate})

title(axh, grpnames{igrp})

end

% save
grph_save('fh', fh, 'fname', 'PSD_BAC1_EMG', 'frmt', {'ai', 'jpg', 'fig'})
prism_data = fh2prism(fh);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare psd across time for WT mice (emg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab data
[psd_data, psd_cfg] = mcu_psdOrg(1);
[stats, cfg] = mcu_psdCBPT(psd_data, psd_cfg);

% extract params
ndays = size(psd_data{1}, 3);
grps = psd_cfg.snames;
clr = psd_cfg.clr;
faxis = psd_cfg.faxis;
txt_yax = psd_cfg.txt_yax;

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

prismData = fh2prism(fh);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare EMG and AS at baseline for wt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grab data (as)
[psd_data, psd_cfg] = mcu_psdOrg(2, false);
nstates = length(psd_cfg.sstates);

% grab data (emg)
[psd_data2, psd_cfg2] = mcu_psdOrg(2, true);
nstates = length(psd_cfg.sstates);

% organize data
psd_data{2} = psd_data2{1};
psd_data{1}(3, : ,:) = [];

% open figure
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [1, 2];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'PSD Analysis', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% psd per state
for istate = 1 : 2
    axh = nexttile(th, istate, [1, 1]); cla; hold on

    psd_state = cell(1, 2);
    for igrp = 1 : length(psd_data)
        psd_state{igrp} = squeeze(psd_data{igrp}(istate, :, :));
    end

    % analyze psd difference between groups for this state
    % [stats, ~] = mcu_psdCBPT(psd_state, psd_cfg);

    mcu_psdPlot([], psd_state, psd_cfg, 'axh', axh, 'ptype', 'freq',...
        'clr', psd_cfg.clr{istate})

    title(axh, psd_cfg.snames{istate})
end

% save
grph_save('fh', fh, 'fname', 'PSD_BSL_EMGvAS', 'frmt', {'ai', 'jpg', 'fig'})

prism_data = fh2prism(fh);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for each group
grps = {'wt', 'mcu'};
idx_rmDays = [];              % remove bac on and off
clear grppaths
for igrp = 1 : length(grps)

    mnames = mcu_sessions(grps(igrp));  
    for imouse = 1 : length(mnames)
        basepaths = mcu_sessions(mnames{imouse});
        basepaths(idx_rmDays) = [];
        grppaths{igrp}(imouse, :) = string(basepaths)';
    end
end

% Band ~ Group * Session + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize for lme
frml = 'Band ~ Group * Day + BoutLength + (1|Mouse) + (1|Day)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, true);

% select state
istate = categorical({'High EMG'});
plot_tbl = lme_tbl(lme_tbl.State == istate, :);

% % select group
% igrp = categorical({'WT'});
% plot_tbl = plot_tbl(plot_tbl.Group == igrp, :);
% 
% % update formula
% lme_cfg.frml = 'Band ~ Day + (1|Mouse) + (1|Day)';

% run lme
lme = fitlme(plot_tbl, lme_cfg.frml);

% plot
fh = mcu_lmePlot(plot_tbl, lme, 'ptype', 'line');

% save
frml = [char(lme.Formula), '_', char(istate)];
frml = frml2char(frml, 'rm_rnd', false);
th = get(gcf, 'Children'); % if you're sure it's a tiled layout
title(th, frml, 'interpreter', 'none')
grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})

prism_data = fh2prism(fh);



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











