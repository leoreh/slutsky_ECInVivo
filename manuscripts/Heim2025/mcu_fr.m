% mcu_fr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate fr per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uses bout times from psd (after cleaning bouts)

flg_emg = false;

% go over each mouse and analyze all experiment days
grps = [mcu_sessions('wt'), mcu_sessions('mcu')];

% go over baseline for wt vs. mcu. during experiment, only psdEmg should be
grps = {'mcu_bsl'; 'wt_bsl'};

% vars
if flg_emg
    vars = {'spikes'; 'psdEmg'; 'sleep_statesEmg'};
    fname_fr = 'frEmg';
else
    vars = {'spikes'; 'psd'; 'sleep_states'};
    fname_fr = 'fr';
end

% iterate
for igrp = 1 : length(grps)

    % get group basepaths
    queryStr = grps{igrp};
    basepaths = mcu_sessions(queryStr);
    nfiles = length(basepaths);
    
    % load state vars
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);

    for ifile = 1 : nfiles

        % files
        basepath = basepaths{ifile};
        [~, basename] = fileparts(basepath);
        cd(basepath)
        
        % bout time
        btimes = v(ifile).psd.bouts.times;

        % get mfr per bout
        fr = calc_fr(v(ifile).spikes.times, 'basepath', basepath,...
            'graphics', false, 'binsize', 60, 'saveVar', fname_fr, 'forceA', true,...
            'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf],...
            'btimes', btimes);             
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grps = {'wt_bsl'; 'mcu_bsl'};

clear grppaths
for igrp = 1 : length(grps)
    grppaths{igrp} = string(mcu_sessions(grps{igrp})');
end


% FR per unit, WT vs MCU for RS vs FS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'FR ~ Group * UnitType + (1|Mouse)';

% organize for lme
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
lme = fitlme(lme_tbl, lme_cfg.frml, 'FitMethod', 'REML');

% plot
fh = mcu_lmePlot(lme_tbl, lme);
grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})

prism_data = fh2prism(fh);
exlTbl = lme2exl(lme, true);


% FR per unit, WT vs MCU across states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'FR ~ Group * State + (1|Mouse)';

% organize for lme
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
iunit = categorical({'pPYR'});
plot_tbl = lme_tbl(lme_tbl.UnitType == iunit, :);
lme = fitlme(plot_tbl, lme_cfg.frml, 'FitMethod', 'REML');

% plot
fh = mcu_lmePlot(plot_tbl, lme, 'ptype', 'bar');

% save
frml = [char(lme.Formula), '_', char(iunit)];
frml = frml2char(frml, 'rm_rnd', false);
th = get(gcf, 'Children');
title(th, frml, 'interpreter', 'none')
axh = get(th, 'Children');
ylim(axh(3), [0 3])
grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})

prism_data = fh2prism(fh);
exlTbl = lme2exl(lme, true);



% Create contrast to test MCU vs WT during NREM
% Using the contrast matrix through coefTest or contrastTest
contrasts = zeros(1, length(lme.CoefficientNames));
% contrasts(strcmp(lme.CoefficientNames, 'Group_MCU-KO')) = 1;
contrasts(strcmp(lme.CoefficientNames, 'State_NREM')) = 1;
% contrasts(strcmp(lme.CoefficientNames, 'Group_MCU-KO:State_AW')) = 1;

% Test the contrast
[p, F, df] = coefTest(lme, contrasts)



% FR per unit per bout, WT vs MCU across states 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'FR ~ Group * State * BoutLength + (1|Mouse) + (1|UnitID)';

% organize for lme
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
iunit = categorical({'pPYR'});
fr_tbl = lme_tbl(lme_tbl.UnitType == iunit, :);
lme = fitlme(fr_tbl, lme_cfg.frml, 'FitMethod', 'REML');

% plot
fh = mcu_lmePlot(fr_tbl, lme, 'ptype', 'line');
% save
frml = [char(lme.Formula), '_', char(iunit)];
frml = frml2char(frml, 'rm_rnd', false);
th = get(gcf, 'Children');
title(th, frml, 'interpreter', 'none')
axh = get(th, 'Children');
ylim(axh(3), [0 3])
grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})


% Create contrast to test MCU vs WT during NREM
% Using the contrast matrix through coefTest or contrastTest
contrastInteraction = zeros(1, length(lme.Coefficients.Estimate));
contrastInteraction(strcmp(lme.CoefficientNames, 'Group_MCU-KO:State_NREM')) = 1;
contrastInteraction(strcmp(lme.CoefficientNames, 'Group_MCU-KO:State_REM')) = 1;

% Test the contrast
[p, F, df] = coefTest(lme, contrastInteraction)








% investigate relation BL and FR
[r_mat, p_mat, fh] = mcu_FRvBL(fr_tbl);

% limit bouts so that distribution is the same
[idx_eq, fh] = mcu_eqBLen(fr_tbl);

eq_tbl = fr_tbl(idx_eq, :);
lme = fitlme(eq_tbl, lme_cfg.frml);
mcu_lmePlot(eq_tbl, lme)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis during baclofen
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

% FR ~ Group * Day + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize for lme
frml = 'FR ~ Group * Day + (Day|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, true);

% select unit
iunit = categorical({'pPV'});
plot_tbl = lme_tbl(lme_tbl.UnitType == iunit, :);

% run lme
lme = fitlme(plot_tbl, lme_cfg.frml);

% plot
fh = mcu_lmePlot(plot_tbl, lme, 'ptype', 'bar');

% save
frml = [char(lme.Formula), '_', char(iunit)];
frml = frml2char(frml, 'rm_rnd', false);
th = get(gcf, 'Children');
title(th, frml, 'interpreter', 'none')
grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})

[prism_data] = fh2prism(fh);
exlTbl = lme2exl(lme, true);

% test interaction of wt and mcu across all days
interactionIndices = contains(lme.CoefficientNames, 'Group_MCU-KO:Day'); % Logical vector for interaction terms
contrastVector = zeros(1, length(lme.Coefficients.Estimate));
contrastVector(interactionIndices) = 1; % Test only interaction terms

% Run the test
[p, h, stat] = coefTest(lme, contrastVector);



% FR ~ Day * UnitType + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'FR ~ Day + (1|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths(2), frml, true);

% run lme
iunit = categorical({'pPYR'});
plot_tbl = lme_tbl(lme_tbl.UnitType == iunit, :);
lme = fitlme(plot_tbl, lme_cfg.frml);

% plot
fh = mcu_lmePlot(plot_tbl, lme, 'ptype', 'bar');

frml = [char(lme.Formula), '_', char(iunit), '_MCU-KO'];
frml = frml2char(frml, 'rm_rnd', false);
th = get(gcf, 'Children');
title(th, frml, 'interpreter', 'none')
grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})


th = get(gcf, 'Children');
axh = get(th, 'Children');
axes(axh(3))
yyaxis right
ylim([0, 9])
yyaxis left
ylim([0, 5])

grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})


% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate mean and standard error
mfr = grpstats(lme_tbl, {'Group', 'State'}, {'mean', 'sem'}, 'DataVars', 'FR');

% Extract data for plotting
states = categories(mfr.State);  % Get state labels
groups = categories(mfr.Group);  % Group labels (WT, MCU-KO)

% Prepare data for RS units
frData = reshape(mfr.mean_FR, length(states), length(groups));
frErr = reshape(mfr.sem_FR, length(states), length(groups));

% Plot RS Units
figure;
bar(frData);
hold on;
ngroups = size(frData, 1);
nbars = size(frData, 2);
% Calculate the positions for the error bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = (1:ngroups) - 0.3 + (i-1) * 0.3;
end
% Add error bars
errorbar(x', frData, frErr, 'k', 'linestyle', 'none');
xlabel('State');
ylabel('Firing Rate (Hz)');
legend(groups, 'Location', 'northwest');
hold off;








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mea example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot mea
load('E:\Data\MEA\baclofen\190813_022600\190813_022600.fr.mat')

fh = figure;
tlayout = [1, 1];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'FR MEA', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

axh = nexttile(th, 1, [1, 1]); cla; hold on
hold on
unit_idx = [2, 6, 23, 30];
frMat = fr.norm;
xvals = fr.tstamps * 3 / 60 / 60;
for iu = 1 : size(frMat, 1)
    frMat(iu, :) = smooth(frMat(iu, :), 15);
end
plot(xvals, frMat(unit_idx, :))
plot_stdShade('dataMat', frMat', 'xVal', xvals, ...
    'axh', axh, 'clr', [0.3 0.3 0.3], 'alpha', 0.3);
xlabel('Time (hr)')
ylim([0 4])
xlim([0 18])
xticks([0 : 6 : 24])
legend({'Unit 1', 'Unit 2', 'Unit 3', 'Unit 4', 'MFR'}, 'Location', 'northwest')

grph_save('fh', fh, 'fname', 'FR_MEA', 'frmt', {'ai', 'jpg'})
