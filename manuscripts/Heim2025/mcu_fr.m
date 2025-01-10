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
frml = 'FR ~ Group + UnitType + (1|Mouse)';

% organize for lme
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
fr_tbl = lme_tbl;
lme = fitlme(fr_tbl, lme_cfg.frml);

% plot
mcu_lmePlot(fr_tbl, lme)

% copy results to excel
exlTbl = lme2exl(lme);

mcu_lmeContrasts(lme)


% FR per unit, WT vs MCU across states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'FR ~ Group * State + (1|Mouse)';

% organize for lme
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
iunit = 2;
fr_tbl = lme_tbl(lme_tbl.UnitType == iunit, :);
lme = fitlme(fr_tbl, lme_cfg.frml);

% plot
mcu_lmePlot(fr_tbl, lme)


% FR per unit per bout, WT vs MCU across states 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'FR ~ Group * State + BoutLength + (1|Mouse) + (1|UnitID)';

% organize for lme
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
iunit = 1;
fr_tbl = lme_tbl(lme_tbl.UnitType == iunit, :);
lme = fitlme(fr_tbl, lme_cfg.frml);

% plot
mcu_lmePlot(fr_tbl, lme)

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
clear grppaths
for igrp = 1 : length(grps)

    mnames = mcu_sessions(grps(igrp));
    
    for imouse = 1 : length(mnames)
        basepaths = mcu_sessions(mnames{imouse});
        grppaths{igrp}(imouse, :) = string(basepaths)';
    end
end

% FR ~ Group * Session + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize for lme
frml = 'FR ~ Group * Day + (1|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, true);

% run lme
iunit = 2;
fr_tbl = lme_tbl(lme_tbl.UnitType == iunit, :);
lme = fitlme(fr_tbl, lme_cfg.frml);

% plot
mcu_lmePlot(fr_tbl, lme)


% copy results to excel
exlTbl = lme2exl(lme);




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




