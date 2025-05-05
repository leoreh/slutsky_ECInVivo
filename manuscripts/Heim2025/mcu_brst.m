% mcu_fr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate burstiness per file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go over each mouse and analyze all experiment days
grps = [mcu_sessions('wt'), mcu_sessions('mcu')];
vars = {'spikes'};

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

        % brst (mea)
        brst = spktimes_meaBrst(v(ifile).spikes.times, 'binsize', [],...
            'isiThr', 0.02, 'minSpks', 2, 'bins', [0 Inf],...
            'saveVar', true, 'flg_force', true, 'flg_all', false);

        % spike timing metrics
        st = spktimes_metrics('spktimes', v(ifile).spikes.times, 'sunits', [],...
            'bins', [0 Inf], 'flg_force', true, 'saveVar', true, 'flg_all', false);

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


% Burst ~ Group + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'Burst ~ Group + UnitType + (1|Mouse)';

% organize for lme
var_field = 'royer';
[lme_tbl, lme_cfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flg_emg', false, 'var_field', var_field, 'vCell', {});

% run lme
lme = fitlme(lme_tbl, lme_cfg.frml);

% plot
fh = lme_plot(lme_tbl, lme, 'ptype', 'bar');
grph_save('fh', fh, 'fname', [frml, '_spkprct'], 'frmt', {'ai', 'jpg'})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis during baclofen
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

% Burstiness of wt vs mcu across days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize for lme
frml = 'Burst ~ Group * Day + (1|Mouse)';

% organize for lme
var_field = 'lidor';
[lme_tbl, lme_cfg] = lme_org('grppaths', grppaths, 'frml', frml,...
    'flg_emg', false, 'var_field', var_field, 'vCell', {});

% select unit
iunit = categorical({'pPYR'});
plot_tbl = lme_tbl(lme_tbl.UnitType == iunit, :);

% run lme
lme = fitlme(plot_tbl, lme_cfg.frml);

% plot
fh = lme_plot(plot_tbl, lme, 'ptype', 'line');


% update formula
lme_cfg.frml = 'Burst ~ Day + (1|Mouse)';

% select group
igrp = categorical({'WT'});
plot_tbl = plot_tbl(plot_tbl.Group == igrp, :);

% run lme
lme = fitlme(plot_tbl, lme_cfg.frml);

% plot
fh = lme_plot(plot_tbl, lme, 'ptype', 'line');

frml = [char(lme.Formula), '_', char(iunit)];
frml = frml2char(frml, 'rm_rnd', false);
th = get(gcf, 'Children');
title(th, frml, 'interpreter', 'none')
grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})




