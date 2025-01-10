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
        % brst = spktimes_meaBrst(v(ifile).spikes.times, 'binsize', [],...
        %     'isiThr', 0.02, 'minSpks', 2, 'bins', [0 Inf],...
        %     'saveVar', true, 'flg_force', true, 'flg_all', false);

        % spike timing metrics
        st = spktimes_metrics('spktimes', v(ifile).spikes.times, 'sunits', [],...
            'bins', [0 Inf], 'flg_force', true, 'saveVar', true, 'flg_all', false);

    end
end




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

% Brst ~ Group * Session + (1|Mouse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize for lme
frml = 'Burst ~ Group * Day + (1|Mouse)';
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, true);

% run lme
iunit = 2;
brst_tbl = lme_tbl(lme_tbl.UnitType == iunit, :);
lme = fitlme(brst_tbl, lme_cfg.frml);

% plot
mcu_lmePlot(brst_tbl, lme)




