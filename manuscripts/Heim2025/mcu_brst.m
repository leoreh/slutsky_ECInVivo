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
        brst = spktimes_meaBrst(v(ifile).spikes.times, 'binsize', [], 'isiThr', 0.02,...
            'minSpks', 2, 'saveVar', true, 'force', true, 'bins', [0 Inf]);

    end
end




v = basepaths2vars('basepaths', basepaths, 'vars', {'st_brst'});


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
[fr_tbl, fr_cfg] = mcu_frOrg(grppaths, frml, true);

% run lme
iunit = 1;
lme_tbl = fr_tbl(fr_tbl.UnitType == iunit, :);
lme = fitlme(lme_tbl, fr_cfg.frml);

% copy results to excel
exlTbl = lme2exl(lme);