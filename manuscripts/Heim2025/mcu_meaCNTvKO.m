
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE BASELINE VS RECOVERY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Load data from both groups
grps = {'mea_bac'; 'mea_mcuko'};
grpLbls = {'Control'; 'MCU-KO'};
vars = {'frr', 'st_metrics', 'st_brst'};
vars = {'frr', 'st_brst'};

% Create varMap based on time point
clear varMap
varMap.uGood      = 'frr.uGood';
% varMap.BrMiz      = 'st.mizuseki';
% varMap.BrRoy      = 'st.royer';
varMap.BSpks      = 'brst.bspks';

% Override the source prefix (in main analysis)
mdlPrfx = 'frr.mdlF';

cnt = 1;
clear tblCell
for iGrp = 1 : length(grps)
    basepaths = mcu_sessions(grps{iGrp});
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
   
    % Create consistent UnitIDs for this group (reset between groups). Use
    % frr for grabbing the number of units.
    frr = catfields([v(:).frr], 1);
    [nUnits, ~] = size(frr.fr);
    unitIDs = (1:nUnits)';

    for iCol = 1 : 3

        % Fr maps to different sources based on time point
        switch iCol
            case 1
                varMap.Fr  = [mdlPrfx, '.frBsl'];
                tagAll.Time = 'BSL';
            case 2
                varMap.Fr  = [mdlPrfx, '.frTrough'];
                tagAll.Time = 'BAC 1h';

            case 3
                varMap.Fr  = [mdlPrfx, '.frSs'];
                tagAll.Time = 'BAC 24h';
        end

        % Prepare tag structures for v2tbl
        tagAll.Group = grpLbls{iGrp};
        tagFiles.Name = get_mname(basepaths);

        tempTbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles',...
            tagFiles, 'tagAll', tagAll, 'idxCol', iCol);
        
        % Override the UnitID to ensure consistency across time points
        % Add group offset to make UnitIDs unique between groups
        groupOffset = (iGrp - 1) * 10000;  % 10000 units per group
        tempTbl.UnitID = groupOffset + unitIDs;

        tblCell{cnt} = tempTbl;
        cnt = cnt + 1;
    end
end
tbl = vertcat(tblCell{:});
tbl = sortrows(tbl, 'Group');


% -------------------------------------------------------------------------
% Bursty
varRsp = 'BSpks';

% Organize for analysis
lmeData = tbl(tbl.uGood, :);
lmeData.uGood = [];
lmeData = rmmissing(lmeData);

% Normalize to baseline
nData = tbl_transform(lmeData, 'varsInc', {varRsp}, 'flgZ', false,...
    'skewThr', 2, 'varsGrp', {'Group'}, 'varNorm', 'Time',...
    'flgLog', false, 'flgNorm', false);
nData.BSpks = nData.BSpks + min(nData.BSpks(nData.BSpks > 0)) / 2;

% Remove unused categorical levels
% nData(nData.Time == 'BAC 24h', :) = [];
nData(nData.Time == 'BAC 1h', :) = [];
if ismember('Time', nData.Properties.VariableNames)
    actualTimeLevels = unique(nData.Time);
    nData.Time = categorical(nData.Time, actualTimeLevels);
end
grpstats(lmeData, {'Group', 'Time'});

% run lme
frml = [varRsp, ' ~ Group * Time + (1|Name)'];
lmeCfg.contrasts = [1 : 7];
lmeCfg.distribution = 'Gamma';
[lmeStats, lmeMdl] = lme_analyse(nData, frml, lmeCfg);

% plot
hFig = lme_plot(nData, lmeMdl, 'lmeStats', lmeStats,...
    'ptype', 'bar', 'axShape', 'square', 'idxRow', [4, 3]);

% Update labels
hAx = gca;
ylabel(hAx, 'Burstiness P(Spk\inBrst)', 'Interpreter', 'tex')
xlabel(hAx, '')
title(hAx, '')
hAx.Legend.Location = 'southeast';
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square',...
    'axHeight', 300)



% -------------------------------------------------------------------------
% Firing Rate

% Organize for analysis
lmeData = tbl(tbl.uGood, :);
lmeData.uGood = [];
lmeData = rmmissing(lmeData);

% Normalize to baseline
nData = tbl_transform(lmeData, 'varsInc', {varRsp}, 'flgZ', false,...
    'skewThr', 2, 'varsGrp', {'Group'}, 'varNorm', 'Time',...
    'flgLog', false, 'flgNorm', false);
nData.Fr = nData.Fr + min(nData.Fr(nData.Fr > 0)) / 2;

% Remove unused categorical levels
% nData(nData.Time == 'BSL', :) = [];
nData(nData.Time == 'BAC 1h', :) = [];
if ismember('Time', nData.Properties.VariableNames)
    actualTimeLevels = unique(nData.Time);
    nData.Time = categorical(nData.Time, actualTimeLevels);
end
grpstats(lmeData, {'Group', 'Time'});

% run lme
varRsp = 'Fr';
frml = [varRsp, ' ~ Group * Time + (1|Name)'];
lmeCfg.contrasts = [];
lmeCfg.distribution = 'Gamma';
[lmeStats, lmeMdl] = lme_analyse(nData, frml, lmeCfg);

% plot
hFig = lme_plot(nData, lmeMdl, 'lmeStats', lmeStats,...
    'ptype', 'bar', 'axShape', 'square', 'idxRow', [1 : 7]);

% Update labels
hAx = gca;
ylabel(hAx, 'Firing Rate (Hz)', 'Interpreter', 'tex')
xlabel(hAx, '')
title(hAx, '')
hAx.Legend.Location = 'southeast';
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square',...
    'axHeight', 300)

