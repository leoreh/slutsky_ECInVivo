function mcu_spikes(alt, lmeData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCU_SPIKES - Analyze in vivo MCU spike data
%
% Inputs:
%   alt - Analysis type: 1 for baseline, 2 for baclofen
%   lmeData - Optional pre-loaded data. If empty, data will be loaded
%
% Usage:
%   mcu_spikes(1)  % Run baseline analysis, load data automatically
%   mcu_spikes(2)  % Run baclofen analysis, load data automatically
%   mcu_spikes(1, myData)  % Run baseline analysis with pre-loaded data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data if not provided
if nargin < 2 || isempty(lmeData)
    switch alt
        case 1
            lmeData = load_data(alt);
        case 2
            lmeData = load_data(alt);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch alt
    case 1
        % BASELINE
        run_bsl(lmeData);
    case 2
        % BACLOFEN
        run_bac(lmeData);
    otherwise
        error('Invalid alt value');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASELINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_bsl(lmeData)

% Select Params
varRsp = 'FR';
% varRsp = 'BRoy';
% varRsp = 'BLidor';
% varRsp = 'BSpks';
varRsp = 'PRC';

clear lmeCfg
lmeCfg.contrasts = 'all';
lmeCfg.contrasts = [1 : 9];
lmeCfg.distribution = 'Normal';
if strcmp(varRsp, 'FR')
    fPrfx = 'FR';
    lblY = 'Firing Rate (Hz)';
    idxRow = [3, 6];
elseif strcmp(varRsp, 'BRoy')
    lblY = 'Burstiness Index (a.u.)';
    fPrfx = 'Burst';
    idxRow = [6, 3, 9];
elseif strcmp(varRsp, 'PRC')
    lblY = 'Population Coupling (a.u.)';
    fPrfx = 'PRC';
    idxRow = [];
end

% Fit
frml = [varRsp, ' ~ Group * UnitType + (1|Name)'];
[lmeStats, lmeMdl] = lme_analyse(lmeData, frml, lmeCfg);

% plot
hFig = lme_plot(lmeData, lmeMdl, 'lmeStats', lmeStats,...
    'idxRow', [], 'ptype', 'bar', 'axShape', 'square');
hAx = gca;
[barIdx, barLbl] = lme_sigLines(lmeStats, lmeMdl, lmeData,...
    'idxRow', idxRow);
plot_sigLines(hAx, barIdx, barLbl, 'flgNS', true);

% Update labels
ylabel(hAx, lblY)
xlabel(hAx, '')
title(hAx, '')
hAx.Legend.Visible = "off";
hAx.Legend.Location = 'best';
hAx.Legend.Position(1) = hAx.Position(1) + hAx.Position(3) - hAx.Legend.Position(3);
hAx.Legend.Position(2) = hAx.Position(2);
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'square', 'axHeight', 300)

% Save
fname = lme_frml2char(frml, 'rmRnd', true, 'resNew', '');
% fname = [fname, '_uAlt', num2str(2)];
[cfg] = mcu_cfg();
fPath = fullfile(cfg.savepath, fPrfx);
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', lmeData, 'lmeStats', lmeStats, 'lmeMdl', lmeMdl, 'fPath', fPath)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BACLOFEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_bac(lmeData)

% Select Params
unitType = 'pINT';
varRsp = 'FR';
% varRsp = 'BRoy';
% varRsp = 'BLidor';

clear lmeCfg
lmeCfg.contrasts = 'all';
lmeCfg.distribution = 'Gamma';
idxRow = [];
if strcmp(varRsp, 'FR')
    lblRsp = 'Firing Rate (% BSL)';
    fPrfx = 'FR';
    if strcmp(unitType, 'pINT')
        lmeCfg.contrasts = [1 : 19];
        idxRow = [12, 16, 18, 6, 7];
    else
        lmeCfg.contrasts = [1 : 19, 21];
        idxRow = [12 : 15, 20, 6, 18];
    end
elseif strcmp(varRsp, 'BRoy')
    fPrfx = 'Burst';
    lblRsp = 'Burstiness Index (% BSL)';
    if strcmp(unitType, 'pPYR')
        lmeCfg.contrasts = [1 : 19];
        idxRow = [4, 7, 16, 18];
    else
        idxRow = [1];
    end
end
if strcmp(varRsp, 'BLidor')
    lmeCfg.distribution = 'Normal';
end

% Organize
iUnit = categorical({unitType});
plotTbl = lmeData(lmeData.UnitType == iUnit, :);
plotTbl = tbl_transform(plotTbl, 'varNorm', 'Day',...
    'varsGrp', {'Group', 'UnitType'}, 'flgNorm', true);

% run lme
frml = [varRsp, ' ~ Group * Day + (Day|Name)'];
[lmeStats, lmeMdl] = lme_analyse(plotTbl, frml, lmeCfg);

% plot
hFig = lme_plot(plotTbl, lmeMdl, 'lmeStats', lmeStats,...
    'idxRow', [], 'ptype', 'bar', 'axShape', 'wide');
hAx = gca;

% Generate significance lines
[barIdx, barLbl] = lme_sigLines(lmeStats, lmeMdl, lmeData,...
    'idxRow', idxRow);
plot_sigLines(hAx, barIdx, barLbl, 'flgNS', true);

% Graphics
ylabel(hAx, lblRsp)
xlabel(hAx, '')
title(hAx, unitType)
plot_axSize('hFig', hFig, 'szOnly', false, 'axShape', 'wide', 'axHeight', 300);
hAx.Legend.Location = 'northeast';
hAx.Legend.Position(1) = hAx.Position(1) + hAx.Position(3) - hAx.Legend.Position(3);
hAx.Legend.Position(2) = hAx.Position(2);
if strcmp(unitType, 'pINT')
    hAx.Legend.Visible = "off";
end

% save
fname = lme_frml2char(frml, 'rmRnd', true, 'resNew', '',...
    'sfx', [' _', unitType, '_Norm']);
fname = [fname, '_uAlt', num2str(2)];
[cfg] = mcu_cfg();
fPath = fullfile(cfg.savepath, fPrfx);
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
    'lmeData', plotTbl, 'lmeStats', lmeStats, 'lmeMdl', lmeMdl, 'fPath', fPath);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD_DATA - Unified data loading function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lmeData = load_data(alt)

% Get labels and parameters
[cfg] = mcu_cfg();
cfg.varMap.UnitType = 'units.clean';

% Parse additional parameters
if alt == 2
    idxRm = [2, 6];
end

% Initialize
clear tblCell v tagAll tagFiles

if alt == 1
    % BASELINE DATA LOADING
    grps = {'wt_bsl'; 'mcu_bsl'};

    for iGrp = 1 : length(grps)
        basepaths = mcu_sessions(grps{iGrp});
        nFiles = length(basepaths);
        v{iGrp} = basepaths2vars('basepaths', basepaths, 'vars', cfg.vars);

        tagAll.Group = cfg.lbl.grp{iGrp};
        tagFiles.Name = get_mname(basepaths);
        tblCell{iGrp} = v2tbl('v', v{iGrp}, 'varMap', cfg.varMap, 'tagAll',...
            tagAll, 'tagFiles', tagFiles, 'idxCol', 2);
    end

elseif alt == 2
    % BACLOFEN DATA LOADING
    grps = {'wt', 'mcu'};

    for iGrp = 1 : length(grps)
        mNames = mcu_sessions(grps{iGrp});

        for iMouse = 1 : length(mNames)
            % Get all paths for this mouse and remove specified days
            tmpPaths = mcu_sessions(mNames{iMouse});
            tmpPaths(idxRm) = [];

            % Load data for all days for this mouse at once
            v{iGrp, iMouse} = basepaths2vars('basepaths', tmpPaths, 'vars', cfg.vars);

            % Prepare tag structures for this mouse
            tagAll.Group = cfg.lbl.grp{iGrp};
            tagAll.Name = mNames{iMouse};
            tagFiles.Day = cfg.lbl.day(1:length(tmpPaths));

            % Create table for this mouse using new flexible approach
            tblCell{iGrp, iMouse} = v2tbl('v', v{iGrp, iMouse}, 'varMap', cfg.varMap, ...
                'tagFiles', tagFiles, 'tagAll', tagAll, 'idxCol', 2);
        end
    end
end

% Combine all tables
lmeData = vertcat(tblCell{:});

% Clean up and organize
lmeData.UnitType = categorical(cfg.lbl.unit(lmeData.UnitType + 1)');
lmeData = rmmissing(lmeData);
lmeData.Group = reordercats(lmeData.Group, cfg.lbl.grp);
lmeData.UnitType = reordercats(lmeData.UnitType, cfg.lbl.unit);

% Additional cleanup for baclofen data
if alt == 2
    lmeData.Day = reordercats(lmeData.Day, cfg.lbl.day);
    lmeData.BRoy = lmeData.BRoy + min(lmeData.BRoy(lmeData.BRoy > 0)) / 2;
    lmeData.BSpks = lmeData.BSpks + min(lmeData.BSpks(lmeData.BSpks > 0)) / 2;
    lmeData.BMiz = lmeData.BMiz + min(lmeData.BMiz(lmeData.BMiz > 0)) / 2;
    lmeData.FR = lmeData.FR + min(lmeData.FR(lmeData.FR > 0)) / 2;
end

end


% % -------------------------------------------------------------------------
% % PLOT CORRELATIONS
% % Prepare data
% varsInc = {'FR', 'BSpks', 'BRoy', 'BMiz'};
% lData = tbl_transform(lmeData, 'varsInc', varsInc, 'flgZ', false,...
%     'skewThr', 0.1, 'varsGrp', {'Group'}, 'flgLog', true);
% cfg = mcu_cfg();
%
% varsInc = [varsInc, {'BLidor', 'PRC'}];
% [hFig, ~, hGrid] = plot_corrHist(lData, 'varsInc', varsInc,...
%     'grpIdx', 'Group', 'clrGrp', cfg.clr.grp, 'thrOut', 100);
% plot_axSize('hFig', hFig, 'szOnly', false,...
%     'axWidth', 1200, 'axHeight', 600, 'flgPos', true);