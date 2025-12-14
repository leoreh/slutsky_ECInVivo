function mcu_spikes(alt, tblUnit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCU_SPIKES - Analyze in vivo MCU spike data
%
% Inputs:
%   alt - Analysis type: 1 for baseline, 2 for baclofen
%   tblUnit - Optional pre-loaded data. If empty, data will be loaded
%
% Usage:
%   mcu_spikes(1)  % Run baseline analysis, load data automatically
%   mcu_spikes(2)  % Run baclofen analysis, load data automatically
%   mcu_spikes(1, myData)  % Run baseline analysis with pre-loaded data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data if not provided
if nargin < 2 || isempty(tblUnit)
    switch alt
        case 1
            tblUnit = load_data(alt);
        case 2
            tblUnit = load_data(alt);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch alt
    case 1
        % BASELINE
        run_bsl(tblUnit);
    case 2
        % BACLOFEN
        run_bac(tblUnit);
    otherwise
        error('Invalid alt value');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASELINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_bsl(tblUnit)

% Select Params
varRsp = 'FR';
% varRsp = 'BRoy';
% varRsp = 'BLidor';
% varRsp = 'BSpks';
% varRsp = 'PRC';

clear cfgLme
cfgLme.contrasts = 'all';
cfgLme.contrasts = [1 : 9];
cfgLme.distribution = 'Normal';
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
[lmeStats, lmeMdl] = lme_analyse(tblUnit, frml, cfgLme);

% plot
hFig = lme_plot(tblUnit, lmeMdl, 'lmeStats', lmeStats,...
    'idxRow', [], 'ptype', 'bar', 'axShape', 'square');
hAx = gca;
[barIdx, barLbl] = lme_sigLines(lmeStats, lmeMdl, tblUnit,...
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
% fname = lme_frml2char(frml, 'rmRnd', true, 'resNew', '');
% % fname = [fname, '_uAlt', num2str(2)];
% [cfg] = mcu_cfg();
% fPath = fullfile(cfg.savepath, fPrfx);
% lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat', 'xlsx'},...
%     'tblUnit', tblUnit, 'lmeStats', lmeStats, 'lmeMdl', lmeMdl, 'fPath', fPath)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BACLOFEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_bac(tblUnit)

% Select Params
unitType = 'FS';
varRsp = 'FR';
% varRsp = 'BRoy';
% varRsp = 'BLidor';

clear cfgLme
cfgLme.contrasts = 'all';
cfgLme.distribution = 'Gamma';
idxRow = [];
if strcmp(varRsp, 'FR')
    lblRsp = 'Firing Rate (% BSL)';
    fPrfx = 'FR';
    if strcmp(unitType, 'FS')
        cfgLme.contrasts = [1 : 19];
        idxRow = [12, 16, 18, 6, 7];
    else
        cfgLme.contrasts = [1 : 19, 21];
        idxRow = [12 : 15, 20, 6, 18];
    end
elseif strcmp(varRsp, 'BRoy')
    fPrfx = 'Burst';
    lblRsp = 'Burstiness Index (% BSL)';
    if strcmp(unitType, 'pPYR')
        cfgLme.contrasts = [1 : 19];
        idxRow = [4, 7, 16, 18];
    else
        idxRow = [1];
    end
end
if strcmp(varRsp, 'BLidor')
    cfgLme.distribution = 'Normal';
end

% Organize
iUnit = categorical({unitType});
tblLme = tblUnit(tblUnit.UnitType == iUnit, :);
tblLme = tbl_transform(tblLme, 'varNorm', 'Day',...
    'varsGrp', {'Group', 'UnitType'}, 'flgNorm', true);

% run lme
frml = [varRsp, ' ~ Group * Day + (Day|Name)'];
[lmeStats, lmeMdl] = lme_analyse(tblLme, frml, cfgLme);

% plot
hFig = lme_plot(tblLme, lmeMdl, 'lmeStats', lmeStats,...
    'idxRow', [], 'ptype', 'bar', 'axShape', 'wide');
hAx = gca;

% Generate significance lines
[barIdx, barLbl] = lme_sigLines(lmeStats, lmeMdl, tblUnit,...
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
    'tblUnit', tblLme, 'lmeStats', lmeStats, 'lmeMdl', lmeMdl, 'fPath', fPath);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD_DATA - Unified data loading function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tblUnit = load_data()

basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];

% Load table
tblUnit = mcu_unitTbl('basepaths', basepaths);

% Remove bad units
tblUnit(tblUnit.UnitType == 'Other', :) = [];
tblUnit.UnitType = removecats(tblUnit.UnitType, 'Other');

% Remove bac on / off
tblUnit(tblUnit.Day == 'BAC_ON', :) = [];
tblUnit(tblUnit.Day == 'BAC_OFF', :) = [];
tblUnit.Day = removecats(tblUnit.Day, {'BAC_ON', 'BAC_OFF'});

% Assert minimum value > zero



tblUnit.BRoy = tblUnit.BRoy + min(tblUnit.BRoy(tblUnit.BRoy > 0)) / 2;
tblUnit.BSpks = tblUnit.BSpks + min(tblUnit.BSpks(tblUnit.BSpks > 0)) / 2;
tblUnit.BMiz = tblUnit.BMiz + min(tblUnit.BMiz(tblUnit.BMiz > 0)) / 2;
tblUnit.FR = tblUnit.FR + min(tblUnit.FR(tblUnit.FR > 0)) / 2;


end


% % -------------------------------------------------------------------------
% % PLOT CORRELATIONS
% % Prepare data
% varsInc = {'FR', 'BSpks', 'BRoy', 'BMiz'};
% lData = tbl_transform(tblUnit, 'varsInc', varsInc, 'flgZ', false,...
%     'skewThr', 0.1, 'varsGrp', {'Group'}, 'flgLog', true);
% cfg = mcu_cfg();
%
% varsInc = [varsInc, {'BLidor', 'PRC'}];
% [hFig, ~, hGrid] = plot_corrHist(lData, 'varsInc', varsInc,...
%     'grpIdx', 'Group', 'clrGrp', cfg.clr.grp, 'thrOut', 100);
% plot_axSize('hFig', hFig, 'szOnly', false,...
%     'axWidth', 1200, 'axHeight', 600, 'flgPos', true);