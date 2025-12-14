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
unitType = 'RS';
varRsp = 'FR';
% varRsp = 'BRoy';
% varRsp = 'BLidor';

% Configuration
clear cfgLme
cfgLme.contrasts = 'all';
cfgLme.dist = 'Gamma';
cfgLme.contrasts = [1 : 19];

% Select data and normalize to baseline
iUnit = categorical({unitType});
tblLme = tblUnit(tblUnit.UnitType == iUnit, :);
tblLme = tbl_transform(tblLme, 'varNorm', 'Day',...
    'varsGrp', {'Name', 'UnitType'}, 'flgNorm', true);

% Check best model
statsPark = lme_parkTest(tblLme, frml);
statsDist = lme_compareDists(tblLme, frml);

% run lme
frml = [varRsp, ' ~ Group * Day + (Day|Name)'];
[lmeStats, lmeMdl] = lme_analyse(tblLme, frml, cfgLme);

hFig = tblGUI_bar(tblLme, 'yVar', 'FR', 'xVar', 'Day', 'grpVar', 'Group');



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
vars = tblUnit.Properties.VariableNames;
for iVar = 1:numel(vars)
    if isnumeric(tblUnit.(vars{iVar}))
        if any(tblUnit.(vars{iVar}) == 0)
            halfMin = min(tblUnit.(vars{iVar})(tblUnit.(vars{iVar}) > 0)) / 2;
            tblUnit.(vars{iVar}) = tblUnit.(vars{iVar}) + halfMin;
            fprintf('Warning: %s contained zeros. Added half-min value (%.5f).\n', vars{iVar}, halfMin);
        end
    end
end


end

