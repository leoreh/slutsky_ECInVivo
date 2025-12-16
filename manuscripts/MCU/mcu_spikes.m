
%  MCU_SPIKES - Analyze in vivo MCU spike data

%% ========================================================================
%  LOAD_DATA 
%  ========================================================================

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

% Assert no zero values
tblLme = tbl_transform(tblLme, 'flg0', true, 'verbose', true);


%% ========================================================================
%  BACLOFEN (y ~ Group * Day + (Day|Name))
%  ========================================================================

% Select Params
unitType = 'RS';
varRsp = 'FR';
varRsp = 'BRoy';

% Configuration
clear cfgLme
cfgLme.dist = 'Gamma';
cfgLme.contrasts = 1 : 19;

% Select data
tblLme = tblUnit(tblUnit.UnitType == unitType, :);

% run lme
frml = [varRsp, ' ~ Group * Day + (Day|Name)'];
[lmeStats, lmeMdl] = lme_analyse(tblLme, frml, cfgLme);

% Plot
hFig = tblGUI_bar(tblLme, 'yVar', varRsp, 'xVar', 'Day', 'grpVar', 'Group');


% Check best model
statsPark = lme_parkTest(tblLme, frml);
statsDist = lme_compareDists(tblLme, frml);






%% ========================================================================
%  BASELINE (y ~ Group * UnitType + (1|Name))
%  ========================================================================


% Select Params
varRsp = 'FR';
varRsp = 'BRoy';

% Configuration
clear cfgLme
cfgLme.contrasts = 1 : 9;
cfgLme.dist = 'Gamma';

% Select data
tblLme = tblUnit(tblUnit.Day == 'BSL', :);

% Fit
frml = [varRsp, ' ~ Group * UnitType + (1|Name)'];
[ , lmeMdl] = lme_analyse(tblUnit, frml, cfgLme);

% Plot
hFig = tblGUI_bar(tblLme, 'yVar', varRsp, 'xVar', 'UnitType', 'grpVar', 'Group');







