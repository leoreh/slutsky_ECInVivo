
pathName = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results';
xlsName = 'mcu_suppTbl.xlsx';

%% ========================================================================
% FIGURE 2
% =========================================================================

% A -----------------------------------------------------------------------

% Load
presets = {'steadyState'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

% LME FR
frml = 'fr ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tbl, frml, 'dist', 'log-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('2A (FR)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

% LME P_burst
frml = 'pBspk ~ Group * fr + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tbl, frml, 'dist', 'logit-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('2A (P_burst)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


%% ========================================================================
% FIGURE 3
% =========================================================================

% H -----------------------------------------------------------------------

% Load
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
presets = {'brst'};
tbl = mcu_tblVivo('basepaths', basepaths, 'flgClean', true, ...
    'presets', presets);

% Assert no zero values
tblLme = tbl_trans(tbl, 'flg0', true, 'verbose', true);

% Limit to RS units
uIdx = tblLme.unitType == 'RS';
tblLme = tblLme(uIdx, :);

% Remove WASH
tblLme(tblLme.Day == 'WASH', :) = [];
tblLme.Day = removecats(tblLme.Day, {'WASH'});

% LME
frml = 'fr ~ Group * Day + (Day|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'gamma');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('3H', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

