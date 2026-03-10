
%% ========================================================================
% FIGURE 2
% =========================================================================

% A -----------------------------------------------------------------------

% Load
presets = {'steadyState'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

% FR
varRsp = 'fr';
frml = [varRsp, ' ~ Group + (1|Name)'];
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tbl, frml, 'dist', 'log-normal');

% Fit
varRsp = 'pBspk';
frml = [varRsp, ' ~ Group * fr + (1|Name)'];
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tbl, frml, 'dist', 'logit-normal');