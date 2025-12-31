
%% ========================================================================
%  LOAD
%  ========================================================================

basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];
nFiles = length(basepaths);

vars = {'mea', 'fr', 'brstDyn', 'brst', 'frRcv', 'frRcv_mdl', ...
    'stats', 'ca', 'prc'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

cfg = mcu_cfg;


% -------------------------------------------------------------------------
% DOWNSAMPLE

for iFile = 1 : nFiles

    sig = v(iFile).ca.mito;
    [nUnits, nSig] = size(sig);
    nBins = floor(nSig / 60);
    sig = sig(:, 1 : nBins * 60);
    sig = reshape(sig, nUnits, 60, nBins);
    sig = reshape(sum(sig, 2), nUnits, nBins);
    sig = fr_denoise(sig, 'frameLen', 30, 'flgPlot', false);
    v(iFile).ca.mito = sig;

    sig = v(iFile).ca.cyto;
    sig = sig(:, 1 : nBins * 60);
    sig = reshape(sig, nUnits, 60, nBins);
    sig = reshape(sum(sig, 2), nUnits, nBins);
    sig = fr_denoise(sig, 'frameLen', 30, 'flgPlot', false);
    v(iFile).ca.cyto = sig;

end

% -------------------------------------------------------------------------
% TABLE

% Vars to align
varMap = struct();
varMap.frt = 'fr.fr';
varMap.btRate = 'brstDyn.rate';
varMap.btDur = 'brstDyn.dur';
varMap.btFreq = 'brstDyn.freq';
varMap.btIBI = 'brstDyn.ibi';
varMap.btFrac = 'brstDyn.bfrac';
varMap.caCyto = 'ca.cyto';
varMap.caMito = 'ca.mito';

% Align
[v, t] = mea_tAlign(v, varMap, 'fr.info.idxPert');
xVec = t / 3600;

% Unit vars
varMap.uGood = 'fr.uGood';
varMap.uRcv = 'rcv.uRcv';
varMap.uPert = 'rcv.uPert';

% Tag structures
tagFiles.Name = get_mname(basepaths, 0);
tagFiles.Group = repmat(cfg.lbl.grp(1), 1, nFiles);
tagFiles.Group(contains(tagFiles.Name, 'KO')) = cfg.lbl.grp(2);

% Table
tbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles, ...
    'uOffset', 0);

% -------------------------------------------------------------------------
% COMBINE & CLEAN

% Clean bad units
tbl(~tbl.uPert, :) = [];
tbl(~tbl.uGood, :) = [];

tbl.UnitID = categorical(tbl.UnitID);

%% ========================================================================
%  PLOT
%  ========================================================================

experiments = categories(tbl.Name);
idxExp = tbl.Name == experiments(2);
tblGUI_xy(xVec, tbl(idxExp, :), 'tileVar', 'Group', 'yVar', 'caMito');
