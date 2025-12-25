function tbl = mea_tbl()
% MEA_TBL loads data from MCU-KO and CONTROL experiments in mea and
% organizes it in a table 

%% ========================================================================
%  LOAD
%  ========================================================================

vars = {'mea', 'fr', 'brstDyn', 'brst', 'frRcv', 'frRcv_mdl', 'stats'};
basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
nFiles = length(basepaths);

cfg = mcu_cfg;


% -------------------------------------------------------------------------
% TEMPORAL DYNAMICS

% Vars to align
varMap = struct();
varMap.frt = 'fr.fr';
varMap.btRate = 'brstDyn.rate';
varMap.btDur = 'brstDyn.dur';
varMap.btFreq = 'brstDyn.freq';
varMap.btIBI = 'brstDyn.ibi';
varMap.btFrac = 'brstDyn.bfrac';

% Align
[v, t] = mea_tAlign(v, varMap, 'fr.info.idxPert', true);
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
tblT = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles);


% -------------------------------------------------------------------------
% UNIT FIRING STATS

% Unit vars
varMap = struct();
varMap.uGood = 'fr.uGood';
varMap.uRcv = 'rcv.uRcv';
varMap.uPert = 'rcv.uPert';

% Bursts
varMap.bRate = 'stats.rate';
varMap.bDur = 'stats.brstDur';
varMap.bFreq = 'stats.freq';
varMap.bIBI = 'stats.ibi';
varMap.bFrac = 'stats.bspks';
varMap.bSpks = 'stats.nspks';

% FR
varMap.fr = 'rcv.frBsl';
varMap.frSs = 'rcv.frSs';
varMap.spkDfct = 'rcv.spkDfct';
varMap.rcvTime = 'rcvMdl.rcvTime';
varMap.spktimes = 'mea.spktimes';

% Tag structures
tagFiles.Name = get_mname(basepaths, 0);
tagFiles.Group = repmat(cfg.lbl.grp(1), 1, nFiles);
tagFiles.Group(contains(tagFiles.Name, 'KO')) = cfg.lbl.grp(2);

% Table
tbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles, 'idxCol', 1);

% -------------------------------------------------------------------------
% COMBINE & CLEAN

% JOIN TABLES
tbl = outerjoin(tbl, tblT, 'MergeKeys', true);

% Clean bad units
tbl(~tbl.uPert, :) = [];
tbl(~tbl.uGood, :) = [];

% Limit spktimes 
winExp = [0, 9 * 60]  * 60;
tbl.spktimes = cellfun(@(x) x(x >= winExp(1) & x <= winExp(2)), ...
    tbl.spktimes, 'UniformOutput', false);

% Convert time to hours
tbl.rcvTime = tbl.rcvTime / 3600;

% Normalize time
% winNorm = [1, v(1).fr.info.idxPert - 5];
% tblNorm = tbl_tNorm(tbl, {'frt', 'btRate', 'btDur', 'btFreq', 'btIBI', 'btFrac'}, ...
%     winNorm, {'Name'});

end     % EOF