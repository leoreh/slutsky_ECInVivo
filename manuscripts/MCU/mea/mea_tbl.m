function [tbl, xVec, basepaths, v] = mea_tbl(basepaths, v)
% MEA_TBL loads data from MCU-KO and CONTROL experiments in mea and
% organizes it in a table 

%% ========================================================================
%  LOAD
%  ========================================================================

if nargin < 1 || isempty(basepaths)
    basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];
end
nFiles = length(basepaths);

if nargin < 2 || isempty(v)
    vars = {'mea', 'fr', 'brstDyn', 'brst', 'frRcv', 'frRcv_mdl', ...
        'stats', 'ca', 'prc'};
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
end

cfg = mcu_cfg;

% -------------------------------------------------------------------------
% TEMPORAL DYNAMICS

flgMito = false;

if flgMito

    % % Calculate Ca2+ stats per unit during windows
    % winBsl = [1 : v(1).fr.info.idxPert - 5];
    % cytoBsl = mean(tbl.caCyto(:, winBsl), 2, 'omitnan');
    % mitoBsl = mean(tbl.caMito(:, winBsl), 2, 'omitnan');
    %
    % tbl = addvars(tbl, cytoBsl, mitoBsl);

    varMap = struct();
    varMap.caCyto = 'ca.cyto';
    varMap.caMito = 'ca.mito';
    xVec = v(1).ca.time / 3600;
    
    for iFile = 1 : nFiles
        v(iFile).ca.cyto = v(iFile).ca.cyto(:, [1 : 32395]);
        v(iFile).ca.mito = v(iFile).ca.mito(:, [1 : 32395]);
        v(iFile).ca.time = v(iFile).ca.time(:, [1 : 32395]);
    end

else
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

end

% Unit vars
varMap.uGood = 'fr.uGood';
varMap.uRcv = 'rcv.uRcv';
varMap.uPert = 'rcv.uPert';

% Tag structures
tagFiles.Name = get_mname(basepaths, 0);
tagFiles.Group = repmat(cfg.lbl.grp(1), 1, nFiles);
tagFiles.Group(contains(tagFiles.Name, 'KO')) = cfg.lbl.grp(2);

% Table
tblT = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles, ...
    'uOffset', 0);


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
varMap.rcvErr = 'rcvMdl.rcvErr';
varMap.rcvDiff = 'rcvMdl.rcvDiff';
varMap.spktimes = 'mea.spktimes';

% PRC
varMap.prc = 'prc.prc0_z';

% Table
tbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles, ...
    'idxCol', 1, 'uOffset', 0);

% -------------------------------------------------------------------------
% BURSTINESS STEADY STATE

% Unit vars
varMap = struct();
varMap.uGood = 'fr.uGood';
varMap.uRcv = 'rcv.uRcv';
varMap.uPert = 'rcv.uPert';

% Bursts
varMap.bRateSS = 'stats.rate';
varMap.bDurSS = 'stats.brstDur';
varMap.bFreqSS = 'stats.freq';
varMap.bIBISS = 'stats.ibi';
varMap.bFracSS = 'stats.bspks';
varMap.bSpksSS = 'stats.nspks';

% Table
tbl2 = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles, ...
    'idxCol', 3, 'uOffset', 0);

% -------------------------------------------------------------------------
% COMBINE & CLEAN

% JOIN TABLES
tbl = outerjoin(tbl, tblT, 'MergeKeys', true);
tbl = outerjoin(tbl, tbl2, 'MergeKeys', true);

% Clean bad units
tbl(~tbl.uPert, :) = [];
tbl(~tbl.uGood, :) = [];

% Limit spktimes 
winExp = [0, 9 * 60]  * 60;
tbl.spktimes = cellfun(@(x) x(x >= winExp(1) & x <= winExp(2)), ...
    tbl.spktimes, 'UniformOutput', false);

% Convert time to hours
tbl.rcvTime = tbl.rcvTime / 3600;

% Manually add raw rcvErr
tbl.rcvErr = log2((tbl.frSs + 1 / 3600) ./ (tbl.fr + 1 / 3600));
tbl.rcvFrac = log2((tbl.bFracSS + 1e-6) ./ (tbl.bFrac + 1e-6));
tbl.rcvDur = log2((tbl.bDurSS + 1e-6) ./ (tbl.bDur + 1e-6));



end     % EOF