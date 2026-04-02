function [tbl, xVec, basepaths, v] = mcu_tblMea(varargin)
% MCU_TBLMEA Loads and organizes MEA data for MCU experiments.
%
%   [tbl, xVec, basepaths, v] = MCU_TBLMEA(...) loads data from BAC and
%   MCU-KO directories and organizes it into a single table. Input 'presets'
%   determine which variables are loaded and included in the output table.
%
%   INPUTS:
%       (Key-Value Pairs)
%       'basepaths' - (cell) List of recording folders. If empty, loads defaults.
%       'v'         - (struct) Pre-loaded data struct.
%       'presets'   - (cell) List of data types to include:
%                     'time'        : Include temporal dynamics (burstDyn).
%                     'spktimes'    : Include spike times (mea.spktimes).
%                     'steadyState' : Include steady-state recovery metrics (ss_).
%                     'frNet'       : Include network metrics (dim, mcc, cc).
%                     Note: Core firing stats (rate, dur, etc) are always included.
%
%   OUTPUTS:
%       tbl         - (table) Combined data table.
%       xVec        - (vector) Time vector for dynamics (hours).
%       basepaths   - (cell) Paths used.
%       v           - (struct) Raw loaded data.
%
%   See also: V2TBL, BASEPATHS2VARS, MCU_BASEPATHS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'basepaths', {});
addParameter(p, 'v', []);
addParameter(p, 'presets', {});
addParameter(p, 'flgOtl', true);
parse(p, varargin{:});

basepaths = p.Results.basepaths;
v         = p.Results.v;
presets   = p.Results.presets;
flgOtl    = p.Results.flgOtl;

% Default Basepaths
if isempty(basepaths)
    basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];
end

nFiles = length(basepaths);
cfg = mcu_cfg;


%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Define Variables to Load based on Presets
% Always include core vars
vars = {'fr', 'rcv', 'stats'};

if ismember('time', presets)
    vars = [vars, {'burstDyn'}];
end

if ismember('rcv', presets)
    vars = [vars, {'rcv_mdl'}];
end

if ismember('spktimes', presets)
    vars = [vars, {'mea', 'burst'}];
end

if ismember('frNet', presets)
    vars = [vars, {'frNet'}];
end

% Load if v is not provided
if isempty(v)
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
end


%% ========================================================================
%  INITIALIZE VARMAP & TAGS
%  ========================================================================

varMap = struct();

% Core Unit Properties (Always Included)
% --------------------------------------
varMap.uGood     = 'fr.uGood';

% Full Recovery Metrics
% -----------
if ismember('rcv', presets)
    rcvFile = 'rcv';
    varMap.fr          = [rcvFile, '.frBsl'];
    varMap.ss_fr        = [rcvFile, '.frSs'];
    varMap.frAcute     = [rcvFile, '.frAcute'];
    varMap.frTrough    = [rcvFile, '.frTrough'];
    varMap.pertDepth   = [rcvFile, '.pertDepth'];
    varMap.rcvErr      = [rcvFile, '.rcvErr'];
    varMap.rcvBsl      = [rcvFile, '.rcvBsl'];
    varMap.rcvGain     = [rcvFile, '.rcvGain'];
    varMap.rcvWork     = [rcvFile, '.rcvWork'];
    varMap.rcvDiff     = [rcvFile, '.rcvDiff'];
    varMap.bslTime     = [rcvFile, '.bslTime'];
    varMap.rcvTime     = [rcvFile, '.rcvTime'];
    varMap.rcvSlope    = [rcvFile, '.rcvSlope'];
    varMap.normSlope   = [rcvFile, '.normSlope'];
    varMap.spkDfct     = [rcvFile, '.spkDfct'];
    varMap.uRcv        = [rcvFile, '.uRcv'];
    varMap.uPert       = [rcvFile, '.uPert'];
else
    varMap.fr        = 'rcv.frBsl';
    varMap.ss_fr      = 'rcv.frSs';
    varMap.frAcute   = 'rcv.frAcute';
    varMap.frTrough  = 'rcv.frTrough';
    varMap.uRcv      = 'rcv.uRcv';
    varMap.rcvBsl    = 'rcv.rcvBsl';
    varMap.spkDfct   = 'rcv.spkDfct';
    varMap.uPert     = 'rcv.uPert';
end

varMap.br        = 'stats.br';
varMap.bDur      = 'stats.dur';
varMap.bFreq     = 'stats.freq';
varMap.bIBI      = 'stats.ibi';
varMap.pBurst    = 'stats.pBurst';
varMap.bSize     = 'stats.bSize';
varMap.frTot     = 'stats.fr';
varMap.frBurst   = 'stats.frBurst';
varMap.frSingle  = 'stats.frSingle';

% Spike Times
% -----------
if ismember('spktimes', presets)
    varMap.spktimes  = 'mea.spktimes';
    varMap.brstTimes = 'burst.spktimes';
end

% File Tags
% ---------
tagFiles.sbjID = get_mname(basepaths, 0);
tagFiles.genotype = repmat(cfg.lbl.grp(1), 1, nFiles);
tagFiles.genotype(contains(tagFiles.sbjID, 'ko')) = cfg.lbl.grp(2);


%% ========================================================================
%  BUILD TABLES
%  ========================================================================

% Main Table
% ----------
tbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles, ...
    'idxCol', 1, 'uOffset', 0);


% Steady State Table
% ---------------------------------
if ismember('steadyState', presets)
    mapSS = struct();
    mapSS.ss_br    = 'stats.br';
    mapSS.ss_bDur  = 'stats.dur';
    mapSS.ss_bFreq = 'stats.freq';
    mapSS.ss_bIBI  = 'stats.ibi';
    mapSS.ss_pBurst = 'stats.pBurst';
    mapSS.ss_bSize  = 'stats.bSize';
    mapSS.ss_frBurst = 'stats.frBurst';
    mapSS.ss_frSingle = 'stats.frSingle';

    tblSS = v2tbl('v', v, 'varMap', mapSS, 'tagFiles', tagFiles, ...
        'idxCol', 3, 'uOffset', 0);

    % Join
    tbl = outerjoin(tbl, tblSS, 'MergeKeys', true);
end


% Temporal Dynamics
% -----------------
xVec = [];
if ismember('time', presets)
    % Create a map for alignment
    mapDyn = struct();
    mapDyn.t_fr    = 'fr.fr';
    mapDyn.t_br    = 'burstDyn.br';
    mapDyn.t_bDur  = 'burstDyn.dur';
    mapDyn.t_bFreq = 'burstDyn.freq';
    mapDyn.t_bIBI  = 'burstDyn.ibi';
    mapDyn.t_pBurst = 'burstDyn.pBurst';
    mapDyn.t_bSize  = 'burstDyn.bSize';
    mapDyn.t_frTot = 'burstDyn.fr';
    mapDyn.t_frBurst = 'burstDyn.frBurst';
    mapDyn.t_frSingle = 'burstDyn.frSingle';

    % Align dynamics
    [v, t] = mea_tAlign(v, mapDyn, 'fr.info.idxPert');
    xVec = t / 3600;

    % Use separate v2tbl call with idxCol=[] to get full vectors
    tblDyn = v2tbl('v', v, 'varMap', mapDyn, 'tagFiles', tagFiles, ...
        'idxCol', [], 'uOffset', 0);

    % Join
    tbl = outerjoin(tbl, tblDyn, 'MergeKeys', true);
end


% Network Metrics (frNet)
% -----------------------
if ismember('frNet', presets)

    % Preprocess v to expand file-level metrics to unit-level for v2tbl
    for iFile = 1:length(v)
        if ~isfield(v(iFile), 'frNet') || isempty(v(iFile).frNet)
            continue;
        end

        frNet = v(iFile).frNet;
        zMet = 'shuffle';

        % Determine number of units from a known unit variable
        % (v(i).fr.uGood is reliable as it's used in the main table)
        nUnits = length(v(iFile).fr.uGood);

        % Expand Dimensionality (take 1st chunk/index)
        valDim = frNet.dim(1);
        v(iFile).frNet.dimExp = repmat(valDim, nUnits, 1);

        % Expand Mean Correlation
        valMcc = frNet.corr.(zMet).mcc(1);
        v(iFile).frNet.mccExp = repmat(valMcc, nUnits, 1);
    end

    mapNet = struct();
    mapNet.dim = 'frNet.dimExp';
    mapNet.mcc = 'frNet.mccExp';
    mapNet.funcon  = ['frNet.corr.', zMet, '.funconAvg'];
    mapNet.funcon_fish  = ['frNet.corr.', 'fisher', '.funconAvg'];
    % mapNet.funcon_raw  = ['frNet.corr.', 'raw', '.funcon'];

    tblNet = v2tbl('v', v, 'varMap', mapNet, 'tagFiles', tagFiles, ...
        'idxCol', 1, 'uOffset', 0);

    % Join
    tbl = outerjoin(tbl, tblNet, 'MergeKeys', true);
end


%% ========================================================================
%  CLEAN & POST-PROCESS
%  ========================================================================

tbl(~tbl.uGood, :) = [];
tbl = removevars(tbl, 'uGood');

tbl.unitID = categorical(tbl.unitID);

% Post-Process Time
if ismember('rcv', presets)
    tbl.rcvTime = tbl.rcvTime / 3600;
end

% Post-Process spktimes
if ismember('spktimes', presets)
    winExp = [0, 33600];
    tbl.spktimes = cellfun(@(x) x(x >= winExp(1) & x <= winExp(2)), ...
        tbl.spktimes, 'UniformOutput', false);
    tbl.brstTimes = cellfun(@(x) x(x >= winExp(1) & x <= winExp(2)), ...
        tbl.brstTimes, 'UniformOutput', false);
end

% Outlier Removal. Note: Manual inspection of FR traces confirmed these
% units are unstable.
if flgOtl

    % Clean units that were not perturbed
    tbl(~tbl.uPert, :) = [];
    tbl = removevars(tbl, 'uPert');

    % Clean residual outliers
    frml = 'ss_fr ~ (frBurst + frSingle) + (1 | sbjID)';

    % Control
    tblWt  = tbl(tbl.genotype == 'Control', :);
    lmeMdl = lme_analyse(tblWt, frml, 'dist', 'log-normal', 'verbose', false);
    res    = residuals(lmeMdl, 'ResidualType', 'Pearson');
    otlWt  = abs(res) > 3;
    idWt   = tblWt.unitID(otlWt);

    % MCU-KO
    tblMcu = tbl(tbl.genotype == 'MCU-KO', :);
    lmeMdl = lme_analyse(tblMcu, frml, 'dist', 'log-normal', 'verbose', false);
    res    = residuals(lmeMdl, 'ResidualType', 'Pearson');
    otlMcu = abs(res) > 3;
    idMcu  = tblMcu.unitID(otlMcu);

    % Remove
    otlIDs = [idMcu; idWt];
    otlIdx = ismember(tbl.unitID, otlIDs);
    tbl(otlIdx, :) = [];
end


end     % EOF