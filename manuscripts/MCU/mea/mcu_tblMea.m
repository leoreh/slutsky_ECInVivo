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
%                     'time'        : Include temporal dynamics (brstDyn).
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
parse(p, varargin{:});

basepaths = p.Results.basepaths;
v         = p.Results.v;
presets   = p.Results.presets;

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
    vars = [vars, {'brstDyn'}];
end

if ismember('rcv', presets)
    vars = [vars, {'rcv_mdl'}];
end

if ismember('spktimes', presets)
    vars = [vars, {'mea'}];
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
    varMap.frSs        = [rcvFile, '.frSs'];
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
    varMap.frSs      = 'rcv.frSs';
    varMap.frTrough  = 'rcv.frTrough';
    varMap.uRcv      = 'rcv.uRcv';
    varMap.rcvBsl    = 'rcv.rcvBsl';
    varMap.spkDfct   = 'rcv.spkDfct';
    varMap.uPert     = 'rcv.uPert';
end

varMap.bRate     = 'stats.eventRate';
varMap.bDur      = 'stats.dur';
varMap.bFreq     = 'stats.freq';
varMap.bIBI      = 'stats.ibi';
varMap.pBspk     = 'stats.pBspk';
varMap.nBspk     = 'stats.nBspk';
varMap.frTot     = 'stats.frTot';
varMap.frBspk    = 'stats.frBspk';
varMap.frSspk    = 'stats.frSspk';

% Spike Times
% -----------
if ismember('spktimes', presets)
    varMap.spktimes = 'mea.spktimes';
end

% File Tags
% ---------
tagFiles.Name = get_mname(basepaths, 0);
tagFiles.Group = repmat(cfg.lbl.grp(1), 1, nFiles);
tagFiles.Group(contains(tagFiles.Name, 'ko')) = cfg.lbl.grp(2);


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
    mapSS.ss_bRate = 'stats.eventRate';
    mapSS.ss_bDur  = 'stats.dur';
    mapSS.ss_bFreq = 'stats.freq';
    mapSS.ss_bIBI  = 'stats.ibi';
    mapSS.ss_pBspk = 'stats.pBspk';
    mapSS.ss_nBspk = 'stats.nBspk';

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
    mapDyn.t_bRate = 'brstDyn.eventRate';
    mapDyn.t_bDur  = 'brstDyn.dur';
    mapDyn.t_bFreq = 'brstDyn.freq';
    mapDyn.t_bIBI  = 'brstDyn.ibi';
    mapDyn.t_pBspk = 'brstDyn.pBspk';
    mapDyn.t_nBspk = 'brstDyn.nBspk';
    mapDyn.t_frTot = 'brstDyn.frTot';
    mapDyn.t_frBspk = 'brstDyn.frBspk';
    mapDyn.t_frSspk = 'brstDyn.frSspk';

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
    mapNet.funcon  = ['frNet.corr.', zMet, '.funcon'];
    % mapNet.funcon_fish  = ['frNet.corr.', 'fisher', '.funcon'];
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

tbl.UnitID = categorical(tbl.UnitID);

% Post-Process Time
if ismember('rcv', presets)
    tbl.rcvTime = tbl.rcvTime / 3600;
end

% Post-Process spktimes
if ismember('spktimes', presets)
    winExp = [0, 33600];
    tbl.spktimes = cellfun(@(x) x(x >= winExp(1) & x <= winExp(2)), ...
        tbl.spktimes, 'UniformOutput', false);
end


end     % EOF