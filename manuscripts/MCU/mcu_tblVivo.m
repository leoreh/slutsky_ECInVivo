function [tbl, basepaths, v, xVec] = mcu_tblVivo(varargin)

% MCU_TBLVIVO Loads and processes unit data for the MCU project.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                loads defaults.
%   v            (struct) Pre-loaded data structure.
%   varMap       (struct) Variable mapping for table creation.
%   flgClean     (logical) Remove bad units and bac on / off.
%   presets      (cell) List of data types to include ('swv', 'prc', 'frNet').
%
% OUTPUT:
%   Tbl          (table) Unit table with metadata.
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'basepaths', {}, @iscell);
addParameter(p, 'v', [], @isstruct);
addParameter(p, 'varMap', [], @isstruct);
addParameter(p, 'flgClean', false, @islogical);
addParameter(p, 'presets', {}, @iscell);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
v         = p.Results.v;
varMap    = p.Results.varMap;
flgClean  = p.Results.flgClean;
presets   = p.Results.presets;

% Output
xVec = [];

%% ========================================================================
%  VAR MAP
%  ========================================================================

% Set varMap
cfg = mcu_cfg;
if isempty(varMap)
    varMap = cfg.varMap;
end
vars = cfg.vars;

% Presets
if ismember('swv', presets)
    vars = [vars, 'swv_metrics'];
    varMap.wv_tp = 'swv.tp';
    varMap.wv_asym = 'swv.asym';
    varMap.wv_hpk = 'swv.hpk';
end

if ismember('burst', presets)
    vars = [vars, 'burstStats'];
    varMap.br        = 'stats.br';
    varMap.bDur      = 'stats.dur';
    varMap.bFreq     = 'stats.freq';
    varMap.bIBI      = 'stats.ibi';
    varMap.pBurst    = 'stats.pBurst';
    varMap.bSize     = 'stats.bSize';
    varMap.frBurst   = 'stats.frBurst';
    varMap.frSingle  = 'stats.frSingle';
end

if ismember('spktimes', presets)
    vars = [vars, 'spikes'];
    varMap.spktimes = 'spikes.times';
end

if ismember('prc', presets)
    vars = [vars, 'prc'];
    varMap.PRC = 'prc.prc0_norm';
end

if ismember('frNet', presets)
    vars = [vars, 'drift', 'frNet'];
    zMet = 'shuffle';
    varMap.drift = 'drft.drftExp';
    varMap.dim   = 'frNet.dimExp';
    varMap.mcc   = 'frNet.mccExp';
    varMap.funcon_shf  = 'frNet.funcon_shf';
    varMap.funcon_fish  = 'frNet.funcon_fish';
    varMap.funcon_raw  = 'frNet.funcon_raw';
end

if ismember('rippSpks', presets)
    vars = [vars, 'rippSpks', 'rippSpkLfp'];
    spkType = '';
    varMap.frRipp = ['rippSpks.', spkType, '.frEvt'];
    varMap.frRand = ['rippSpks.', spkType, '.frCtrl'];
    varMap.frZ = ['rippSpks.', spkType, '.frZ'];
    varMap.frMod = ['rippSpks.', spkType, '.frMod'];
    varMap.pFire = ['rippSpks.', spkType, '.pFire'];
    varMap.cRipp = ['rippSpks.', spkType, '.cEvt'];
    varMap.rankMean = ['rippSpks.', spkType, '.rankMean'];
    varMap.rankVar = ['rippSpks.', spkType, '.rankVar'];
    varMap.frActive = ['rippSpks.', spkType, '.frActive'];
    varMap.cActive = ['rippSpks.', spkType, '.cActive'];
    varMap.frPre = ['rippSpks.', spkType, '.frPre'];
    varMap.frPost = ['rippSpks.', spkType, '.frPost'];
    varMap.asym = ['rippSpks.', spkType, '.asym'];
    varMap.com = ['rippSpks.', spkType, '.com'];
    varMap.peth = ['rippSpks.', spkType, '.peth'];

    phaseType = '';
    varMap.theta = ['spkLfp.', phaseType, '.theta'];
    varMap.mrl = ['spkLfp.', phaseType, '.mrl'];
    varMap.ppc = ['spkLfp.', phaseType, '.ppc1'];
end

if ismember('ripp', presets)
    vars = {'ripp'};
    varMap = struct();
    varMap.dur          = 'ripp.dur';
    varMap.amp          = 'ripp.amp';
    varMap.freq         = 'ripp.freq';
    varMap.freqEvent    = 'ripp.freqEvent';
    varMap.energy       = 'ripp.energy';
    varMap.state        = 'ripp.state';
    varMap.spkGain      = 'ripp.spkGain';
    varMap.skew         = 'ripp.skew';
    varMap.frac         = 'ripp.spks.RS.frac';
    varMap.asym         = 'ripp.spks.RS.asym';
    varMap.com          = 'ripp.spks.RS.com';
end

if ismember('rippMaps', presets)
    vars = [vars, 'rippMaps'];
    varMap = struct();
    varMap.t_lfp        = 'rippMaps.lfp';
    varMap.t_filt       = 'rippMaps.filt';
    varMap.t_amp        = 'rippMaps.amp';
    varMap.t_freq       = 'rippMaps.freq';
    varMap.t_z          = 'rippMaps.z';
    % population PETH (RS/FS/MU) is not precomputed here; the 3D raster in
    % rippSpkMaps.mat can supply a per-type PETH if a table later needs it.
end

if ismember('rippStates', presets)
    vars = {'rippStates'};
    varMap = struct();
    varMap.Rate         = 'rippStates.Rate';
    varMap.Density      = 'rippStates.Density';
    varMap.Duration     = 'rippStates.Duration';
    varMap.State        = 'rippStates.State';
end

if ismember('acg', presets)
    % acg_narrow: narrow autocorrelogram (100 ms, 0.5 ms bins).
    % st_metrics is already in cfg.vars; xVec is extracted post-load.
    varMap.acg_narrow = 'st.acg_narrow';
    varMap.acg_wide = 'st.acg_wide';

end

%% ========================================================================
%  LOAD DATA
%  ========================================================================

if isempty(basepaths)
    basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu'), mcu_basepaths('lh137')];
end
if isempty(v)
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
end

% Post-process frNet expansion
if ismember('frNet', presets) 
    for iFile = 1:length(v)
        if ~isfield(v(iFile), 'frNet') || isempty(v(iFile).frNet)
            continue;
        end

        nUnits = length(v(iFile).units.type);
        
        % Expand Drift (take 1st chunk/index)
        valDrft = v(iFile).drft.drate;
        v(iFile).drft.drftExp = repmat(valDrft, nUnits, 1);

        % Expand Dimensionality 
        valDim = v(iFile).frNet.dim(1);
        v(iFile).frNet.dimExp = repmat(valDim, nUnits, 1);

        % Expand Mean Correlation       
        valMcc = v(iFile).frNet.corr.(zMet).mcc(1);
        v(iFile).frNet.mccExp = repmat(valMcc, nUnits, 1);

        % Average funcon
        v(iFile).frNet.funcon_shf = mean(v(iFile).frNet.corr.shuffle.funcon, 2, 'omitnan');
        v(iFile).frNet.funcon_fish = mean(v(iFile).frNet.corr.fisher.funcon, 2, 'omitnan');
        v(iFile).frNet.funcon_raw = mean(v(iFile).frNet.corr.raw.funcon, 2, 'omitnan');

    end
end

% Post-process ripp states
if ismember('rippStates', presets)
    for iFile = 1:length(v)
        badIdx = ismember(v(iFile).rippStates.State, {'WAKE', 'N/REM', 'REM'});
        v(iFile).rippStates(badIdx, :) = [];
    end
end

% Post-process ripp: keep only accepted events. QA now MARKS events via
% .accepted (NREM/valid state, low EMG, above the MUA-gain gate) rather than
% removing them, so the per-event table is filtered here. Old .ripp.mat files
% without .accepted are left untouched (they already hold only the survivors).
if ismember('ripp', presets)
    for iFile = 1:length(v)
        if isfield(v(iFile), 'ripp') && isstruct(v(iFile).ripp) ...
                && isfield(v(iFile).ripp, 'accepted')
            v(iFile).ripp = evt_subset(v(iFile).ripp, v(iFile).ripp.accepted);
        end
    end
end

% Extract xVec for ripples
if ismember('rippSpks', presets)
    xVec = v(1).rippSpks.tstamps;
end
if ismember('rippMaps', presets)
    xVec = v(1).rippMaps.tstamps;
end

% Extract xVec for ACG (lag axis in ms, same for every recording)
if ismember('acg', presets)
    xVec.narrow = v(1).st.info.acg_narrow_tstamps * 1000;   % [s] → [ms]
    xVec.wide = v(1).st.info.acg_wide_tstamps * 1000;   % [s] → [ms]
end

%% ========================================================================
%  TABLE
%  ========================================================================

% Metadata
tagFiles = struct();
tagFiles.sbjID = get_mname(basepaths);
[~, fileNames] = fileparts(basepaths);
tagFiles.fileID = fileNames;

% Table
tbl = v2tbl('v', v, 'varMap', varMap, 'tagAll',...
    struct(), 'tagFiles', tagFiles, 'idxCol', []);


%% ========================================================================
%  PROCESS METADATA
%  ========================================================================

% Group metadata
tbl.genotype = ones(height(tbl), 1) * 1;
tbl.genotype(ismember(tbl.sbjID, cfg.miceMCU), :) = 2;
tbl.genotype = categorical(tbl.genotype, [1, 2], cfg.lbl.grp);

% Day metadata
fileTbl = unique(tbl(:, {'sbjID', 'fileID'}), 'rows');
fileGrp = findgroups(fileTbl.sbjID);
dayCell = splitapply(@(x) {(1:numel(x))'}, fileTbl.fileID, fileGrp);
fileTbl.day = vertcat(dayCell{:});
tbl = join(tbl, fileTbl, 'Keys', {'sbjID', 'fileID'});
tbl.day = categorical(tbl.day, [1 : 7], cfg.lbl.day);

% Reorder columns
tblVars = tbl.Properties.VariableNames;
if any(contains(tblVars, "unitType"))
    varOrder = {'genotype', 'sbjID', 'fileID', 'day', 'unitID', 'unitType'};
else
    varOrder = {'genotype', 'sbjID', 'fileID', 'day', 'unitID'};
end
tbl = movevars(tbl, varOrder, 'Before', 1);

% Assert category order
tbl.genotype = reordercats(tbl.genotype, cfg.lbl.grp);
tbl.day = reordercats(tbl.day, cfg.lbl.day);
tbl.unitID = categorical(tbl.unitID);
if any(contains(tblVars, "unitType"))
    tbl.unitType = reordercats(tbl.unitType, cfg.lbl.unit);
end

if flgClean
    % Remove bad and FS units
    tbl(tbl.unitType == 'Other', :) = [];
    tbl(tbl.unitType == 'FS', :) = [];
    tbl.unitType = removecats(tbl.unitType, {'Other', 'FS'});
    tbl.unitType = [];

    % Remove bac on, bac off, and washout
    if any(tbl.day == 'BAC_OFF')
        tbl(tbl.day == 'BAC_ON', :) = [];
        tbl(tbl.day == 'BAC_OFF', :) = [];
        tbl(tbl.day == 'WASH', :) = [];
        tbl.day = removecats(tbl.day, {'BAC_ON', 'BAC_OFF', 'WASH'});
    end
end

end     % EOF
