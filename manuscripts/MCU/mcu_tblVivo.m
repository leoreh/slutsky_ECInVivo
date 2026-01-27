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
%  LOAD DATA
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

if ismember('brst', presets)
    vars = [vars, 'brstStats'];
    varMap.bRate     = 'stats.rate';
    varMap.bDur      = 'stats.dur';
    varMap.bFreq     = 'stats.freq';
    varMap.bIBI      = 'stats.ibi';
    varMap.pBspk     = 'stats.pBspk';
    varMap.nBspk     = 'stats.nBspk';
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
    varMap.frRipp = ['rippSpks.', spkType, '.frRipp'];
    varMap.frRand = ['rippSpks.', spkType, '.frRand'];
    varMap.frZ = ['rippSpks.', spkType, '.frZ'];
    varMap.frMod = ['rippSpks.', spkType, '.frMod'];
    varMap.pFire = ['rippSpks.', spkType, '.pFire'];
    varMap.cRipp = ['rippSpks.', spkType, '.cRipp'];
    varMap.rankMean = ['rippSpks.', spkType, '.rankMean'];
    varMap.rankVar = ['rippSpks.', spkType, '.rankVar'];
    varMap.frActive = ['rippSpks.', spkType, '.frActive'];
    varMap.cActive = ['rippSpks.', spkType, '.cActive'];
    varMap.frPre = ['rippSpks.', spkType, '.frPre'];
    varMap.frPost = ['rippSpks.', spkType, '.frPost'];
    varMap.asym = ['rippSpks.', spkType, '.asym'];
    varMap.com = ['rippSpks.', spkType, '.com'];

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
end

if ismember('rippStates', presets)
    vars = {'rippStates'};
    varMap = struct();
    varMap.Rate         = 'rippStates.Rate';
    varMap.Density      = 'rippStates.Density';
    varMap.Duration     = 'rippStates.Duration';
    varMap.State        = 'rippStates.State';
end

if ismember('rippPeth', presets)
    vars = {'units', 'rippPeth'};
    varMap = struct();
    varMap.UnitType = 'units.type';
    varMap.pethSU = 'rippPeth.pethSU';
end

% Load
if isempty(basepaths)
    basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
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

% Post-process ripp peth
if ismember('rippPeth', presets) 
    for iFile = 1:length(v)
        v(iFile).rippPeth.pethSU = squeeze(mean(v(iFile).rippPeth.su.ripp, 2, 'omitnan'));
    end
    xVec = v(iFile).rippPeth.su.tstamps;
end


% Metadata
tagFiles = struct();
tagFiles.Name = get_mname(basepaths);
[~, fileNames] = fileparts(basepaths);
tagFiles.File = fileNames;

% Table
tbl = v2tbl('v', v, 'varMap', varMap, 'tagAll',...
    struct(), 'tagFiles', tagFiles, 'idxCol', []);

%% ========================================================================
%  PROCESS METADATA
%  ========================================================================

% Group metadata
tbl.Group = ones(height(tbl), 1) * 1;
tbl.Group(ismember(tbl.Name, cfg.miceMCU), :) = 2;
tbl.Group = categorical(tbl.Group, [1, 2], cfg.lbl.grp);

% Day metadata
fileTbl = unique(tbl(:, {'Name', 'File'}), 'rows');
fileGrp = findgroups(fileTbl.Name);
dayCell = splitapply(@(x) {(1:numel(x))'}, fileTbl.File, fileGrp);
fileTbl.Day = vertcat(dayCell{:});
tbl = join(tbl, fileTbl, 'Keys', {'Name', 'File'});
tbl.Day = categorical(tbl.Day, [1 : 7], cfg.lbl.day);

% Reorder columns
tblVars = tbl.Properties.VariableNames;
if any(contains(tblVars, "UnitType"))
    varOrder = {'Group', 'Name', 'File', 'Day', 'UnitID', 'UnitType'};
else
    varOrder = {'Group', 'Name', 'File', 'Day', 'UnitID'};
end
tbl = movevars(tbl, varOrder, 'Before', 1);

% Assert category order
tbl.Group = reordercats(tbl.Group, cfg.lbl.grp);
tbl.Day = reordercats(tbl.Day, cfg.lbl.day);
tbl.UnitID = categorical(tbl.UnitID);
if any(contains(tblVars, "UnitType"))
    tbl.UnitType = reordercats(tbl.UnitType, cfg.lbl.unit);
end

if flgClean
    % Remove bad units
    tbl(tbl.UnitType == 'Other', :) = [];
    tbl.UnitType = removecats(tbl.UnitType, 'Other');

    % Remove bac on / off
    tbl(tbl.Day == 'BAC_ON', :) = [];
    tbl(tbl.Day == 'BAC_OFF', :) = [];
    tbl.Day = removecats(tbl.Day, {'BAC_ON', 'BAC_OFF'});
end

end     % EOF
