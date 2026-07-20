function [stStates, brstStates] = spk_states(varargin)
% SPK_STATES Computes per-unit spike metrics separately in each state.
%
%   [stStates, brstStates] = SPK_STATES(varargin)
%
%   SUMMARY:
%       Per-session producer for the state-resolved spike tables. Maps the
%       two metric families over the vigilance-state bouts via spk_byCond
%       and writes one file each, so a later analysis loads only the family
%       it needs:
%
%           <basename>.stStates.mat    acg and isi metrics
%           <basename>.brstStates.mat  burst and rate metrics
%
%       Both tables carry the same (uid, state) keys in the same order, so
%       they join on those two columns. Burst DETECTION stays state-agnostic
%       - a burst is a burst, and letting the ISI thresholds float with the
%       state would confound every comparison downstream. Only the summary
%       is conditioned.
%
%       Conditions are whatever is passed. The default is the sleep states,
%       but bouts / lbls also take time chunks (n2chunks) or drug epochs
%       without any change here.
%
%   INPUTS:
%       varargin - Parameter/Value pairs:
%           'basepath' - (Char) Recording path. {pwd}
%           'spktimes' - (Cell) {nUnits x 1} spike times [s]. Loaded from
%                               <basename>.spikes.cellinfo.mat when empty.
%           'burst'    - (Struct) burst_detect output. Loaded from
%                               <basename>.burst.mat when empty.
%           'bouts'    - (Cell) {nCond x 1} of [nBouts x 2] [s]. Taken from
%                               ss.bouts.times when empty.
%           'lbls'     - (Cell) {nCond x 1} condition names. Taken from
%                               ss.info.names when empty.
%           'flgSave'  - (Log)  Save both tables. {true}
%
%   OUTPUTS:
%       stStates   - (Table) [nUnits * nCond x ...] acg / isi metrics.
%       brstStates - (Table) [nUnits * nCond x ...] burst / rate metrics.
%
%   DEPENDENCIES:
%       spk_byCond, spktimes_metrics, burst_stats, backup_file.
%
%   HISTORY:
%       Created: 260719
%
%   See also: SPK_BYCOND, SPKTIMES_METRICS, BURST_STATS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'spktimes', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'burst', [], @(x) isstruct(x) || isempty(x));
addParameter(p, 'bouts', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'lbls', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'flgSave', true, @islogical);

parse(p, varargin{:});
basepath = p.Results.basepath;
spktimes = p.Results.spktimes;
burst    = p.Results.burst;
bouts    = p.Results.bouts;
lbls     = p.Results.lbls;
flgSave  = p.Results.flgSave;

[~, basename] = fileparts(basepath);


%% ========================================================================
%  LOAD
%  ========================================================================

if isempty(spktimes)
    s = load(fullfile(basepath, [basename, '.spikes.cellinfo.mat']), 'spikes');
    spktimes = s.spikes.times;
end

if isempty(burst)
    s = load(fullfile(basepath, [basename, '.burst.mat']), 'burst');
    burst = s.burst;
end

if isempty(bouts)
    s = load(fullfile(basepath, [basename, '.sleep_states.mat']), 'ss');
    bouts = s.ss.bouts.times;

    % accusleep names carry a trailing label for unclassified bins ('BIN')
    % that has no bout array, so names always runs one longer than bouts
    lbls = s.ss.info.names(1 : numel(bouts));
end


%% ========================================================================
%  MAP
%  ========================================================================

stStates = spk_byCond(@(b) spktimes_metrics(spktimes, b), bouts, lbls);

brstStates = spk_byCond(@(b) burst_stats(burst, spktimes, 'winCalc', b, ...
    'flgPool', true), bouts, lbls);


%% ========================================================================
%  SAVE
%  ========================================================================

if flgSave
    stFile = fullfile(basepath, [basename, '.stStates.mat']);
    brstFile = fullfile(basepath, [basename, '.brstStates.mat']);

    backup_file(stFile);
    backup_file(brstFile);

    save(stFile, 'stStates')
    save(brstFile, 'brstStates')
end

end     % EOF
