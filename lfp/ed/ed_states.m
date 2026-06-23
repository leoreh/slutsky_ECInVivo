function [stateIdx, edStates] = ed_states(ed, boutTimes, varargin)
% ED_STATES Assign a vigilance state to each ED and tabulate per-bout rate.
%
%   [stateIdx, edStates] = ED_STATES(ed, boutTimes, varargin)
%
%   SUMMARY:
%       Thin forwarder to ripp_states, whose math is event-agnostic: it maps
%       each event peak into the supplied sleep bouts and builds a per-bout
%       table of rate / density / duration / state. Run quietly and re-saved
%       under the ed naming convention.
%
%   INPUTS:
%       ed          - (Struct) Requires .times [N x 2] and .peakTime [N x 1] (s).
%       boutTimes   - (Cell) {nStates x 1} of [start end] bout matrices.
%       varargin    - Parameter/Value pairs:
%           'basepath' - (Char) Save location. {pwd}
%           'flgSave'  - (Log)  Save <basename>.edStates.mat? {false}
%
%   OUTPUTS:
%       stateIdx    - (Cat)   [N x 1] vigilance state per event.
%       edStates    - (Table) Per-bout [Rate, Density, Duration, State, Start, End].
%
%   DEPENDENCIES:
%       ripp_states.
%
%   HISTORY:
%       Created: 22 Jun 2026

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'ed', @isstruct);
addRequired(p, 'boutTimes', @iscell);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', false, @islogical);

parse(p, ed, boutTimes, varargin{:});
ed        = p.Results.ed;
boutTimes = p.Results.boutTimes;
basepath  = p.Results.basepath;
flgSave   = p.Results.flgSave;

[~, basename] = fileparts(basepath);

%% ========================================================================
%  FORWARD TO RIPP_STATES
%  ========================================================================

[stateIdx, edStates] = ripp_states(ed.times, ed.peakTime, boutTimes, ...
    'basepath', basepath, 'flgPlot', false, 'flgSave', false);

if flgSave
    save(fullfile(basepath, [basename, '.edStates.mat']), 'edStates', '-v7.3');
end

end     % EOF
