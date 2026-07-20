function [sig, fs, edCh] = ed_sigLoad(basepath, varargin)
% ED_SIGLOAD Load the ED detection signal: one raw channel from the .lfp.
%
%   [sig, fs, edCh] = ED_SIGLOAD(basepath, varargin)
%
%   SUMMARY:
%       One source, one channel. The channel comes from ed_pickCh unless the
%       caller names it.
%
%       This used to read the sleep_sig eeg instead, and that was wrong in a
%       way worth recording: sleep_sig.eeg is the MEAN of whichever channels
%       were picked for sleep scoring, which is a different set in every mouse
%       - in raMCU1 an average across two shanks, in raMCU2 one that included
%       the weakest channel of fifteen, and low-passed at 450 Hz in four mice
%       but unfiltered in the fifth. A discharge is a laminar event, so
%       averaging across sites attenuates it by a per-mouse amount that has
%       nothing to do with the biology, and the pipeline was comparing
%       measurements rather than mice. See dev/ed_pipeline_rebuild.md.
%
%       Requires a binary <basename>.lfp and a session.mat. Every session in
%       the MCU cohort has both.
%
%   INPUTS:
%       basepath - <char> session directory.
%       varargin - Parameter/Value:
%           'basename' - <char>   file stem. {folder name}
%           'edCh'     - <num>    explicit channel; else ed_pickCh. {[]}
%           'win'      - <vec>    window [start end] (s). {[0 Inf]}
%           'session'  - <struct> session metadata. {loaded if empty}
%           'passband' - <vec>    band ed_pickCh scores in (Hz). {[60 150]}
%
%   OUTPUTS:
%       sig      - <vec>  [n x 1] windowed signal (microvolts).
%       fs       - <num>  sampling frequency [Hz].
%       edCh     - <num>  the channel used (1-indexed).
%
%   DEPENDENCIES:
%       basepaths2vars, ed_pickCh, evt_loadCh.
%
%   HISTORY:
%       Created: 260622
%       Updated: 260706 (lfp branch via the shared evt_loadCh).
%       Updated: 260721 (rebuilt on the .lfp alone; met.chMode and the
%                sleep_sig eeg branch are gone, and so is the EMG output -
%                once the three shape criteria were applied the EMG metric
%                removed nothing measurable, so nothing read it.)

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'edCh', [], @isnumeric);
addParameter(p, 'win', [0 Inf], @isnumeric);
addParameter(p, 'session', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'passband', [60 150], @isnumeric);
parse(p, basepath, varargin{:});
edCh     = p.Results.edCh;
win      = p.Results.win;
session  = p.Results.session;
passband = p.Results.passband;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end

%% ========================================================================
%  CHANNEL + SIGNAL
%  ========================================================================
if isempty(session)
    v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'}, ...
        'flgPrnt', false);
    session = v.session;
end

if isempty(edCh)
    edCh = ed_pickCh(basepath, 'basename', basename, 'session', session, ...
        'win', win, 'passband', passband);
end

[sig, fs] = evt_loadCh(basepath, basename, session, edCh, win, []);

end     % EOF
