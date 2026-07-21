function [ed, hFig] = ed_curate(basepath, varargin)
% ED_CURATE Curate discharges by waveform TYPE (stage 2).
%
%   [ed, hFig] = ED_CURATE(basepath, varargin)
%
%   SUMMARY:
%       Loads <basename>.ed.mat and its per-event waveforms, applies met.qa as
%       a POOL, groups the pool into waveform clusters, and lets you accept
%       whole clusters. The GUI itself is evt_curate, shared with the ripple
%       pipeline; everything here is loading.
%
%       A 24 h recording proposes thousands of candidates and holds a few dozen
%       discharges, so the unit of curation is a TYPE, not an event. Read
%       evt_curate for what the controls do.
%
%       Headless (flgGui = false) applies only met.qa, for a batch run that has
%       no human. That mask is NOT an answer - it is the pool.
%
%   INPUTS:
%       basepath - <char> session directory (must hold <basename>.ed.mat).
%       varargin - Parameter/Value:
%           'basename' - <char>   file stem. {folder name}
%           'met'      - <struct> config; reads .qa and .clust.
%                                 {ed_methods('default')}
%           'flgGui'   - <log>    open the GUI (true) or filter headless.{true}
%           'Visible'  - <char>   'on' | 'off', for headless GUI tests. {'on'}
%           'verbose'  - <log>    print progress? {true}
%
%   OUTPUTS:
%       ed   - <struct> the loaded events, with .accepted as of the call.
%       hFig - <handle> the GUI figure ([] when headless).
%
%   DEPENDENCIES:
%       evt_files, evt_curate, ed_methods.
%
%   HISTORY:
%       260720 created as the ED twin of ripp_curate (threshold knobs over a
%              kept-vs-removed mean waveform).
%       260721 rebuilt on waveform clustering: accept TYPES, not events.
%       260721b shape rejections became sticky and state rejections reversible.
%       260722 the GUI moved to lfp/events/evt_curate and is now shared with the
%              ripple pipeline; this file is the ED loader over it. Two things
%              changed in the move: the vigilance-state scope is recorded in
%              info.qa (info.clustStates is still READ, for files saved before
%              this), and the filtered trace joins the raw one as a second Y
%              option in the view.

%% ========================================================================
%  ARGUMENTS + LOAD
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'met', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'flgGui', true, @islogical);
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(char(x), {'on', 'off'})));
addParameter(p, 'verbose', true, @islogical);
parse(p, basepath, varargin{:});
met     = p.Results.met;
flgGui  = p.Results.flgGui;
verbose = p.Results.verbose;

basename = p.Results.basename;
if isempty(basename), [~, basename] = fileparts(basepath); end
if isempty(met), met = ed_methods('default'); end

files = evt_files(basepath, basename, 'ed');
if ~isfile(files.evt)
    error('ed_curate:noFile', ...
        '%s not found; run detection first (ed_wrapper).', files.evt);
end
S  = load(files.evt, 'ed');
ed = S.ed;

cfg = struct('met', met, 'file', files.evt, 'var', 'ed', ...
    'basepath', basepath, 'basename', basename, 'lbl', 'ED', ...
    'flgGui', flgGui, 'Visible', char(p.Results.Visible), 'onSaved', []);

%% ========================================================================
%  CURATE
%  ========================================================================
% Waveforms are the clustering's whole input, so they are loaded only when
% there is a human to look at them.
wv = struct();
tst = [];
if flgGui
    [wv, tst] = loadMaps(files.maps, ed);
end

[ed.accepted, hFig] = evt_curate(ed, wv, tst, cfg);

if verbose
    how = 'headless';
    if flgGui, how = 'GUI'; end
    fprintf('[ED_CURATE] %s : %d / %d accepted (%s)\n', basename, ...
        nnz(ed.accepted), numel(ed.accepted), how);
end

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function [wv, tst] = loadMaps(file, ed)
% Per-event waveforms behind the clustering and the view. Passed whole: the
% flanks beyond the clustering window are where a discharge separates from a
% sharp wave, and the view opens on met.clust.win with the rest a zoom away.
% Unlike the ripple loader this does not detrend - evt_clust detrends its own
% window, and nothing here crops away the flanks it needs.
if ~isfile(file)
    error('ed_curate:noMaps', ...
        'no edMaps file; re-run detection with flgSave.');
end
S = load(file, 'edMaps');
if size(S.edMaps.lfp, 1) ~= numel(ed.peakTime)
    error('ed_curate:staleMaps', ...
        'edMaps does not match the event list; re-run detection.');
end

tst = S.edMaps.tstamps;
wv = struct('lfp', double(S.edMaps.lfp));
if isfield(S.edMaps, 'filt')
    wv.filt = double(S.edMaps.filt);
end

end     % loadMaps
