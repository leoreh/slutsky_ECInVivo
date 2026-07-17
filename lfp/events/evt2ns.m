function evt2ns(evtTimes, peakTimes, varargin)
% EVT2NS Export events to a NeuroScope event file.
%
%   EVT2NS(evtTimes, peakTimes, varargin)
%
%   SUMMARY:
%       Writes detected events to the NeuroScope .evt format for external
%       visualization, modality-agnostic (ripples, ED, ...). Creates one file:
%           <basename>.<fileTag>.evt
%       NeuroScope loads it via File -> Load Event File(s) and draws every mark
%       as a dashed line spanning the display; Events -> Next / Previous Event
%       then steps the window between marks. Each event contributes three marks
%       (start, peak, stop), each carrying a description string.
%
%       The palette groups by description, so each distinct string gets its own
%       colour and its own toggle. Passing 'accepted' therefore splits the two
%       curation classes apart: rejected events are labelled 'Rejected ...' and
%       can be shown or hidden on their own, which is what makes the export
%       useful for judging the detector rather than only its survivors.
%
%       Format (NeuroScope user manual, ch. 4): one line per mark,
%       '<ms>\t<description>', ascending in time. The file annotates the .lfp /
%       .dat beside it and shares their time base, so the inputs must be
%       absolute (session) seconds, not window-relative.
%
%   INPUTS:
%       evtTimes    - (Mat)  [N x 2] Event start / stop (s, absolute).
%       peakTimes   - (Vec)  [N x 1] Event peak (s, absolute).
%       varargin    - Parameter/Value pairs:
%           'basepath' - (Char) Target directory. {pwd}
%           'basename' - (Char) File stem (may differ from the folder).
%           'fileTag'  - (Char) NeuroScope id, shown in the events palette.
%                        Exactly three characters and not all digits, which
%                        NeuroScope enforces on load (neuroscopedoc.cpp:
%                        'if(name.length() != 3 ...) return INCORRECT_FILE').
%                        A 4-letter modality tag such as 'ripp' is therefore
%                        rejected outright; ripples pass 'rip'. {'evt'}
%           'lbl'      - (Char) Description stem for accepted events. {'Event'}
%           'accepted' - (Vec)  [N x 1] Curation mask. Rejected events take the
%                        'Rejected' stem instead. Empty labels all N with lbl.
%                        {[]}
%
%   OUTPUTS:
%       None (writes <basename>.<fileTag>.evt).
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 260706 (generalised from ripp2ns; fileTag param; no cd).
%       Updated: 260715 (.evt rather than .res / .clu - no sampling rate and no
%                spike-group binding to get wrong; seconds in; acceptance
%                labels; wired into ripp_wrapper).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'evtTimes', @isnumeric);
addRequired(p, 'peakTimes', @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'fileTag', 'evt', @(x) ischar(x) && numel(x) == 3);
addParameter(p, 'lbl', 'Event', @ischar);
addParameter(p, 'accepted', [], @(x) isempty(x) || islogical(x));
parse(p, evtTimes, peakTimes, varargin{:});
basepath = p.Results.basepath;
basename = p.Results.basename;
fileTag  = p.Results.fileTag;
lbl      = p.Results.lbl;
accepted = p.Results.accepted;

if isempty(basename)
    [~, basename] = fileparts(basepath);
end
evtfile = fullfile(basepath, [basename, '.', fileTag, '.evt']);

nEvents = numel(peakTimes);
if isempty(accepted)
    accepted = true(nEvents, 1);
end
accepted = accepted(:);

%% ========================================================================
%  BUILD MARKS
%  ========================================================================
% Stack start / peak / stop onto one timeline (ms). The description doubles as
% the palette key: the curation class picks the stem, the mark picks the
% suffix, and the two are pasted per mark.
marks = [evtTimes(:, 1); peakTimes(:); evtTimes(:, 2)] * 1000;

lblStem = repmat({lbl}, nEvents, 1);
lblStem(~accepted) = {'Rejected'};
lblMark = [repmat({'start'}, nEvents, 1); ...
           repmat({'peak'},  nEvents, 1); ...
           repmat({'stop'},  nEvents, 1)];
desc = strcat(repmat(lblStem, 3, 1), {' '}, lblMark);

% Ascending time; the palette navigation walks the file in order.
[marks, idxSort] = sort(marks);
desc = desc(idxSort);

%% ========================================================================
%  WRITE FILE
%  ========================================================================
fid = fopen(evtfile, 'w');
if fid == -1, error('evt2ns:open', 'failed to open %s', evtfile); end
for iMark = 1 : numel(marks)
    fprintf(fid, '%.3f\t%s\n', marks(iMark), desc{iMark});
end
fclose(fid);

end     % EOF
