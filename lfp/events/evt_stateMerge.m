function s = evt_stateMerge(s)
% EVT_STATEMERGE Collapse AccuSleep labels to the states an event report uses.
%
%   s = EVT_STATEMERGE(s)
%
%   SUMMARY:
%       QWAKE joins WAKE, LSLEEP joins NREM. Everything else passes through.
%
%       This is a REPORTING grouping applied after detection, not a change to
%       the scoring: <basename>.sleep_states.mat keeps all six labels, and
%       anything asking a question about light sleep still has them. It exists
%       because an event count is a small number - a mouse carries dozens of
%       discharges over 24 h - and splitting them across six states leaves
%       several with one or two events and an exposure too short to divide by.
%       QWAKE and LSLEEP are also the two labels the scorer is least certain
%       about, being the transitions either side of the states they merge into.
%
%       Merging LABELS is the easy half. A merged RATE is the summed count over
%       the summed exposure, never the mean of the two rates - see ed_tbl,
%       which aggregates counts and durations before it divides.
%
%   INPUT:
%       s - <cat|str> state labels, as a categorical or a string/cellstr.
%
%   OUTPUT:
%       s - same type, with the merged labels and no empty categories.
%
%   HISTORY:
%       260721 created so the map lives in one place rather than in ed_tbl and
%              ed_wvTbl separately.

PAIRS = {'QWAKE', 'WAKE'; 'LSLEEP', 'NREM'};

flgCat = iscategorical(s);
if ~flgCat
    was = s;
    s = categorical(cellstr(string(s)));
end

for iPair = 1 : size(PAIRS, 1)
    src = PAIRS{iPair, 1};
    dst = PAIRS{iPair, 2};
    cats = categories(s);
    if ~ismember(src, cats)
        continue
    end
    if ismember(dst, cats)
        s = mergecats(s, {src, dst}, dst);
    else
        s = renamecats(s, src, dst);    % mergecats needs both to exist
    end
end
s = removecats(s);

if ~flgCat
    s = reshape(string(s), size(was));
    if iscellstr(was), s = cellstr(s); end %#ok<ISCLSTR>
end

end     % EOF
