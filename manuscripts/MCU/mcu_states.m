% MCU_STATES  Time spent in each vigilance state per mouse (baseline).
%
% A scratch pad, not a function. Run the section, then copy each Prism block.
%
% Baseline sessions across the three genotypes (same set as mcu_ed): one point
% per mouse per state. AccuSleep scores six states; for reporting they are
% collapsed by MERGE below - the one knob in this file. Two numbers per state:
% absolute time (hours) and percent of the recording.


%% ========================================================================
%  TIME IN STATE  (baseline, one point per mouse x state)
%  ========================================================================
% ss.bouts already holds, per raw state and over the whole recording, both the
% total scored seconds (totDur) and that state's share of the recording
% (prctDur, one shared denominator). Merging is therefore a grouped SUM -
% additive, unlike a rate - so a merged state's time is the summed time of its
% parts and its percent the summed percent. This is why evt_stateMerge, which
% merges an event's state LABEL for a rate, is not used here.
%
% MERGE: fieldname = reported state, value = the raw AccuSleep states summed
% into it. This is the only thing to edit. Raw names are WAKE QWAKE LSLEEP NREM
% REM N/REM (v(1).ss.info.names). Examples:
%   - keep quiet wake apart : give it its own field, merge.QWAKE = {'QWAKE'},
%                             and drop 'QWAKE' from merge.WAKE.
%   - pool the small states : merge.Other = {'QWAKE', 'LSLEEP', 'N/REM'}.
% A raw state named in no field simply does not appear, and the percentages sum
% to under 100 by that much - as they already do for unscored / artifact time.

% Files - baseline, three genotypes (as in mcu_ed)
basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl'), ...
    mcu_basepaths('ra')];

% Load scored states
v = basepaths2vars('basepaths', basepaths, 'vars', {'sleep_states'});

% Reporting states and their constituents
merge = struct();
merge.WAKE = {'WAKE', 'QWAKE'};
merge.NREM = {'NREM', 'LSLEEP'};
merge.REM  = {'REM', 'N/REM'};

% Guard: a mistyped raw name would silently drop that state's time
rawNames  = v(1).ss.info.names(1 : size(v(1).ss.bouts.totDur, 2));
mergeVals = struct2cell(merge);
tokens    = [mergeVals{:}];
bad       = tokens(~ismember(tokens, rawNames));
if ~isempty(bad)
    warning('mcu_states:badState', 'merge names not in ss.info.names: %s', ...
        strjoin(unique(bad), ', '));
end

% Build tidy table: one row per session x reported state
rptStates = fieldnames(merge);
mnames    = get_mname(basepaths);
rows      = cell(numel(basepaths), 1);
for iFile = 1 : numel(basepaths)

    bouts  = v(iFile).ss.bouts;
    nState = size(bouts.totDur, 2);
    names  = v(iFile).ss.info.names(1 : nState);        % drop trailing BIN
    totDur = bouts.totDur(1, :);                         % [s] per raw state
    prcDur = bouts.prctDur(1, :);                        % [%] of recording

    durAbs = zeros(numel(rptStates), 1);
    durPct = zeros(numel(rptStates), 1);
    for iRpt = 1 : numel(rptStates)
        idx = ismember(names, merge.(rptStates{iRpt}));
        durAbs(iRpt) = sum(totDur(idx)) / 3600;          % -> hours (/60=min)
        durPct(iRpt) = sum(prcDur(idx));                 % of recording
    end

    rows{iFile} = table( ...
        repmat(categorical(mnames(iFile)), numel(rptStates), 1), ...
        categorical(rptStates, rptStates), durAbs, durPct, ...
        'VariableNames', {'sbjID', 'state', 'durAbs', 'durPct'});
end
tbl = vertcat(rows{:});
tbl.genotype = mcu_geno(tbl.sbjID);

% Quick look (optional triage)
guiTbl_bar(tbl, 'xVar', 'state', 'yVar', 'durPct', 'grpVar', 'genotype');

% -> PRISM, grouped layout: row per state, column per genotype, one subcolumn
% per mouse. flgSort false keeps the category order, so columns come out
% Control / MCU-KO / CAG-MCU-KO and rows in the merge order above. Copy one,
% paste into Prism, then run the other.
tbl2prism(tbl, 'yVar', 'durAbs', 'grpVar', 'genotype', ...
    'rowVar', 'state', 'flgSort', false);       % absolute [h]

tbl2prism(tbl, 'yVar', 'durPct', 'grpVar', 'genotype', ...
    'rowVar', 'state', 'flgSort', false);       % percent of recording

% EOF
