function tests = test_rippCurate
% TEST_RIPPCURATE Unit tests for ripp_curate over the shared evt_curate GUI.
%
%   tests = TEST_RIPPCURATE
%
%   SUMMARY:
%       Builds a synthetic ripple session in a temp folder - oscillations,
%       step artifacts and background - so the cluster curation can be
%       exercised without touching a real recording or a curated mask.
%
%       The fixture is rebuilt PER TEST, in its own directory: ripp_curate
%       resumes a saved curation, so a test that presses Save would otherwise
%       decide what the next test opens into.
%
%       What is tested here and not in test_edCurate: the ripple side of the
%       shared engine - the generated metric knobs, the two cost caps (nFit,
%       nView), the ripple save contract (.spks stripped, products invalidated)
%       and the two waveforms in the view.
%
%       Run with: runtests('test_rippCurate')
%
%   HISTORY:
%       260722 created with the ripp_curate cluster rebuild.

tests = functiontests(localfunctions);
end


%% ========================================================================
%  FIXTURE
%  ========================================================================
function setup(tc)
rng(11);
nA = 120; nB = 90; nN = 90;         % ripple, step artifact, background
fs = 1250;
tst = -0.1 : 1 / fs : 0.1;
nS = numel(tst);

% family A: a 150 Hz burst on a sharp wave. family B: a step deflection with
% no oscillation - the population curation exists to remove.
wvA = sin(2 * pi * 150 * tst) .* exp(-(tst / 0.015) .^ 2) * 120 ...
    - 180 * exp(-(tst / 0.030) .^ 2);
wvB = 400 * tanh(tst / 0.002);

% background is band-limited, not white: white noise carries more variance
% than either real shape and would take the principal components with it
bg = filtfilt(ones(1, 25) / 25, 1, 200 * randn(nN, nS)')' * 5;

wv = [repmat(wvA, nA, 1) .* (0.8 + 0.4 * rand(nA, 1)); ...
    repmat(wvB, nB, 1) .* (0.8 + 0.4 * rand(nB, 1)); bg];
wv = wv + filtfilt(ones(1, 9) / 9, 1, 30 * randn(size(wv))')' * 3;
nEv = size(wv, 1);

% a crude band-pass, only so the view has a second Y and the clustering has a
% sibling field to ignore
filt = wv - filtfilt(ones(1, 31) / 31, 1, wv')';

ripp = struct();
ripp.peakTime = sort(rand(nEv, 1) * 3600);
ripp.times    = [ripp.peakTime - 0.02, ripp.peakTime + 0.02];
ripp.amp      = max(abs(wv), [], 2);
ripp.dur      = [40 + randn(nA, 1) * 5; 25 + randn(nB, 1) * 5; ...
    30 + randn(nN, 1) * 8];
ripp.freq     = [150 + randn(nA, 1) * 8; 100 + randn(nB, 1) * 20; ...
    110 + randn(nN, 1) * 25];
ripp.freqPeak = ripp.freq;
ripp.peakProm = [3 + randn(nA, 1) * 0.5; randn(nB, 1) * 0.5; ...
    0.5 + randn(nN, 1) * 0.5];
ripp.emg      = abs(randn(nEv, 1)) * 0.2;         % all inside [-Inf 1]
ripp.spkGain  = [3 + randn(nA, 1); 1.5 + randn(nB, 1) * 0.3; ...
    1.4 + randn(nN, 1) * 0.3];
ripp.state    = categorical(repmat({'NREM'}, nEv, 1));
ripp.accepted = true(nEv, 1);
ripp.info = struct('fs', fs, 'win', [0 Inf], 'rippCh', 1, ...
    'passband', [80 250], 'detectMet', 3, 'zMet', 'nrem', 'otlThr', 8, ...
    'met', 'test');

rippMaps = struct('tstamps', tst, 'lfp', single(wv), 'filt', single(filt));

tc.TestData.dir = tempname;
mkdir(tc.TestData.dir);
tc.TestData.name = 'ripptest';
save(fullfile(tc.TestData.dir, 'ripptest.ripp.mat'), 'ripp', '-v7.3');
save(fullfile(tc.TestData.dir, 'ripptest.rippMaps.mat'), 'rippMaps', '-v7.3');

tc.TestData.wv = wv;
tc.TestData.tst = tst;
tc.TestData.truth = [ones(nA, 1); 2 * ones(nB, 1); 3 * ones(nN, 1)];
tc.TestData.nEv = nEv;
end


function teardown(tc)
close all force
try
    rmdir(tc.TestData.dir, 's');    % a lingering matfile lock is not a failure
catch
end
end


%% ========================================================================
%  EVT_CLUST - the ripple window and the fit cap
%  ========================================================================
function test_clustSeparatesRippleFromStep(tc)
% The oscillation and the step must not land in one cluster - that separation
% is the whole reason curation moved to shapes.
met = ripp_methods('default');
cid = evt_clust(tc.TestData.wv, tc.TestData.tst, 'win', met.clust.win, ...
    'nClust', 6, 'wSize', met.clust.wSize);

truth = tc.TestData.truth;
tc.verifyNotEqual(mode(cid(truth == 1)), mode(cid(truth == 2)), ...
    'ripples and step artifacts collapsed into one cluster');
tc.verifyGreaterThan(mean(cid(truth == 1) == mode(cid(truth == 1))), 0.7);
tc.verifyGreaterThan(mean(cid(truth == 2) == mode(cid(truth == 2))), 0.7);
end


function test_clustNFitLabelsEveryEvent(tc)
% The fit cap must bound the fit, not the labelling: every event still gets a
% cluster, because the mask is built from labels and an unlabelled event can
% never be accepted.
[cid, cInfo] = evt_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 5, ...
    'nFit', 100);

tc.verifyEqual(cInfo.nFit, 100, 'the fit cap was ignored');
tc.verifyEqual(numel(cid), tc.TestData.nEv);
tc.verifyFalse(any(isnan(cid)), 'projected events lost their labels');
tc.verifyEqual(sort(unique(cid))', 1 : cInfo.nClust);
end


function test_clustNFitIsDeterministic(tc)
% The subsample is drawn from evt_clust's own seeded stream, so a Re-cluster
% with nothing changed is a no-op rather than a reshuffle - and the caller's
% RNG must come back untouched.
a = evt_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 5, 'nFit', 120);
b = evt_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 5, 'nFit', 120);
tc.verifyEqual(a, b, 'a subsampled fit is not reproducible');

rng(42); r1 = rand();
rng(42); evt_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 4, 'nFit', 120);
r2 = rand();
tc.verifyEqual(r1, r2, 'evt_clust leaked its RNG seed to the caller');
end


function test_clustNFitBelowFloorIsNoCap(tc)
% A cap smaller than the fit floor is a misconfiguration, not a request for a
% 20-event fit. It must be read as no cap rather than silently returning a
% partition estimated from 20 points.
ref = evt_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 5);
for v = [0 -1 5]
    [cid, cInfo] = evt_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 5, ...
        'nFit', v);
    tc.verifyEqual(cInfo.nFit, tc.TestData.nEv, ...
        sprintf('nFit = %g was taken as a cap', v));
    tc.verifyEqual(cid, ref);
end
end


function test_clustNFitOffIsUnchanged(tc)
% An empty cap must leave the partition exactly as it was before the cap
% existed - that is what keeps the ED pipeline bit-identical.
a = evt_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 5);
b = evt_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 5, 'nFit', []);
c = evt_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 5, 'nFit', 1e6);
tc.verifyEqual(b, a);
tc.verifyEqual(c, a, 'a cap above the pool size changed the partition');
end


%% ========================================================================
%  HEADLESS
%  ========================================================================
function test_headlessAppliesGate(tc)
% Headless applies met.qa - the full spec, states included - and writes it.
met = ripp_methods('default');
ripp = ripp_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'met', met, 'flgGui', false, 'verbose', false);

S = load(fullfile(tc.TestData.dir, 'ripptest.ripp.mat'), 'ripp');
tc.verifyEqual(S.ripp.accepted, evt_gate(ripp, met.qa));
tc.verifyFalse(isfield(S.ripp, 'clustId'), ...
    'the headless gate invented cluster labels');
end


function test_headlessKeepsSavedPartition(tc)
% A batch gate makes no shape judgement, so it must leave one alone. Nulling
% the partition would let one wrapper run destroy a session's manual curation
% with no warning and no way back.
[~, h1] = openGui(tc);
st = h1.UserData;
st.chk(2).Value = false;
st.chk(2).ValueChangedFcn([], []);
press(h1, 'Save');
S1 = load(fullfile(tc.TestData.dir, 'ripptest.ripp.mat'), 'ripp');
close(h1, 'force');

ripp_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', false, 'verbose', false);

S2 = load(fullfile(tc.TestData.dir, 'ripptest.ripp.mat'), 'ripp');
tc.verifyEqual(S2.ripp.clustId, S1.ripp.clustId, 'the partition was nulled');
tc.verifyEqual(S2.ripp.info.clustSel, S1.ripp.info.clustSel, ...
    'the cluster choice was nulled');
end


function test_headlessLoadsNoWaveforms(tc)
% A batch gate must not pay for the map read. Deleting the maps file must not
% stop it.
delete(fullfile(tc.TestData.dir, 'ripptest.rippMaps.mat'));
met = ripp_methods('default');
ripp = ripp_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'met', met, 'flgGui', false, 'verbose', false);
tc.verifyEqual(nnz(ripp.accepted), nnz(evt_gate(ripp, met.qa)));
end


%% ========================================================================
%  GUI
%  ========================================================================
function test_guiOpensAndClusters(tc)
% The GUI must build, cluster, and expose one checkbox per cluster.
[~, hFig] = openGui(tc);
st = hFig.UserData;

tc.verifyGreaterThan(st.nClust, 1);
tc.verifyEqual(numel(st.chk), st.nClust);
tc.verifyEqual(nnz(~isnan(st.cid)), nnz(st.pool));
tc.verifyEqual(numel(st.cid), tc.TestData.nEv);
tc.verifyTrue(all(arrayfun(@(h) h.Value, st.chk)), ...
    'clusters do not start accepted');
end


function test_knobsComeFromQaRanges(tc)
% One knob per FINITE bound in met.qa.ranges, labelled with its metric. The
% GUI never names a metric in code, so adding one to ripp_methods is enough.
[~, hFig] = openGui(tc);
st = hFig.UserData;

tc.verifyEqual(sort({st.knob.fld}), {'emg', 'spkGain'});
tc.verifyEqual(st.knob(strcmp({st.knob.fld}, 'emg')).iBnd, 2, ...
    'the emg knob is not on its upper bound');
tc.verifyEqual(st.knob(strcmp({st.knob.fld}, 'spkGain')).iBnd, 1);
tc.verifyEqual(knobOf(hFig, 'emg').Value, 1);
tc.verifyEqual(knobOf(hFig, 'spkGain').Value, 1);
end


function test_knobShrinksPool(tc)
% Raising a knob must shrink the pool on Re-cluster. The knob alone only
% updates the label.
[~, hFig] = openGui(tc);
n0 = nnz(hFig.UserData.pool);

hKnob = knobOf(hFig, 'spkGain');
hKnob.Value = 2.5;
press(hFig, 'Re-cluster');

tc.verifyLessThan(nnz(hFig.UserData.pool), n0);
tc.verifyEqual(nnz(~isnan(hFig.UserData.pool)), tc.TestData.nEv);
tc.verifyTrue(all(hFig.UserData.evt.spkGain(hFig.UserData.pool) >= 2.5));
end


function test_viewCarriesBothWaveforms(tc)
% The raw and the band-passed trace must both be selectable, opening on the
% raw one - whether a cluster oscillates is usually what settles it.
[~, hFig] = openGui(tc);
ud = hFig.UserData.hPanel.UserData;

tc.verifyTrue(isfield(ud, 'setDataFcn'));
tc.verifyEqual(ud.yVar, 'lfp');
tc.verifyTrue(ismember('lfp', ud.ddYVar.Items));
tc.verifyTrue(ismember('filt', ud.ddYVar.Items), ...
    'the filtered trace is not offered as a Y variable');
tc.verifyTrue(ismember('cluster', ud.ddPlotBy.Items));
tc.verifyTrue(ismember('state', ud.ddPlotBy.Items));
tc.verifyTrue(ismember('status', ud.ddGrpBy.Items));
end


function test_viewCapDoesNotTouchTheMask(tc)
% nView moves ROWS only. The mask, the checkbox counts and what Save writes are
% all computed on every event.
met = ripp_methods('default');
met.clust.nView = 3;
met.clust.nClust = 5;
[~, hFig] = ripp_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'met', met, 'flgGui', true, 'Visible', 'off', 'verbose', false);
tc.addTeardown(@() close(hFig, 'force'));

st = hFig.UserData;
ud = st.hPanel.UserData;
tc.verifyLessThanOrEqual(height(ud.dataTbl), 3 * st.nClust);
tc.verifyLessThan(height(ud.dataTbl), nnz(st.acc), 'the cap did nothing');
tc.verifySubstring(st.lblKeep.Text, 'view sampled', ...
    'the view was thinned without saying so');

press(hFig, 'Save');
S = load(fullfile(tc.TestData.dir, 'ripptest.ripp.mat'), 'ripp');
tc.verifyEqual(nnz(S.ripp.accepted), nnz(st.acc), ...
    'the view cap reached the saved mask');
end


function test_autoCountReachesTheRule(tc)
% 0 in the box asks for evt_clust's rule, not for the shipped count the box
% opened on. The label says "0 = auto"; it has to be true.
[~, hFig] = openGui(tc);
st = hFig.UserData;
tc.verifyEqual(st.edK.Value, 20, 'the box did not open on the shipped count');

st.edK.Value = 0;
press(hFig, 'Re-cluster');
st = hFig.UserData;
tc.verifyEqual(st.nClust, round(0.65 * sqrt(nnz(st.pool))), ...
    'the box fell back to the preset count instead of the rule');
end


function test_rejectAllSurvivesReopen(tc)
% Rejecting every cluster is an answer, not an absence of one. Reopening must
% not silently re-tick everything and hand back the whole pool.
[~, h1] = openGui(tc);
arrayfun(@(h) set(h, 'Value', false), h1.UserData.chk);
h1.UserData.chk(1).ValueChangedFcn([], []);
press(h1, 'Save');
close(h1, 'force');

[~, h2] = openGui(tc);
tc.verifyFalse(any(arrayfun(@(h) h.Value, h2.UserData.chk)), ...
    'a reject-all curation came back as accept-all');
tc.verifyEqual(nnz(h2.UserData.acc), 0);
end


function test_savedSpecMatchesTheMask(tc)
% A knob typed into after a fit only PREVIEWS - the labels, and so the mask,
% still come from the old pool. Save must therefore record the ranges the fit
% used, or info.qa claims a filter the saved mask does not obey.
[~, hFig] = openGui(tc);
hKnob = knobOf(hFig, 'spkGain');
hKnob.Value = 2.5;                      % not followed by Re-cluster
press(hFig, 'Save');

S = load(fullfile(tc.TestData.dir, 'ripptest.ripp.mat'), 'ripp');
tc.verifyEqual(S.ripp.info.qa.ranges.spkGain, [1 Inf], ...
    'the saved spec claims a filter the mask does not obey');
tc.verifyGreaterThanOrEqual(nnz(evt_gate(S.ripp, S.ripp.info.qa)), ...
    nnz(S.ripp.accepted), 'the mask is not a subset of its own spec');
end


function test_saveWritesLabelsAndChoice(tc)
% Ticking one cluster must accept exactly its events and record the choice.
[~, hFig] = openGui(tc);
st = hFig.UserData;
arrayfun(@(h) set(h, 'Value', false), st.chk);
st.chk(1).Value = true;
nExpect = nnz(st.cid == 1);

press(hFig, 'Save');

S = load(fullfile(tc.TestData.dir, 'ripptest.ripp.mat'), 'ripp');
tc.verifyEqual(nnz(S.ripp.accepted), nExpect);
tc.verifyEqual(S.ripp.info.clustSel, 1);
tc.verifyEqual(numel(S.ripp.clustId), numel(S.ripp.accepted));
tc.verifyEqual(S.ripp.info.qa.ranges.emg, [-Inf 1]);
end


function test_saveStripsSpks(tc)
% .spks is accepted-aligned, so a moved mask makes it wrong. ripp_analyze
% rebuilds it; leaving it would let a reader take stale rows for current ones.
f = fullfile(tc.TestData.dir, 'ripptest.ripp.mat');
S = load(f, 'ripp');
S.ripp.spks = struct('frac', rand(tc.TestData.nEv, 1));
ripp = S.ripp;
save(f, 'ripp', '-v7.3');

[~, hFig] = openGui(tc);
press(hFig, 'Save');

S2 = load(f, 'ripp');
tc.verifyFalse(isfield(S2.ripp, 'spks'), 'a stale .spks survived curation');
end


function test_saveInvalidatesProducts(tc)
% A moved mask must clear the accepted-aligned analyze products, so they cannot
% be read stale before ripp_analyze reruns.
stale = fullfile(tc.TestData.dir, 'ripptest.rippSpks.mat');
x = 1;        
save(stale, 'x');

[~, hFig] = openGui(tc);
st = hFig.UserData;
arrayfun(@(h) set(h, 'Value', false), st.chk);
st.chk(1).Value = true;
press(hFig, 'Save');

tc.verifyFalse(isfile(stale), 'a stale analyze product survived a mask change');
end


function test_saveKeepsProductsWhenAsked(tc)
% flgInvalidate = false is the batch path: ripp_analyze overwrites the products
% next, so deleting them would only cost a rewrite.
stale = fullfile(tc.TestData.dir, 'ripptest.rippSpks.mat');
x = 1;        
save(stale, 'x');

[~, hFig] = ripp_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', true, 'Visible', 'off', 'flgInvalidate', false, ...
    'verbose', false);
tc.addTeardown(@() close(hFig, 'force'));
st = hFig.UserData;
arrayfun(@(h) set(h, 'Value', false), st.chk);
st.chk(1).Value = true;
press(hFig, 'Save');

tc.verifyTrue(isfile(stale));
end


function test_resumesSavedCuration(tc)
% Reopening must resume: the same labels and the same accepted clusters, not a
% fresh clustering.
[~, h1] = openGui(tc);
st1 = h1.UserData;
st1.chk(2).Value = false;
st1.chk(2).ValueChangedFcn([], []);
cidWas = st1.cid;
nAcc   = nnz(h1.UserData.acc);
press(h1, 'Save');
close(h1, 'force');

[~, h2] = openGui(tc);
st2 = h2.UserData;
tc.verifyEqual(st2.cid, cidWas, 'labels were not restored');
tc.verifyEqual(nnz(st2.acc), nAcc, 'accepted set was not restored');
tc.verifyFalse(st2.chk(2).Value, 'cluster 2 came back ticked');
tc.verifySubstring(st2.lblPool.Text, 'restored');
end


function test_resumesKnobValues(tc)
% The saved spec must put the knobs back where they were, or a reopen silently
% widens the pool the next Re-cluster fits over.
[~, h1] = openGui(tc);
hKnob = knobOf(h1, 'spkGain');
hKnob.Value = 2.5;
press(h1, 'Re-cluster');
press(h1, 'Save');
close(h1, 'force');

[~, h2] = openGui(tc);
tc.verifyEqual(knobOf(h2, 'spkGain').Value, 2.5, ...
    'the knob came back at its shipped value');
end


function test_statesDefaultTickedAfterHeadless(tc)
% Every state box must open TICKED after a headless gate. A saved selection of
% nothing is read as "none saved", not as "every state was rejected".
ripp_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', false, 'verbose', false);

[~, hFig] = openGui(tc);
tc.verifyTrue(all(arrayfun(@(h) h.Value, hFig.UserData.chkState)), ...
    'states opened unticked after a headless gate');
end


function test_stateScopeIsSaved(tc)
% The state scope is recorded in the saved qa spec, so a headless replay of it
% reproduces the same mask.
[~, hFig] = openGui(tc);
st = hFig.UserData;
st.chkState(1).Value = false;           % the fixture has one state
st.chkState(1).ValueChangedFcn([], []);
press(hFig, 'Save');

S = load(fullfile(tc.TestData.dir, 'ripptest.ripp.mat'), 'ripp');
tc.verifyEqual(nnz(S.ripp.accepted), 0, ...
    'unticking the only state should accept nothing');
tc.verifyEqual(nnz(evt_gate(S.ripp, S.ripp.info.qa)), 0, ...
    'the saved spec does not reproduce the state scope');
end


function test_showDoesNotChangeMask(tc)
% The show dropdown moves rows only. Cycling it cannot alter what Save writes.
[~, hFig] = openGui(tc);
st = hFig.UserData;
arrayfun(@(h) set(h, 'Value', false), st.chk);
st.chk(1).Value = true;
st.chk(1).ValueChangedFcn([], []);
nExpect = nnz(st.cid == 1);

for v = {'accepted', 'removed', 'both'}
    st.ddShow.Value = v{1};
    st.ddShow.ValueChangedFcn([], []);
end

press(hFig, 'Save');
S = load(fullfile(tc.TestData.dir, 'ripptest.ripp.mat'), 'ripp');
tc.verifyEqual(nnz(S.ripp.accepted), nExpect);
end


function test_rejectThenRecluster(tc)
% Rejecting then re-clustering must fit ONLY what is still accepted, must keep
% the accepted set, and must not let the rejected shape back in.
[~, hFig] = openGui(tc);
st = hFig.UserData;
st.chk(1).Value = false;
st.chk(1).ValueChangedFcn([], []);
nAcc = nnz(hFig.UserData.acc);

press(hFig, 'Re-cluster');
st2 = hFig.UserData;
tc.verifyEqual(nnz(st2.pool), nAcc, 'clustered set is not the accepted set');
tc.verifyEqual(nnz(st2.acc), nAcc, 'the accepted set moved');
tc.verifyTrue(all(arrayfun(@(h) h.Value, st2.chk)));

press(hFig, 'Reset to filter');
tc.verifyGreaterThan(nnz(hFig.UserData.pool), nAcc, 'Reset did not widen');
end


function test_refusesStaleLabels(tc)
% Labels from a different event list must be refused, not drawn against the
% wrong events - a re-detection has to fall through to a fresh clustering.
f = fullfile(tc.TestData.dir, 'ripptest.ripp.mat');
S = load(f, 'ripp');
S.ripp.clustId = (1 : 5)';              % wrong length
S.ripp.info.clustSel = 1;
ripp = S.ripp;
save(f, 'ripp', '-v7.3');

[~, hFig] = openGui(tc);
tc.verifyEqual(numel(hFig.UserData.cid), tc.TestData.nEv);
tc.verifyGreaterThan(hFig.UserData.nClust, 1);
end


function test_rebuildsMissingMaps(tc)
% Waveforms are mandatory for the GUI, and a session whose maps are missing or
% belong to other events must say so rather than open a view of nothing.
S = load(fullfile(tc.TestData.dir, 'ripptest.rippMaps.mat'), 'rippMaps');
rippMaps = S.rippMaps;
rippMaps.lfp = rippMaps.lfp(1 : 10, :);         % no longer the event list
rippMaps.filt = rippMaps.filt(1 : 10, :);
save(fullfile(tc.TestData.dir, 'ripptest.rippMaps.mat'), 'rippMaps', '-v7.3');

% no session / lfp file in the fixture, so the rebuild path cannot complete -
% the point is that it is ATTEMPTED rather than the stale maps being drawn
tc.verifyError(@() ripp_curate(tc.TestData.dir, ...
    'basename', tc.TestData.name, 'flgGui', true, 'Visible', 'off', ...
    'verbose', false), ?MException);
end


%% ========================================================================
%  HELPERS
%  ========================================================================
function [ripp, hFig] = openGui(tc)
% The GUI on the fixture, invisible, torn down with the test.
[ripp, hFig] = ripp_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', true, 'Visible', 'off', 'verbose', false);
% tolerant, because a resume test closes the figure itself before reopening
tc.addTeardown(@() closeIfOpen(hFig));
end


function closeIfOpen(hFig)
if isgraphics(hFig), close(hFig, 'force'); end
end


function press(hFig, txt)
% Click a button by its label, the way a user would.
btn = findall(hFig, 'Type', 'uibutton', 'Text', txt);
btn.ButtonPushedFcn([], []);
end


function h = knobOf(hFig, fld)
% The metric edit box for FLD. evt_curate builds one knob per finite bound in
% met.qa.ranges rather than naming metrics in code, so a test looks its box up
% by the metric it belongs to.
st = hFig.UserData;
h = st.knob(strcmp({st.knob.fld}, fld)).h;
end
