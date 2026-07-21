function tests = test_edCurate
% TEST_EDCURATE Unit tests for ed_clust and the ed_curate cluster GUI.
%
%   tests = TEST_EDCURATE
%
%   SUMMARY:
%       Builds a synthetic session in a temp folder - two waveform families
%       plus noise - so the clustering and the GUI can be exercised without
%       touching a real recording or a curated mask.
%
%       Run with: runtests('test_edCurate')
%
%   HISTORY:
%       260721 created with the cluster curation rebuild.

tests = functiontests(localfunctions);
end


%% ========================================================================
%  FIXTURE
%  ========================================================================
function setupOnce(tc)
rng(7);
nA = 60; nB = 40; nN = 100;         % sharp, slow, noise
fs = 1250;
tst = -0.1 : 1 / fs : 0.1;
nS = numel(tst);

% family A: a sharp biphasic complex. family B: a slow negative deflection.
wvA = -400 * exp(-((tst + 0.004) / 0.003) .^ 2) ...
    + 1200 * exp(-(tst / 0.005) .^ 2);
wvB = -600 * exp(-((tst - 0.02) / 0.045) .^ 2);

% the third family is background, not white noise: real LFP is band-limited,
% and white noise would carry more variance than either real shape and take
% the principal components with it
bg = filtfilt(ones(1, 25) / 25, 1, 260 * randn(nN, nS)')' * 5;

wv = [repmat(wvA, nA, 1) .* (0.8 + 0.4 * rand(nA, 1)); ...
    repmat(wvB, nB, 1) .* (0.8 + 0.4 * rand(nB, 1)); bg];
wv = wv + filtfilt(ones(1, 9) / 9, 1, 40 * randn(size(wv))')' * 3;
nEv = size(wv, 1);

ed = struct();
ed.peakTime = sort(rand(nEv, 1) * 3600);
ed.times    = [ed.peakTime - 0.01, ed.peakTime + 0.01];
ed.fastZ    = [30 + randn(nA, 1) * 3; 18 + randn(nB, 1) * 3; ...
    16 + randn(nN, 1) * 2];
ed.isoZ     = 40 + randn(nEv, 1) * 5;
ed.posZ     = [12 + randn(nA, 1); -2 + randn(nB, 1); randn(nN, 1)];
ed.amp      = max(abs(wv), [], 2);
ed.dur      = [8 + randn(nA, 1); 22 + randn(nB, 1); 15 + randn(nN, 1)];
ed.state    = categorical(repmat({'NREM'}, nEv, 1));
ed.accepted = true(nEv, 1);
ed.info = struct('fs', fs, 'sigDur', 3600, 'win', [0 Inf], 'edCh', 1, ...
    'passband', [60 150], 'thr', 8, 'limDur', [4 200 40], 'met', 'test');

edMaps = struct('tstamps', tst, 'lfp', single(wv), 'filt', single(wv));

tc.TestData.dir = fullfile(tempdir, 'edtest');
if isfolder(tc.TestData.dir), rmdir(tc.TestData.dir, 's'); end
mkdir(tc.TestData.dir);
tc.TestData.name = 'edtest';
save(fullfile(tc.TestData.dir, 'edtest.ed.mat'), 'ed', '-v7.3');
save(fullfile(tc.TestData.dir, 'edtest.edMaps.mat'), 'edMaps', '-v7.3');

tc.TestData.wv = wv;
tc.TestData.tst = tst;
tc.TestData.truth = [ones(nA, 1); 2 * ones(nB, 1); 3 * ones(nN, 1)];
end


function teardownOnce(tc)
close all force
try
    rmdir(tc.TestData.dir, 's');    % a lingering matfile lock is not a failure
catch
end
end


%% ========================================================================
%  ED_CLUST
%  ========================================================================
function test_clustSeparatesFamilies(tc)
% The two planted waveform families must not land in one cluster.
cid = ed_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 3);

tc.verifyEqual(numel(cid), size(tc.TestData.wv, 1));
tc.verifyTrue(all(ismember(cid, 1 : 3)));

truth = tc.TestData.truth;
domA = mode(cid(truth == 1));
domB = mode(cid(truth == 2));
tc.verifyNotEqual(domA, domB, ...
    'sharp and slow families collapsed into one cluster');

% and each family should be mostly in its dominant cluster
tc.verifyGreaterThan(mean(cid(truth == 1) == domA), 0.7);
tc.verifyGreaterThan(mean(cid(truth == 2) == domB), 0.7);
end


function test_clustDefaultK(tc)
% The default count applies without being named, and is reported back.
[cid, cInfo] = ed_clust(tc.TestData.wv, tc.TestData.tst);
tc.verifyEqual(cInfo.nClust, 12);
tc.verifyEqual(max(cid), cInfo.nClust);
end


function test_clustCapsKToPoolSize(tc)
% A count larger than the pool can support must be capped, not error.
[cid, cInfo] = ed_clust(tc.TestData.wv(1 : 30, :), tc.TestData.tst, ...
    'nClust', 12);
tc.verifyLessThanOrEqual(cInfo.nClust, 10);
tc.verifyGreaterThanOrEqual(cInfo.nClust, 2);
tc.verifyEqual(max(cid), cInfo.nClust);
end


function test_clustTooFewEvents(tc)
% Below the minimum, every label is NaN rather than a meaningless fit. The
% GUI reads that as "too few to cluster".
[cid, cInfo] = ed_clust(tc.TestData.wv(1 : 15, :), tc.TestData.tst);
tc.verifyEqual(cInfo.nClust, 0);
tc.verifyTrue(all(isnan(cid)));
end


function test_clustLabelsBySize(tc)
% Cluster 1 must be the largest, so an index means the same thing on a re-run.
cid = ed_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 4);
cnt = accumarray(cid(~isnan(cid)), 1);
tc.verifyEqual(cnt, sort(cnt, 'descend'), ...
    'clusters are not ordered by size');
end


function test_clustNanWaveform(tc)
% An event at a recording edge has a NaN waveform: it must get a NaN label,
% not poison the fit.
wv = tc.TestData.wv;
wv(1, :) = NaN;
cid = ed_clust(wv, tc.TestData.tst, 'nClust', 3);
tc.verifyTrue(isnan(cid(1)));
tc.verifyFalse(any(isnan(cid(2 : end))));
end


function test_clustScalarFeatures(tc)
% Scalar measures must be accepted and must not break the labelling.
S = [tc.TestData.wv(:, 1), max(tc.TestData.wv, [], 2)];
cid = ed_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 3, 'scalar', S);
tc.verifyTrue(all(ismember(cid, 1 : 3)));
end


%% ========================================================================
%  ED_CURATE
%  ========================================================================
function test_curateHeadless(tc)
% Headless applies the noise filter only, and writes a mask of that size.
met = ed_methods('default');
ed = ed_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'met', met, 'flgGui', false, 'verbose', false);

S = load(fullfile(tc.TestData.dir, 'edtest.ed.mat'), 'ed');
expect = evt_gate(ed, met.qa);
tc.verifyEqual(S.ed.accepted, expect);
end


function test_curateGuiOpens(tc)
% The GUI must build, cluster, and expose one checkbox per cluster.
[~, hFig] = ed_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', true, 'Visible', 'off');
tc.addTeardown(@() close(hFig, 'force'));

st = hFig.UserData;
tc.verifyGreaterThan(st.nClust, 1);
tc.verifyEqual(numel(st.chk), st.nClust);
tc.verifyEqual(nnz(~isnan(st.cid)), nnz(st.pool));
tc.verifyEqual(numel(st.cid), numel(st.pool));
end


function test_curateSaveAcceptsClusters(tc)
% Ticking a cluster must accept exactly its events and record the choice.
[~, hFig] = ed_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', true, 'Visible', 'off');
tc.addTeardown(@() close(hFig, 'force'));

st = hFig.UserData;
st.chk(1).Value = true;
nExpect = nnz(st.cid == 1);

btn = findall(hFig, 'Type', 'uibutton', 'Text', 'Save');
btn.ButtonPushedFcn([], []);

S = load(fullfile(tc.TestData.dir, 'edtest.ed.mat'), 'ed');
tc.verifyEqual(nnz(S.ed.accepted), nExpect);
tc.verifyEqual(S.ed.info.clustSel, 1);
tc.verifyEqual(numel(S.ed.clustId), numel(S.ed.accepted));
end


function test_curateViewPivots(tc)
% The hosted guiTbl_xy must offer BOTH cluster and state as tile variables,
% and lfp as the Y variable - that is the pivot the GUI exists to give.
[~, hFig] = ed_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', true, 'Visible', 'off');
tc.addTeardown(@() close(hFig, 'force'));

ud = hFig.UserData.hPanel.UserData;
tc.verifyTrue(isfield(ud, 'setDataFcn'));
tc.verifyTrue(ismember('cluster', ud.ddPlotBy.Items));
tc.verifyTrue(ismember('state', ud.ddPlotBy.Items));
tc.verifyTrue(ismember('status', ud.ddGrpBy.Items));
tc.verifyEqual(ud.yVar, 'lfp');
tc.verifyEqual(height(ud.dataTbl), numel(hFig.UserData.cid));
end


function test_curateKnobsChangePool(tc)
% Raising a threshold must shrink the pool on Re-cluster, and the view must
% follow. The knob alone only updates the label.
[~, hFig] = ed_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', true, 'Visible', 'off');
tc.addTeardown(@() close(hFig, 'force'));

n0 = nnz(hFig.UserData.pool);
hFig.UserData.edFast.Value = 28;
btn = findall(hFig, 'Type', 'uibutton', 'Text', 'Re-cluster');
btn.ButtonPushedFcn([], []);

tc.verifyLessThan(nnz(hFig.UserData.pool), n0);
tc.verifyEqual(nnz(~isnan(hFig.UserData.cid)), nnz(hFig.UserData.pool));
end


function test_curateSmallPool(tc)
% A pool too small to cluster must not error: no checkboxes, no tiles, and a
% Save that writes an all-false mask.
S = load(fullfile(tc.TestData.dir, 'edtest.ed.mat'), 'ed');
ed = S.ed;
ed.fastZ(:) = 0;                    % nothing passes the noise filter
ed.fastZ(1 : 12) = 99;              % except twelve, which is below MINEV
ed.isoZ(:) = 99;
save(fullfile(tc.TestData.dir, 'edsmall.ed.mat'), 'ed', '-v7.3');
copyfile(fullfile(tc.TestData.dir, 'edtest.edMaps.mat'), ...
    fullfile(tc.TestData.dir, 'edsmall.edMaps.mat'));

[~, hFig] = ed_curate(tc.TestData.dir, 'basename', 'edsmall', ...
    'flgGui', true, 'Visible', 'off');
tc.addTeardown(@() close(hFig, 'force'));

tc.verifyEqual(hFig.UserData.nClust, 0);
tc.verifyEmpty(hFig.UserData.chk);

btn = findall(hFig, 'Type', 'uibutton', 'Text', 'Save');
btn.ButtonPushedFcn([], []);
S2 = load(fullfile(tc.TestData.dir, 'edsmall.ed.mat'), 'ed');
tc.verifyEqual(nnz(S2.ed.accepted), 0);
end


function test_curateRecusterKeepsTicks(tc)
% Re-clustering at the SAME count must keep the accepted clusters. This is the
% bug that made the first design unusable.
[~, hFig] = ed_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', true, 'Visible', 'off');
tc.addTeardown(@() close(hFig, 'force'));

st = hFig.UserData;
st.chk(1).Value = true;
st.chk(2).Value = true;
kWas = st.nClust;

btn = findall(hFig, 'Type', 'uibutton', 'Text', 'Re-cluster');
btn.ButtonPushedFcn([], []);

st2 = hFig.UserData;
tc.verifyEqual(st2.nClust, kWas);
tc.verifyTrue(st2.chk(1).Value, 'cluster 1 tick lost on re-cluster');
tc.verifyTrue(st2.chk(2).Value, 'cluster 2 tick lost on re-cluster');
end


function test_curateRecusterClearsOnCountChange(tc)
% Changing the count must clear the ticks - index 3 of 12 is not index 3 of 8.
[~, hFig] = ed_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', true, 'Visible', 'off');
tc.addTeardown(@() close(hFig, 'force'));

hFig.UserData.chk(1).Value = true;
hFig.UserData.edK.Value = 5;
btn = findall(hFig, 'Type', 'uibutton', 'Text', 'Re-cluster');
btn.ButtonPushedFcn([], []);

st = hFig.UserData;
tc.verifyEqual(st.nClust, 5);
tc.verifyFalse(any(arrayfun(@(h) h.Value, st.chk)));
tc.verifySubstring(st.lblPool.Text, 'ticks cleared');
end


function test_curateStateGatesAcceptance(tc)
% Unticking a state must remove its events from the mask without touching the
% cluster choice.
[~, hFig] = ed_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', true, 'Visible', 'off');
tc.addTeardown(@() close(hFig, 'force'));

st = hFig.UserData;
st.chk(1).Value = true;
st.chk(1).ValueChangedFcn([], []);
nWith = nnz(ismember(st.cid, 1));
tc.verifyGreaterThan(nWith, 0);

st.chkState(1).Value = false;       % the fixture has one state only
st.chkState(1).ValueChangedFcn([], []);

btn = findall(hFig, 'Type', 'uibutton', 'Text', 'Save');
btn.ButtonPushedFcn([], []);
S = load(fullfile(tc.TestData.dir, 'edtest.ed.mat'), 'ed');
tc.verifyEqual(nnz(S.ed.accepted), 0, ...
    'unticking the only state should accept nothing');
tc.verifyTrue(hFig.UserData.chk(1).Value, 'cluster tick disturbed');
end


function test_curateShowDoesNotChangeMask(tc)
% The show dropdown must move rows only. Cycling it cannot alter what Save
% writes - that separation is the point of having it.
[~, hFig] = ed_curate(tc.TestData.dir, 'basename', tc.TestData.name, ...
    'flgGui', true, 'Visible', 'off');
tc.addTeardown(@() close(hFig, 'force'));

st = hFig.UserData;
st.chk(1).Value = true;
st.chk(1).ValueChangedFcn([], []);
nExpect = nnz(ismember(st.cid, 1));

for v = {'accepted', 'removed', 'both'}
    st.ddShow.Value = v{1};
    st.ddShow.ValueChangedFcn([], []);
end

btn = findall(hFig, 'Type', 'uibutton', 'Text', 'Save');
btn.ButtonPushedFcn([], []);
S = load(fullfile(tc.TestData.dir, 'edtest.ed.mat'), 'ed');
tc.verifyEqual(nnz(S.ed.accepted), nExpect);
end


function test_clustDeterministic(tc)
% The same pool and count must give the same partition, so pressing
% Re-cluster without changing anything is a no-op rather than a reshuffle.
a = ed_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 8);
b = ed_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 8);
tc.verifyEqual(a, b, 'clustering is not reproducible');

% and it must not leave the global RNG stream disturbed
rng(42); r1 = rand();
rng(42); ed_clust(tc.TestData.wv, tc.TestData.tst, 'nClust', 4); r2 = rand();
tc.verifyEqual(r1, r2, 'ed_clust leaked its RNG seed to the caller');
end
