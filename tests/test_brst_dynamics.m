
%% test_brst_dynamics.m
% Verifies the functionality of brst_dynamics.m (Row-based Refactor)

% 1. GENERATE SYNTHETIC DATA
% -------------------------------------------------------------------------
spktimes = cell(2, 1); % Column vector of cells

% UNIT 1: Continuous Activity (Rate Increase)
st1 = [];
for t = 1:1:9
    b = t + (0:4)*0.005;
    st1 = [st1, b];
end
for t = 10:0.5:20
    b = t + (0:4)*0.005;
    st1 = [st1, b];
end
spktimes{1} = st1(:);

% UNIT 2: Gapped Activity (Masking Test)
% Bursts 0-5s, then Silence 5-15s, then Bursts 15-20s.
st2 = [];
% Phase 1: Short bursts
for t = 1:0.5:5
    b = t + (0:2)*0.005;
    st2 = [st2, b];
end
% Phase 2: Long bursts (after gap)
for t = 15:0.5:20
    b = t + (0:9)*0.005;
    st2 = [st2, b];
end
spktimes{2} = st2(:);


% 2. RUN DETECTION
% -------------------------------------------------------------------------
brst = brst_maxInt(spktimes, 'flgPlot', false, 'flgSave', false, ...
    'minSpks', 3, 'minDur', 0, 'maxISI_start', 0.05, 'maxISI_end', 0.05, ...
    'flgForce', true);

% Check orientation of outputs
assert(all(size(brst.detect) == [2 1]), 'brst.detect should be 2x1 (col vector)');


% 3. RUN DYNAMICS
% -------------------------------------------------------------------------
% Using ibiPct=80 to ensure gap masking
dyn = brst_dynamics(brst, spktimes, 'binSize', 0.1, 'kernelSD', 1, ...
    'flgPlot', false, 'ibiPct', 80);

% 4. VALIDATE
% -------------------------------------------------------------------------

% A. MATRIX STRUCTURE (nUnits x nTime)
nTime = length(dyn.time);
assert(all(size(dyn.rate) == [2, nTime]), 'Dyn matrices must be nUnits x nTime');
assert(all(size(dyn.dur) == [2, nTime]), 'Dyn matrices must be nUnits x nTime');


% B. MASKING TEST (Unit 2, Row 2)
% Gap is from t=5 to t=15.
tGap = dyn.time > 6 & dyn.time < 14;

% Using row index 2
nanFrac = mean(isnan(dyn.dur(2, tGap)));

fprintf('Fraction of NaNs in Gap for Unit 2: %.2f\n', nanFrac);
assert(nanFrac > 0.9, 'Gap was not masked correctly (should be NaNs)');


% C. INTERPOLATION TEST (Unit 2, Row 2)
tEarly = dyn.time > 1 & dyn.time < 4;
tLate  = dyn.time > 16 & dyn.time < 19;

durEarly = mean(dyn.dur(2, tEarly), 'omitnan');
durLate  = mean(dyn.dur(2, tLate), 'omitnan');

fprintf('Unit 2 Dur: Early=%.4f, Late=%.4f\n', durEarly, durLate);
assert(durLate > durEarly * 2, 'Property values incorrect after interpolation');

disp('TEST PASSED: brst_dynamics Matrix Orientation Verified.');
