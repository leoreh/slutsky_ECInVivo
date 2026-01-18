function tests = test_spk_sttc
tests = functiontests(localfunctions);
end

function testIdentity(testCase)
% Test that STTC of a train with itself is 1
spk = sorrt([1, 2, 3, 4, 5]');
dt = 0.1;
tLim = [0, 6];

[val, ~] = spk_sttc(spk, spk, dt, tLim);
verifyEqual(testCase, val, 1, 'RelTol', 1e-10);
end

function testEmpty(testCase)
% Test empty inputs return NaN
spk1 = [1, 2];
spk2 = [];
dt = 0.1;
tLim = [0, 10];

val = spk_sttc(spk1, spk2, dt, tLim);
verifyTrue(testCase, isnan(val));

val = spk_sttc([], [], dt, tLim);
verifyTrue(testCase, isnan(val));
end

function testKnownOverlap(testCase)
% Test a simple case we calculated manually
% spk1 = [2], spk2 = [2.2], dt = 0.5, tLim = [0, 10]
% Matches logic -> STTC = 1

spk1 = 2;
spk2 = 2.2;
dt = 0.5;
tLim = [0, 10];

[val, info] = spk_sttc(spk1, spk2, dt, tLim);

verifyEqual(testCase, info.PA, 1);
verifyEqual(testCase, info.PB, 1);

% TA = 2*dt*N / T = 1/10 = 0.1
verifyEqual(testCase, info.TA, 0.1);
verifyEqual(testCase, info.TB, 0.1);

verifyEqual(testCase, val, 1, 'AbsTol', 1e-10);
end

function testNoOverlapAndBoundary(testCase)
% spk1 = [2], spk2 = [10], dt = 1, tLim = [0, 10]
% No match.
% TA: at 2. 2*dt = 2. No boundary overlap. -> 2. Fraction 0.2.
% TB: at 10. 2*dt = 2. Overlap at end: (10+1) - 10 = 1. TA_raw = 2-1 = 1. Fraction 0.1.

spk1 = 2;
spk2 = 10;
dt = 1;
tLim = [0, 10];

[val, info] = spk_sttc(spk1, spk2, dt, tLim);

verifyEqual(testCase, info.PA, 0);
verifyEqual(testCase, info.PB, 0);
verifyEqual(testCase, info.TA, 0.2);
verifyEqual(testCase, info.TB, 0.1);

% Term1 = (0 - 0.1)/(1-0) = -0.1
% Term2 = (0 - 0.2)/(1-0) = -0.2
% STTC = 0.5 * -0.3 = -0.15
verifyEqual(testCase, val, -0.15, 'AbsTol', 1e-10);
end
