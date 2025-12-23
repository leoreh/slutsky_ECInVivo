
% Test script for mea_tAlign.m

clear; clc;

% --- 1. Synthesize Data ---
% Create a mock 'v' struct array with 2 files.
% File 1: Short pre-pert, Long post-pert
% File 2: Long pre-pert, Short post-pert

% Common params
binSize = 1;
nUnits = 2;

% File 1
idxPert1 = 5;
nBins1 = 15;
fr1 = rand(nUnits, nBins1);
dyn1_rate = rand(nUnits, nBins1);
info1.idxPert = idxPert1;
info1.binSize = binSize;

v(1).fr.fr = fr1;
v(1).fr.info = info1;
v(1).dyn.rate = dyn1_rate;

% File 2
idxPert2 = 10;
nBins2 = 14;
fr2 = rand(nUnits, nBins2);
dyn2_rate = rand(nUnits, nBins2);
info2.idxPert = idxPert2;
info2.binSize = binSize;

v(2).fr.fr = fr2;
v(2).fr.info = info2;
v(2).dyn.rate = dyn2_rate;


% --- 2. Define varMap ---
% Case A: Mapping specific numeric arrays
varMap.fr = 'fr.fr';
varMap.rate = 'dyn.rate';

% --- 3. Run mea_tAlign ---
fprintf('Running mea_tAlign...\n');
[vOut, tGlobal] = mea_tAlign(v, 'varMap', varMap, 'refVar', 'fr');

% --- 4. Validation ---

% 4.1 Check Dimensions
% Max Pre should be max(5-1, 10-1) = 9 bins
% Max Post should be max(15-5, 14-10) = 10 bins
% Total bins = 9(pre) + 1(pert) + 10(post) = 20 bins

expectedBins = 20;
expectedIdxPert = 10; % 9 pre + 1

fprintf('Expected Bins: %d, Actual: %d\n', expectedBins, length(tGlobal));
assert(length(tGlobal) == expectedBins, 'Global time vector length incorrect');

% 4.2 Check Padding File 1
% File 1 had idxPert=5 (4 pre). New idxPert is 10 (9 pre).
% Should be padded by 5 NaNs at start.
% File 1 had 10 post. New post is 10. No padding at end.

fr1_aligned = vOut(1).fr.fr;
if all(isnan(fr1_aligned(:, 1:5)), 'all')
    fprintf('File 1 Pre-padding (Start): OK\n');
else
    error('File 1 Pre-padding incorrect');
end

% 4.3 Check Padding File 2
% File 2 had idxPert=10 (9 pre). New is 10. No padding at start.
% File 2 had 4 post. New is 10. Should be padded by 6 NaNs at end.

fr2_aligned = vOut(2).fr.fr;
if all(isnan(fr2_aligned(:, end-5:end)), 'all')
    fprintf('File 2 Post-padding (End): OK\n');
else
    error('File 2 Post-padding incorrect');
end


fprintf('Validation Passed!\n');
