%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
[mousepath, basename] = fileparts(basepath);

% load and assign
varArray = getSessionVars('dirnames', {basename}, 'mousepath', mousepath,...
    'sortDir', false);
assignVars(varArray, 1)

% params
fs = session.extracellular.sr;
fsLfp = session.extracellular.srLfp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% separate ripples 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

csum_blockDur = cumsum(datInfo.nsec);
dark_phase_start = csum_blockDur(3);
inj_time = csum_blockDur(1);

% rippels in NREM
idx_nrem = InIntervals(ripp.peakPos, ss.stateEpochs{4});
idx_pre = ripp.peakPos < inj_time & idx_nrem;
idx_post = ripp.peakPos > inj_time & ripp.peakPos < dark_phase_start & idx_nrem;

% percent ripples that occured in nrem
percent_nrem_pre = sum(idx_pre) / sum(ripp.peakPos < inj_time);
percent_nrem_post = sum(idx_post) /...
    sum(ripp.peakPos > inj_time & ripp.peakPos < dark_phase_start);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
histBins = 200;

fh = figure;
h = histogram(ripp.peakFreq(idx_pre), histBins, 'Normalization', 'probability');
h.FaceColor = 'k';
h.EdgeColor = 'none';
h.FaceAlpha = 0.3;     
hold on
h = histogram(ripp.peakFreq(idx_post), histBins, 'Normalization', 'probability');
h.FaceColor = 'b';
h.EdgeColor = 'none';
h.FaceAlpha = 0.3;     

xlabel('Peak Frequency [Hz]')
ylabel('Probability')


fh = figure;
peakFreq = ripp.peakFreq(idx_pre)
bh = boxplot(ripp.peakFreq
