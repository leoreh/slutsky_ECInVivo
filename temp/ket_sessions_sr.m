% ket_sessions_sr

forceL = true;
normFlag = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ket 10 mg/kg i.p.
basepaths = [{'K:\Data\lh99\lh99_211224_084528'},...
{'G:\Data\lh98\lh98_211224_084528'},...    
{'F:\Data\Processed\lh96\lh96_211206_070400'},...
    {'G:\Data\lh81\lh81_210204_190000'}];
% basepaths = basepaths(1 : 3);
basepaths = basepaths(4);

% saline i.p.
basepaths = [{'F:\Data\Processed\lh96\lh96_211205_072000'},...
    {'G:\Data\lh81\lh81_210204_190000'},...
    {'K:\Data\lh99\lh99_211222_090505'}];

% ket 60 mg/kg i.p.
basepaths = [{'G:\Data\lh84\lh84_210504_225608'},...
    {'G:\Data\lh81\lh81_210206_190000'}];

% load vars from each session
varsFile = ["sr"; "st_metrics"; "swv_metrics";...
    "cell_metrics"; "sleep_states"; "datInfo"; "session"];
varsName = ["sr"; "st"; "swv"; "cm"; "ss";...
    "datInfo"; "session"];
if ~exist('v', 'var') || forceL
    v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName);
end
nsessions = length(basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nunits = [];
units = [];
injIdx = [];
tLen = [];
for isession = 1 : nsessions
    
    % session params
    basepath = basepaths{isession};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    
    fs = v(isession).session.extracellular.sr;
    fsLfp = v(isession).session.extracellular.srLfp;
    spkgrp = v(isession).session.extracellular.spikeGroups.channels;
    nchans = v(isession).session.extracellular.nChannels;
    
    if contains(basename, 'lh99')
        grp = [1, 3 : 4, 7];
        v(isession).sr.strd = v(isession).sr.strd(grp, :);
    else
        grp = [];
    end
  
    % plot fr vs. time   
%     plot_FRtime_session('basepath', pwd, 'grp', grp,...
%         'frBoundries', [0.01 Inf; 0.01 Inf], 'muFlag', false, 'saveFig', false,...
%         'dataType', 'strd')

%     plot_FRtime_session('basepath', pwd, 'grp', grp,...
%         'frBoundries', [0.01 Inf; 0.01 Inf], 'muFlag', false, 'saveFig', false,...
%         'dataType', 'norm')
    
    % timebins
    fileinfo = dir([basename, '.dat']);
    recLen = floor(fileinfo.bytes / 2 / nchans / fs);
    csec = floor(cumsum(v(isession).datInfo.nsamps / fs));
    [~, pntIdx] = min(abs(csec - 12 * 60 * 60));
    timepoints = csec(pntIdx);
    % timepoints = v(isession).datInfo.nsec;
    chunks = n2nchunks('n', recLen, 'nchunks', 8, 'timepoints', timepoints);
    
    [~, injIdx(isession)] = min(abs(v(isession).sr.tstamps - timepoints));
    tLen(isession) = length(v(isession).sr.tstamps);
    
    nunits(isession) = size(v(isession).sr.strd, 1);
end

maxIdx = max(injIdx);
frMat = nan(sum(nunits), maxIdx + max(tLen));

cnt = 1;
dataType = 'strd';
for isession = 1 : nsessions
    startIdx = maxIdx - injIdx(isession) + 1;
    frIdx = startIdx : tLen(isession) + startIdx - 1;
    frMat(cnt : cnt + nunits(isession) - 1, frIdx) = v(isession).sr.(dataType);
    cnt = cnt + nunits(isession);
end

% noramlize
if normFlag
    frNorm = frMat ./ mean(frMat(:, 1 : maxIdx), 2, 'omitnan');
    frNorm(frNorm == Inf) = nan;
    ytxt = 'Norm Firing Rate';
else
    frNorm = frMat;
    ytxt = 'Firing Rate [Hz]';
end

xidx = [1 : length(frNorm)] / 60;
tidx = maxIdx / 60;

xidx2 = [12 : 1 / 60 : 24, 0 : 1 / 60 : 12];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
plot(xidx, mean(frNorm, 1, 'omitnan'))

