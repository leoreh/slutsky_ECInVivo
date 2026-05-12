%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   basepath    str. path to excel file
%   xlsName     str. name of excel file

basepath = 'C:\Users\Maxim\Desktop\New- Not analyzed';
xlsName = 'sample data input format';

% load from excel and organize vars
xlsData = readtable(fullfile(basepath, xlsName));
sig = xlsData.F;
spktimes = xlsData.SpikeTimes(~isnan(xlsData.SpikeTimes)) / 1000;
caEvents = xlsData.CalciumEventFrameNumber(~isnan(xlsData.CalciumEventFrameNumber));
dff = xlsData.CacliumEventAmplitude_df_f_(~isnan(xlsData.CacliumEventAmplitude_df_f_));

% additional params (from excel or manual). binTol refers to the minimal
% distance (in bins) between a spike and event in order for them to be
% considered a pair. The 1st element describes the minimum when a spike
% occurs before an event, and the 2nd element describes the minimum when an
% event occurs before a spike. If there is an event before and after a
% spike, the one after the spike gets precedence in the assignment
fs = 12.98
binTol = [80, 0];       

% manually align spktimes and calcium
manShift = manAlign('spktimes', spktimes, 'sig', sig, 'fs', fs);

% apply shift to spktimes and save original
origSpktimes = spktimes;
spktimes = spktimes - manShift;

% recording duration from calcium fs
tstamps = ([1 : length(sig)] / fs);
recDur = tstamps(end);

% count spikes in bins corresponding to calcium
binsSpks = [0 : 1 / fs : spktimes(end)];
nspks = histcounts(spktimes, binsSpks, 'Normalization', 'count');

% find the peak closest to each spike, but only if it is +/- from binTol.
idxSpks = [find(nspks > 0)]';           % indices to bins with spikes
idxCa = nan(length(idxSpks), 1);        % initialize
for ibin = 1 : length(idxSpks)
       
    % get minimum distance in the forward direction (ca peak occurs before spike)
    diffVal = idxSpks(ibin) - caEvents; 
    [dFor, iFor] = min(diffVal(diffVal > 0));

    % get minimum distance in the backward direction (spike occurs before peak)
    diffVal(diffVal > 0) = -Inf;
    [dRev, iRev] = min(-diffVal);

    % check and assign. the backgward direction gets precedence. 
    if dRev <= binTol(1)
        idxCa(ibin) = caEvents(iRev);
    elseif dFor <= binTol(2)
        idxCa(ibin) = caEvents(iFor);
    else
        idxCa(ibin) = nan;
    end
end

% remove spk bins that have no corresponding ca peak
keepIdx = ~isnan(idxCa);      

% keep only bins with ca peaks
spkCnt = nspks(idxSpks(keepIdx))';
idxSpks = idxSpks(keepIdx);
idxCa = idxCa(keepIdx);

% find the unique number of ca peaks and sum their corresponding number of
% spks
uPeaks = unique(idxCa);
sumSpks = zeros(size(uPeaks));
for ipeak = 1 : length(uPeaks)
    sumSpks(ipeak) = sum(spkCnt(idxCa == uPeaks(ipeak)));
end

% find dff that corresponds to the unique ca peaks 
dffPeaks = zeros(size(uPeaks));
for ipeak = 1 : length(uPeaks)
f=find(caEvents == uPeaks(ipeak));
if length(f)>1,f=f(1); end
dffIdx = f;
%     dffIdx = find(caEvents == uPeaks(ipeak));
    dffPeaks(ipeak) = dff(dffIdx);
end

% organize output mat
outData = [uPeaks, dffPeaks, sumSpks];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% blue is F. black dashed lines are individual spikes. numbers at
% the top above each peak represent the number of spikes corresponding to
% that peak. number in green at the peaks corresponds to the df/f of that
% peak. 

fh = figure;
plot(tstamps, sig);
yLimit = ylim;
hold on
plot([spktimes spktimes], yLimit, '--', 'Color', [0.5 0.5 0.5])
text(uPeaks / fs, ones(1, length(uPeaks)) * yLimit(2), num2str(sumSpks),...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
text(uPeaks / fs, sig(uPeaks), num2str(dffPeaks, '%.2f'), 'Color', [0.2 0.7 0.2],...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
xlabel('Time (s)')
ylabel('F')
