function fixSpkAndRes(varargin)

% processes neurosuite files after /during manual curation. answers two
% issues: (1) re-aligns individual waveforms according to their peak /
% trough depending on the orientation of the mean waveform. (2) enforces a
% dead time on spike times (res) that typically occurs when kilosort tries
% to find spikes that occur simultaneously. re-alignment occurs first
% because it changes spike times. alas recalculates pca. if backup is
% flagged then copies the original files to a separate directory.
% -------------------------------------------------------------------------
% one problem that persists is that if a waveform includes two spikes, the
% larger spike will always be centered even if the cluster originally
% referred to the smaller spike.
% -------------------------------------------------------------------------
% update 04 jan 21; snipping spikes from whitened data reduces variability
% considerabily. Thus, the total isolation distance is greater (3233 vs
% 2078) and there are more su (87 vs 69) w/o fixSpkAndRes. Alternatively,
% fixSpk can be done on temp_wh data. update 25 mar 21; instead of
% re-snipping from temp_wh, circhshift can be used to recenter the peak /
% trough because the beginning and end of a spike are assumed to be zero
%
% INPUT:
%   basepath    string. path to recording folder {pwd}.
%   spkgrp      array where each cell is the electrodes for a spike group.
%               if empty will be loaded from session info file (cell
%               explorer format)
%   nchans      numeric. number of channels in dat file {35}
%   grp         numeric. groups (tetrodes) to work on. refers to filename
%               and not order. e.g. grp = 3 will load clu.3
%   fs          numeric. sampling frequency [hz]{20000}
%   nsamps      numeric. length of waveform in spk file {32}. if empty will
%               be calculated according to fs as 1.6 ms
%   psamp       numeric. peak / trough sample {16}. if empty will be
%               set to half nsamps.
%   stdFactor   numeric. re-alignment will occur only if std @ peak is
%               greater than std @ begining * stdFactor {1.1}. if 0
%               then all clusters will undergo re-alignment.
%   graphics    logical. plot the before / after of each cluster {false}.
%   bkup        logical. create backup of files before save in basepath/ns/bkup {true}.
%   saveFiles   logical. save ns files after processing {true}.
%   resnip      logical. resnip from dat / temp_wh (true) or use circshift
%               {false}
%
% DEPENDENCIES
%   loadNS
%   saveNS
%   plotCCG
%   plotWaveform
%
% TO DO lIST
%   # snip from temp_wh instead of circshift (done)
%   # find missing spk
%
% 22 jun 20 LH  UPDATES
% 24 mar 21 LH  embed load / save NS, the script now requires session
% 31 mar 21 LH  circshift 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'nchans', 35, @isnumeric);
addOptional(p, 'grp', [], @isnumeric);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'dt', 8, @isnumeric);
addOptional(p, 'nsamps', [], @isnumeric);
addOptional(p, 'psamp', [], @isnumeric);
addOptional(p, 'stdFactor', 0, @isnumeric);
addOptional(p, 'graphics', false, @islogical);
addOptional(p, 'bkup', true, @islogical);
addOptional(p, 'saveFiles', true, @islogical);
addOptional(p, 'resnip', false, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
spkgrp      = p.Results.spkgrp;
nchans      = p.Results.nchans;
grp         = p.Results.grp;
fs          = p.Results.fs;
dt          = p.Results.dt;
nsamps      = p.Results.nsamps;
psamp       = p.Results.psamp;
stdFactor   = p.Results.stdFactor;
graphics    = p.Results.graphics;
bkup        = p.Results.bkup;
saveFiles   = p.Results.saveFiles;
resnip      = p.Results.resnip;

% try to load params from session info (cell explorer format)
cd(basepath)
[~, basename] = fileparts(basepath);
sessionfile = fullfile(basepath, [basename '.session.mat']);
if exist(sessionfile, 'file')
    load(sessionfile)
end
if isempty(spkgrp)
    spkgrp = session.extracellular.spikeGroups.channels;
end
ngrp = length(spkgrp);
if isempty(nchans)
    nchans = length([spkgrp{:}]);
end
if isempty(fs)
    fs = session.extracellular.sr;
end
if isempty(nsamps)
    nsamps = ceil(1.6 * 10^-3 * fs);        % 1.6 ms in samples
end
if isempty(psamp)
    psamp = nsamps / 2;
end
if isempty(grp)
    grp = 1 : ngrp;
end

sniplength = ceil(1.6 * 10^-3 * fs);
win = [-(floor(sniplength / 2) - 1) floor(sniplength / 2)];
precision = 'int16'; % for dat file. size of one data point in bytes
nbytes = class2bytes(precision);

% for CCG
dur = 0.01;
binsize = 0.0001;

if resnip
    
    % build regressor for detrending
    s = 0 : sniplength - 1;
    scaleS = s(end);
    a = s./scaleS;
    b = max(a, 0);
    W = b(:);
    W = [reshape(W, sniplength, []), ones(sniplength,1)];
    [Q, R] = qr(W,0);
    
    % memory map to dat file
    fraw = dir([basename '.dat']);
    nsampsRaw = fraw.bytes / nbytes / nchans;
    m = memmapfile(fraw.name, 'Format', {precision, [nchans, nsampsRaw] 'mapped'});
    raw = m.Data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for igrp = grp
    grpchans = spkgrp{igrp};
    nchansGrp = length(spkgrp{igrp});
    
    % load data
    [clu, ~] = loadNS('datatype', 'clu', 'session', session, 'grpid', igrp);
    uclu = unique(clu);
    nclu = length(uclu);
    nspks = length(clu);
    
    res = loadNS('datatype', 'res', 'session', session, 'grpid', igrp);
    spk = loadNS('datatype', 'spk', 'session', session, 'grpid', igrp,...
        'nspks', nspks);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % go over clusters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    txt = '';
    for iclu = 1 : nclu
        
        % skip noise / artifact clusters
        if uclu(iclu) == 0 || uclu(iclu) == 1
            continue
        end
        
        % print progress
        fprintf(repmat('\b', 1, length(txt)))
        txt = sprintf('Working on cluster %d / %d. \t', iclu, nclu);
        fprintf('%s', txt)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % go over spikes and recenter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cluIdx = find(clu == uclu(iclu));
        nspksGrp = length(cluIdx);
        spkGrp = spk(:, :, cluIdx);
        spkMean = mean(spkGrp, 3);
        [~, ampMaxCh] = max(range(spkMean, 2));
        
        % re-align waveforms only if std at peak greater than at edges
        spkStd = std(spkGrp(ampMaxCh, :, :), [], 3);
        if spkStd(psamp) < mean(spkStd(1 : 5)) * stdFactor
            continue
        end
        
        % find wave polarity to decide if align to min or max
        if abs(min(spkMean(ampMaxCh, :))) > max(spkMean(ampMaxCh, :))
            spkPolarity = 'neg';
        else
            spkPolarity = 'pos';
        end
        
        % save original for plotting
        if graphics
            spkOrig = spkGrp;
            resOrig = res(cluIdx);
        end
        
        % go over each spike and align to min / max
        ishift = zeros(nspksGrp, 1);
        for ispk = 1 : nspksGrp
            if resnip
                spkidx = res(cluIdx(ispk));
                
                % fix special case where spike is at end / begining of recording
                if spkidx + win(1) < 1 || spkidx + win(2) > nsampsRaw
                    warning('\nskipping stamp %d because waveform incomplete', iclu)
                    clu(cluIdx(ispk)) = 0;
                    nspks(igrp) = nspks(igrp) - 1;
                    continue
                end
                
                % get waveform and remove best fit
                wv = double(raw.mapped(grpchans, spkidx + win(1) :...
                    spkidx + win(2)));
                wv = [wv' - W * (R \ Q' * wv')]';
            else
                wv = spkGrp(:, :, ispk);
            end
            
            % realign according to min / max
            switch spkPolarity
                case 'neg'
                    [~, ia] = min(wv, [], 2);
                case 'pos'
                    [~, ia] = max(wv, [], 2);
            end
            wvpeak = ia(ampMaxCh) ;
            ishift(ispk) = wvpeak - psamp;
            if ishift(ispk) ~= 0
                if resnip
                    spkidx = spkidx + ishift(ispk);
                    wv = double(raw.mapped(grpchans, spkidx + win(1) :...
                        spkidx + win(2)));
                    wv = [wv' - W * (R \ Q' * wv')]';
                    
                else
                    wv = circshift(wv, -ishift(ispk), 2);
                end
            end
            spkGrp(:, :, ispk) = wv;
        end
        spk(:, :, cluIdx) = spkGrp;
        clear spkGrp
        res(cluIdx) = res(cluIdx) + ishift;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % enforce dead time on spike detection
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        resGrp = res(cluIdx);
        [resGrp, rmIdx, ~] = unique(resGrp);
        rmIdx = [setdiff(1 : nspksGrp, rmIdx)]';
        while min(diff(resGrp)) < dt
            ir = find(diff(resGrp) < dt);
            resGrp(ir) = [];
            rmIdx = [rmIdx; ir];
        end
%         res(cluIdx(rmIdx)) = [];
        clu(cluIdx(rmIdx)) = 0;
%         spk(:, :, cluIdx(rmIdx)) = [];
%         nspks = length(clu);
        
        % graphics
        if graphics
            figure
            txt2 = sprintf('Group %d Cluster %d', igrp, uclu(iclu));
            suptitle(txt2)
            subplot(2, 2, 1)
            plotWaveform('avgwv', mean(spkOrig, 3),...
                'stdwv', std(spkOrig, [], 3), 'c', 'k',...
                'orient', 'vert', 'sbar', false)
            xticks([])
            yticks([])
            title('Before')
            subplot(2, 2, 2)
            plotWaveform('avgwv', mean(spkGrp(:, :, :), 3),...
                'stdwv', std(spkGrp(:, :, :), [], 3), 'c', 'k',...
                'orient', 'vert', 'sbar', false)
            xticks([])
            yticks([])
            title('After')
            subplot(2, 2, 3)
            [ccg, t] = CCG({resOrig / fs}, [], 'duration', dur, 'binSize', binsize);
            plotCCG('ccg', ccg(:, 1, 1), 't', t, 'basepath', basepath,...
                'saveFig', false, 'c', {'k'});
            subplot(2, 2, 4)
            cluIdx = clu == uclu(iclu);
            [ccg, t] = CCG({res(cluIdx) / fs}, [], 'duration', dur, 'binSize', binsize);
            plotCCG('ccg', ccg(:, 1, 1), 't', t, 'basepath', basepath,...
                'saveFig', false, 'c', {'k'});
            drawnow;
        end
    end
    
    % recalculate PCA
    fprintf('calculating PCA... \t')
    nfet = nchansGrp * 3 + nchansGrp + 1;
    fetMat = zeros(nspks, nfet);
    enrgIdx = nchansGrp * 3;
    if ~isempty(spk)
        for ich = 1 : nchansGrp
            [~, pcFeat] = pca(permute(spk(ich, :, :), [3, 2, 1]));
            chEnrgy = sum(abs(permute(spk(ich, :, :), [3, 2, 1])), 2);
            fetMat(:, ich * 3 - 2 : ich * 3) = (pcFeat(:, 1 : 3));
            fetMat(:, enrgIdx + ich) = (chEnrgy);
            clear pcFeat chEnrgy
        end
    end
    fetMat(:, end) = res;
    fet = int32(fetMat');
    fprintf('done. \n')

    % save files
    if saveFiles
        saveNS(clu, 'datatype', 'clu', 'session', session, 'grpid', igrp, 'bkup', bkup);
        saveNS(res, 'datatype', 'res', 'session', session, 'grpid', igrp, 'bkup', bkup);
        saveNS(spk, 'datatype', 'spk', 'session', session, 'grpid', igrp, 'bkup', bkup);
        saveNS(fet, 'datatype', 'fet', 'session', session, 'grpid', igrp, 'bkup', bkup);
    end
end

fprintf('that took %.2f minutes\n', toc / 60)

end

% EOF