function cleanCluByFet(varargin)

% removes rpv spikes on fet / mDist

% alternative option tested was to compare distributions with the
% Kullback–Leibler divergence, e.g.: kl(ifet) = kldiv(bincents', spkpdf' + eps, rpvpdf' + eps, 'sym');
%
% INPUT:
%   basepath    string. path to recording folder {pwd}.
%   grps        numeric. groups (tetrodes) to work on
%   ref         numeric. minimium isi for spike to be defined as rvd [s]{0.002}
%   rpvThr      remove entire cluster if rpv ratio after clean is still
%               greater than rpvTHr
%   spkgrp      array where each cell is the electrodes for a spike group
%   fs          numeric. sample frequency. if empty will try to take from
%               session (cell explorer format)
%
% DEPENDENCIES
%   loadNS
%   saveNS
%   cluDist
%   fixSpkAndRes
%
% TO DO LIST
%   # remove all global vars except res, fet and clu
%   # track total number of spks removed from cluster
%
% 15 mar 21 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'grps', [], @isnumeric);
addOptional(p, 'ref', 0.002, @isnumeric);
addOptional(p, 'rpvThr', 2, @isnumeric);
addOptional(p, 'spkgrp', {}, @iscell);
addOptional(p, 'fs', [], @isnumeric);
addOptional(p, 'graphics', true, @islogical);
addOptional(p, 'manCur', false, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
grps        = p.Results.grps;
ref         = p.Results.ref;
rpvThr      = p.Results.rpvThr;
spkgrp      = p.Results.spkgrp;
fs          = p.Results.fs;
graphics    = p.Results.graphics;
manCur    = p.Results.manCur;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(basepath)
[~, basename] = fileparts(basepath);
if isempty(fs) || isempty(spkgrp)
    session = CE_sessionTemplate(pwd, 'viaGUI', false,...
        'forceDef', true, 'saveVar', true);
    fs = session.extracellular.sr;
    spkgrp = session.extracellular.spikeGroups.channels;
end

npca = 3;
ref = ref * fs;
if isempty(grps)
    grps = 1 : length(spkgrp);
end

% constants
rmvLim = [0.02 0.25];   % ratio of spikes allowed to be removed from one feature and in total
nspksMin = 10000;       % cluster with less spikes will not be cleaned
cdfThr = 0.1;           % threshold for cdf difference from which to clean cluster
rpvRatioThr = 0.35;     % threshold for rpv ratio above which to clean cluster
rpvCrt = 0.5;           % rpv ratio cluster much reach in while loop

% shared vars between nested and parent
global clu
db1 = []; db2 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manual inspection of features and rpvs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if manCur
    dbstop in cleanCluByFet at 94 if manCur
    
    % select spiking grp
    grpid = 4;
    
    % load data
    res = loadNS('datatype', 'res', 'session', session, 'grpid', grpid);
    fet = loadNS('datatype', 'fet', 'session', session, 'grpid', grpid);

    % print to screen rpv ratio of all clusters in group
    clu = loadNS('datatype', 'clu', 'session', session, 'grpid', grpid);
    rpvRatioAll = displayRPV(res, ref)

    % remove rpv spks from all clusters until criterion is reached
    sclu = [];      % selected clusters. if empty will clean all
    clean2criterion(fet, res, ref, rpvCrt, rmvLim, sclu, grpid, false)
            
    % save new clu
    saveNS(clu, 'datatype', 'clu', 'session', session, 'grpid', grpid);      
    
    % calc quality of cluster separation
    distofclus = displayCluDist(fet, [])
        
    % realign spikes to peak / trough
    fixSpkAndRes('grp', grpid, 'dt', 0, 'stdFactor', 0, 'resnip', false);
    
    % rmv spks from selected cluster
    sclu = 5;      % cluster id
    for sclu = 2 : max(unique(clu))
        [cluidx, rmvSpks] = cleanClu(fet, res, sclu, ref, 0.0001, [0.12 0.06]);
    end
        
    % plot fets of specific cluster
    fetclu = fet(clu == cluid, :);
    [rpv, rpvRatioOrig] = getRpv(res(clu == sclu), ref);
    [~, rpvRatio] = getRpv(res(clu == cluid), ref);
    plotFetRpvHist(fetclu, rpv, 'cluid', cluid, 'grpid', grpid,...
        'rmvSpks', rmvSpks, 'rpvRatioOrig', rpvRatioOrig,...
        'rpvRatio', rpvRatio, 'visible', 'on', 'saveFig', false)

    % check fr
    uclu = unique(clu);
    cnt = 1;
    for iclu = 1 : length(uclu)
        if uclu(iclu) == 0 || uclu(iclu) == 1
            continue
        end
        
        cluidx = clu == uclu(iclu);
        spktimes{cnt} = res(cluidx) / fs;
        cnt = cnt + 1;
    end
    fr = calc_fr(spktimes, 'basepath', basepath,...
    'graphics', false, 'binsize', 60, 'saveVar', false,...
    'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf]);   
    fh = figure;
    plot(fr.tstamps / 60 / 60, mean(fr.strd))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over groups and clus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for igrp = 1 : length(grps)
    
    grpid = grps(igrp);
    grpchans = spkgrp{grpid};
    fprintf('\nworking on spkgrp #%d\n', grpid)
    
    % load data from neurosuite files
    res = loadNS('datatype', 'res', 'session', session, 'grpid', grpid);
    fet = loadNS('datatype', 'fet', 'session', session, 'grpid', grpid);
    clu = loadNS('datatype', 'clu', 'session', session, 'grpid', grpid);
    uclu = unique(clu);
    
    % here can write conditions to clean only specific clusters, for
    % example only if have >X nspks
    sclu = uclu;
    
    % rmv spks from all clusters until criterion is reached
    clean2criterion(fet, res, ref, rpvCrt, rmvLim, sclu, grpid, true)
    
    % save new clu
    saveNS(clu, 'datatype', 'clu', 'session', session, 'grpid', grpid);
end
end
% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% go over fets, calc cdf diff and remove spikes from clu
function [cluidx, rmvSpks] = cleanClu(fet, res, cluid, ref, cdfThr, rmvLim)

% OUTPUT:
%   cluidx      updated indices for cluid
%   rmvSpks     mat nfet x 3. 1st column - fet value from which spikes
%               were removed on the left. 2nd column - fet value
%               from which spikes were removed on the right. 3rd
%               column - cumsum number of spikes removed per fet

global clu
cluidx = find(clu == cluid);
fetclu = fet(cluidx, :);
nspks = length(cluidx);
nfet = size(fet, 2);
[rpv, rpvRatio] = getRpv(res(cluidx), ref);
spkRmvLimFet = floor(nspks * rmvLim(1));
spkRmvLimTot = floor(nspks * rmvLim(2));

rmvSpkIdx = [];
rmvSpks = nan(nfet, 3);
for ifet = 1 : nfet
    
    % find bin where number of spikes reaches percent limit
    [spkcount, binedges] = histcounts(fetclu(:, ifet), 100, 'Normalization', 'count');
    [rpvcount] = histcounts(fetclu(rpv, ifet), 'BinEdges', binedges, 'Normalization', 'count');
    leftSpk = cumsum(spkcount);
    rightSpk = fliplr(cumsum(spkcount, 'reverse'));
    leftRpv = cumsum(rpvcount);
    rightRpv = fliplr(cumsum(rpvcount, 'reverse'));
    [~, limLeft] = min(abs(spkRmvLimFet - leftSpk));
    [~, limRight] = min(abs(spkRmvLimFet - rightSpk));
    
    % calculate difference between cdfs of spks and rpvs
    % left side
    spkcdf = leftSpk(1 : limLeft - 1) / nspks;
    rpvcdf = leftRpv(1 : limLeft - 1) / length(rpv);
    [cdfDiff, rmBinIdx] = max(rpvcdf - spkcdf);
    if cdfDiff > cdfThr        % remove only if
        rmvSpks(ifet, 1) = binedges(rmBinIdx);
        rmvSpkIdx = [rmvSpkIdx; find(fetclu(:, ifet) < rmvSpks(ifet, 1))];
    end
    
    % right side
    spkcdf = rightSpk(1 : limRight) / nspks;
    rpvcdf = rightRpv(1 : limRight) / length(rpv);
    [cdfDiff, rmBinIdx] = max(rpvcdf - spkcdf);
    if cdfDiff > cdfThr
        rmvSpks(ifet, 2) = binedges(end - rmBinIdx);
        rmvSpkIdx = [rmvSpkIdx; find(fetclu(:, ifet) > rmvSpks(ifet, 2))];
    end
    rmvSpkIdx = unique(rmvSpkIdx);
    rmvSpks(ifet, 3) = length(rmvSpkIdx);
    
    % check if total number of spikes removed is greater than upper
    % limit
    if rmvSpks(ifet, 3) > spkRmvLimTot
        break
    end
end

% remove spks
clu(cluidx(rmvSpkIdx)) = 0;     % from clu
cluidx(rmvSpkIdx) = [];         % from cluidx
[~, rpvRatioNew] = getRpv(res(cluidx), ref);

% print to screen
fprintf('clu #%d: removed %d (%d%%) spikes; rpv ratio %.2f -> %.2f\n',...
    cluid, rmvSpks(end, 3), round(rmvSpks(end, 3) / nspks  * 100), rpvRatio, rpvRatioNew)

end

% -------------------------------------------------------------------------
% clean all clusters to rpv criteria
function clean2criterion(fet, res, ref, rpvCrt, rmvLim, sclu, grpid, graphics)
global clu
uclu = unique(clu);
if isempty(sclu)
    sclu = uclu;
end
for iclu = 1 : length(sclu)
    cluid = sclu(iclu);
    if cluid == 0 || cluid == 1
        continue
    end
    cluidx = find(clu == cluid);
    nspks = length(cluidx);
    fetcluOrig = fet(cluidx, :);
    
    % check rpv
    [rpv, rpvRatio] = getRpv(res(cluidx), ref);
    rpvRatioOrig = rpvRatio;
    if rpvRatio < rpvCrt
        continue
    end
    
    stepsize = 0.005;
    rmvLimTemp = [stepsize rmvLim(2)];
    nrmvSpk = 0;
    rmvSpks = zeros(size(fet, 2), 3);
    rmvSpks(:, 1) = -Inf;
    rmvSpks(:, 2) = Inf;
    while (nrmvSpk / nspks) < rmvLim(2) &&...           % stop if removed too many spikes
            rpvRatio > rpvCrt                           % stop if reached criterion
        while (nrmvSpk / nspks) < rmvLim(2)
            [cluidx, rmvSpksTemp] = cleanClu(fet, res, cluid, ref, 0.0001, rmvLimTemp);
            
            % update feature point of removal and nrmvSpk
            rowidx = rmvSpksTemp ~= 0;
            % left
            rmvSpks(rowidx(:, 1), 1) =...
                max(rmvSpksTemp(rowidx(:, 1), 1), rmvSpks(rowidx(:, 1), 1));
            % right
            rmvSpks(rowidx(:, 2), 2) =...
                min(rmvSpksTemp(rowidx(:, 2), 2), rmvSpks(rowidx(:, 2), 2));
            % nspks per feature
            rmvSpks(:, 3) = rmvSpks(:, 3) + rmvSpksTemp(:, 3);
            
            % stop if rmvLim no longer effective
            if round(rmvSpksTemp(end, 3) / nspks) * 100 < 1
                break
            end
        end
        [~, rpvRatio] = getRpv(res(cluidx), ref);
        rmvLimTemp(1) = rmvLimTemp(1) + stepsize;
    end
    
    % print to screen
    fprintf('clu #%d: removed total of %d (%d%%) spikes; rpv ratio %.2f -> %.2f\n\n',...
        cluid, rmvSpks(end, 3), round(rmvSpks(end, 3) / nspks  * 100), rpvRatioOrig, rpvRatio)
        
    % plot histogram of features for spks and rpvs
    if graphics
        plotFetRpvHist(fetcluOrig, rpv, 'cluid', cluid, 'grpid', grpid,...
            'rmvSpks', rmvSpks, 'rpvRatioOrig', rpvRatioOrig,...
            'rpvRatio', rpvRatio, 'visible', 'off', 'saveFig', true)
    end
end
end

% -------------------------------------------------------------------------
% plot histogram of features for spks and rpvs
function plotFetRpvHist(fetclu, rpv, varargin)

p = inputParser;
addOptional(p, 'cluid', 0, @isnumeric);
addOptional(p, 'grpid', 0, @isnumeric);
addOptional(p, 'rmvSpks', [], @isnumeric);
addOptional(p, 'rpvRatioOrig', 0, @isnumeric);
addOptional(p, 'rpvRatio', 0, @isnumeric);
addOptional(p, 'visible', 'off', @ischar);
addOptional(p, 'saveFig', false, @islogical);

parse(p, varargin{:})
cluid        = p.Results.cluid;
grpid        = p.Results.grpid;
rmvSpks      = p.Results.rmvSpks;
rpvRatioOrig = p.Results.rpvRatioOrig;
rpvRatio     = p.Results.rpvRatio;
visible      = p.Results.visible;
saveFig      = p.Results.saveFig;

global clu
nfet = size(fetclu, 2);
nspks = size(fetclu, 1);
[nsub] = numSubplots(nfet);
nbins = 100;
normmode = 'probability';
figure('units','normalized','outerposition',[0 0 1 1], 'visible', visible)
for ifet = 1 : nfet
    subplot(nsub(1), nsub(2), ifet)
    histogram(fetclu(:, ifet), nbins,...
        'FaceAlpha', 0.4, 'LineStyle', 'none', 'Normalization', normmode)
    hold on
    histogram(fetclu(rpv, ifet), nbins,...
        'FaceAlpha', 0.4, 'LineStyle', 'none', 'Normalization', normmode)
    yLim = ylim;
    if rmvSpks(ifet, 1) ~= -Inf
        plot([rmvSpks(ifet, 1) rmvSpks(ifet, 1)], yLim, '--k', 'LineWidth', 3)
    end
    if rmvSpks(ifet, 2) ~= Inf
        plot([rmvSpks(ifet, 2) rmvSpks(ifet, 2)], yLim, '--k', 'LineWidth', 3)
    end
    title(sprintf('fet #%d\nremoved %d spikes', ifet, rmvSpks(ifet, 3)))
    if ifet == 1
        legend({'All spks', 'RPVs'})
    end
end
suptitle(sprintf('T#%d clu#%d: removed %d (%d%%) spikes\n rpv ratio %.2f -> %.2f\n\n',...
    grpid, cluid, rmvSpks(end, 3), round(rmvSpks(end, 3) / nspks  * 100),...
    rpvRatioOrig, rpvRatio))

if saveFig
    mkdir(fullfile('graphics', 'clusterClean'))
    figname = fullfile('graphics', 'clusterClean', sprintf('T%d_clu%d', grpid, cluid));
    export_fig(figname)
end

end

% -------------------------------------------------------------------------
% print to screen isolation distance for all clusters in grp
function distOfClus = displayCluDist(fet, sclu)
global clu
uclu = unique(clu);
nfet = size(fet, 2);
fetDist = [fet, clu];
fetDist(clu <= 1, :) = [];
lRat = zeros(length(uclu), 1);
iDist = zeros(length(uclu), 1);
if isempty(sclu)
    sclu = uclu;
end
for iclu = 1 : length(sclu)
    cluid = sclu(iclu);
    if cluid == 0 || cluid == 1
        continue
    end
    vecidx = find(uclu == sclu(iclu));
    [lRat(vecidx, 1), iDist(vecidx, 1), ~] =...
        cluDist(fetDist(:, 1 : nfet), find((fetDist(:, end) == cluid)));
end
distOfClus = [uclu, iDist, lRat];
end

% -------------------------------------------------------------------------
% print to screen rpv ratio for all clusters in grp
function rpvRatioAll = displayRPV(res, ref)
global clu
uclu = unique(clu);
rpvRatioAll = zeros(length(uclu), 2);
rpvRatioAll(:, 1) = uclu;
for iclu = 1 : length(uclu)
    cluid = uclu(iclu);
    if cluid == 0 || cluid == 1
        continue
    end
    cluidx = find(clu == cluid);
    npks = length(cluidx);
    [~, rpvRatioAll(iclu, 2)] = getRpv(res(cluidx), ref);
end
end

% -------------------------------------------------------------------------
% find rpv and calc rpv ratio
function [rpv, rpvRatio] = getRpv(spktimes, ref)
rpv = find(diff([0; spktimes]) < ref);
nspks = length(spktimes);
rpvRatio = length(rpv) / nspks * 100;
end

% -------------------------------------------------------------------------
% additional subroutines

% alternative for cleaning cluster; compute pca only on cluster spikes and
% attemp at gmm to separate rvds
% recalculate pca
% clear fetnew
% for ichan = 1 : length(grpchans)
%     [~, pcFeat] = pca(permute(spk(ichan, :, cluidx), [3, 2, 1]), 'NumComponents', 3);
%     fetnew(:, ichan * 3 - 2 : ichan * 3) = pcFeat;
% end
%
% rpvgrp = ones(nspks, 1);
% rpvgrp(rpv) = 2;
% figure
% [nsub] = numSubplots(nfet / 2);
% kount = 1;
% for ifet = 1 : 2 : nfet
%     subplot(nsub(1), nsub(2), kount)
%     gscatter(fetnew(:, ifet), fetnew(:, ifet + 1), rpvgrp)
%     xlabel(sprintf('fet #%d', ifet))
%     ylabel(sprintf('fet #%d', ifet + 1))
%     kount = kount + 1;
% end
%
% % plot ccg of specific units
% % CCG
% binSize = 0.0001; dur = 0.6; % low res
% iclus = [83, 116];
% for icluster = 1 : length(iclus)
%     cluidx = clu(clu == iclus(icluster));
%     spkclutimes{icluster} = res(cluidx)
% [ccg, t] = CCG({xx.times{:}}, [], 'duration', dur, 'binSize', binSize);
% u = 20;
% plotCCG('ccg', ccg(:, u, u), 't', t, 'basepath', basepath,...
%     'saveFig', false, 'c', {'k'}, 'u', spikes.UID(u));