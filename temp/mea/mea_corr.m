

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = pwd;
cd(basepath)
[~, basename] = fileparts(basepath);

frfile = fullfile(basepath, [basename, '.fr.mat']);
meafile = fullfile(basepath, [basename, '.mea.mat']);
load(frfile)
load(meafile)

% unit stats
spktimes = mea.spktimes(fr.stable);
lastspike = max(cellfun(@max, spktimes, 'uni', true));     % last spike [s]
nspks = cellfun(@length, spktimes, 'uni', true);
[~, maxnspks] = max(nspks);

% limit to high fr neurons
% sunits = nspks > prctile(nspks, 50);
% spktimes = spktimes(sunits);

nunits = length(spktimes);
npairs = nunits * nunits - nunits;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc corr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

countbinsize = 0.1;     % [s]
corrbinsize = 60 * 60;   % [s]

countbins = [0 : countbinsize : lastspike];
corrbins = [0 : corrbinsize : lastspike];
nbins = length(corrbins) - 1;

spkcounts = cellfun(@(x) histcounts(x, countbins), spktimes, 'uni', false);
spkcounts = cell2mat(spkcounts');
maxcounts = max(spkcounts, [], 2);

% limit to hr units in baseline
countidx = countbins >= 2 * 60 &...
    countbins < 60 * 60;


cc = zeros(nbins, nunits, nunits);
for icorr = 1 : nbins
    countidx = countbins >= corrbins(icorr) &...
        countbins < corrbins(icorr + 1);
    cc(icorr, :, :) = corrcoef(spkcounts(:, countidx)', 'rows', 'all');
end
% set diagonal to maximum to avoid saturation
for ibin = 1 : nbins
    tmp = squeeze(cc(ibin, :, :));
    tmp(1 : nunits + 1 : numel(tmp)) = nan;
    tmp(1 : nunits + 1 : numel(tmp)) = max(tmp, [], 'all');
    cc(ibin, :, :) = tmp;
end

% sort
[~, maxcorrunits] = max(sum(squeeze(cc(1, :, :)), 'omitnan'));
scc = cc;
[~, sortidx] = sort(squeeze(cc(1, maxcorrunits, :)), 'descend');
for icorr = 1 : nbins
    scc(icorr, :, :) = scc(icorr, sortidx, sortidx);
end

for ibin = 1 : nbins
    corrstrength(ibin) = sum(cc(ibin, :, :), 'all');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
nbins = size(cc, 1);
cLimit = [min(cc, [], 'all'), max(cc, [], 'all')];
cLimit = [prctile(cc(:), 5), prctile(cc(:), 95)];

% firing rate
subplot(3, ceil(nbins / 2), [1 : ceil(nbins / 2)])
frMat = fr.strd(fr.stable, :);
plot(fr.tstamps / 60 / 60, mean(frMat, 1, 'omitnan'), 'k', 'LineWidth', 2)
hold on
axis tight
xlabel('Time [~h]')
ylabel('Firing rate [Hz]')
set(gca, 'box', 'off')
title(basename)
legend(sprintf('nunits = %d', size(frMat, 1)))

for ibin = 1 : nbins - 1
    subplot(3, ceil(nbins / 2), ceil(nbins / 2) + ibin)
    imagesc([1, nunits], [1, nunits], squeeze(scc(ibin, :, :)), cLimit)
    title(sprintf('%dhr', ibin))
    xticks([])
    yticks([])
    if ibin == 1
        xlabel('Units')
        ylabel('Units')
    end
    if ibin == nbins - 1
        subplot(3, ceil(nbins / 2), ceil(nbins / 2) + ibin + 1)
        imagesc([cLimit cLimit])
        set(gca, 'YColor', 'w', 'XColor', 'w')
        colorbar('eastoutside')
    end    
end
