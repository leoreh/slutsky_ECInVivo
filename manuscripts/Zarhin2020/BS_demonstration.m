% BS_demonstration

% select mouse
i = 1;      % grp
j = 25;     % mouse

% params
graphics = true;
loaddata = true;
saveFig = false;

basepath{1} = 'E:\Data\Others\DZ\IIS\WT';
basepath{2} = 'E:\Data\Others\DZ\IIS\APPPS1';
basepath{3} = 'E:\Data\Others\DZ\IIS\APPKi';
basepath{4} = 'E:\Data\Others\DZ\IIS\FADx5';
cd(basepath{i})
filename = dir('*.abf');
files = natsort({filename.name});
nfiles = 1 : length(files);
[~, grpname] = fileparts(basepath{i});
[~, basename] = fileparts(files{nfiles(j)});

fs = 1250;
binsize = (2 ^ nextpow2(5 * fs));
smf = 7;
marg = 0.05;
thr = [0 0];
ch = 1;

if loaddata
    % load lfp
%     lfp = getLFP('basepath', basepath{i}, 'ch', ch, 'chavg', {},...
%         'fs', 1250, 'interval', [0 inf], 'extension', 'abf', 'pli', true,...
%         'savevar', false, 'force', false, 'basename', basename);
%     sig = double(lfp.data(:, ch));
    
    minDur = 2;                         % minimum event duration [s]
    maxDur = 3000000;                   % maximum event duration [s]
    interDur = 4;                       % minimum time between events [bins]
    binsize = 2 * fs;                 % [s] * [fs] = [samples]
    nclust = 2;                         % two clusters: burst and suppression
    
    ibins = [1 : binsize : length(sig)];
    ibins(end) = length(sig);
    
    % divide signal to bins
    sigre = sig(1 : end - (mod(length(sig), binsize) + binsize));
    sigmat = reshape(sigre, binsize, (floor(length(sig) / binsize) - 1));
    
    % last bin
    siglastbin = sig(length(sigre) + 1 : length(sig));
    
    sumvec = sum(abs(sigmat));
    sumvec = [sumvec, sum(abs(siglastbin))];
    maxvec = max(abs(sigmat));
    maxvec = [maxvec, max(abs(siglastbin))];
    stdvec = std(sigmat);
    stdvec = [stdvec, std(siglastbin)];
    % concatenate and lognorm
    varmat = [stdvec; maxvec; sumvec];
    if size(varmat, 1) < size(varmat, 2)
        varmat = varmat';
    end
    varmat = log10(varmat);
    
    options = statset('MaxIter', 50);
    gm = fitgmdist(varmat, nclust,...
        'options', options, 'Start', 'plus', 'Replicates', 50);
    % cluster
    [gi, ~, ~, ~, ~] = cluster(gm, varmat);
    if mode(bs.gm.mu(1, :) > bs.gm.mu(2, :))
        gi(gi == 2) = 0;
        gi(gi == 1) = 2;
        gi(gi == 0) = 1;
    end
end

if graphics
    fh = figure('Visible', 'on');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    suptitle(basename)
    
    k = 1;
    for i = 1 : 3
        subplot(6, 2, [k k + 2])
        h = histogram(varmat(:, i), 100, 'Normalization', 'Probability');
        h.EdgeColor = 'none';
        h.FaceColor = 'k';
        h.FaceAlpha = 1;
        ylabel('Probability [%]')
        set(gca, 'TickLength', [0 0])
        box off
        k = k + 4;
        switch i
            case 1; xlabel('std [log10(\sigma)')
            case 2; xlabel('max [log10(V)]')
            case 3; xlabel('sum [log(mV)]')
        end
    end
    
    % clustering - first two variables
    subplot(6, 2, [2 : 2 : 6])
    options = statset('MaxIter', 50);
    gm = fitgmdist(varmat(:, [1, 2]), nclust,...
        'options', options, 'Start', 'plus', 'Replicates', 50);
    gscatter(varmat(:, 1), varmat(:, 2), gi, 'rk', '.', 3);
    axis tight
    hold on
    gmPDF = @(x1, x2)reshape(pdf(gm, [x1(:) x2(:)]), size(x1));
    ax = gca;
    fcontour(gmPDF, [ax.XLim ax.YLim], 'HandleVisibility', 'off')
    xlabel(['std [log10(\sigma)]'])
    ylabel(['max [log10(V)]'])
    legend off
    set(gca, 'TickLength', [0 0])
    box off
    
    subplot(6, 2, [8 : 2 : 12])
    options = statset('MaxIter', 50);
    gm = fitgmdist(varmat(:, [2, 3]), nclust,...
        'options', options, 'Start', 'plus', 'Replicates', 50);
    gscatter(varmat(:, 2), varmat(:, 3), gi, 'rk', '.', 3);
    axis tight
    hold on
    gmPDF = @(x1, x2)reshape(pdf(gm, [x1(:) x2(:)]), size(x1));
    ax = gca;
    fcontour(gmPDF, [ax.XLim ax.YLim], 'HandleVisibility', 'off')
    xlabel(['max [log10(V)]'])
    ylabel('sum [log(mV)]')
    legend off
    set(gca, 'TickLength', [0 0])
    box off
    
    if saveFig
        figname = 'BS demonstration';
        export_fig(figname, '-tif', '-transparent')
        % savePdf(figname, basepath, ff)
    end
    
end