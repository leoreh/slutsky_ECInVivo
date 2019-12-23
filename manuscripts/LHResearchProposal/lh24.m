%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and arrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

force = false;
u = ([20, 27]);
%%% CCH

    % high res
% binSize = 0.0001; % [s]
% dur = 0.01;
% [ccg2, t2] = CCG({spikes.times{u}}, [], 'duration', dur, 'binSize', binSize);
% 
% % low res
% binSize = 0.001; % [s]
% dur = 0.06;
% [ccg1, t1] = CCG({spikes.times{u}}, [], 'duration', dur, 'binSize', binSize);


if force
    basepath = 'E:\Data\Dat\lh24';
    cd(basepath)
    load('lh24.Raw1.Info.mat')
    load('lh24.fet.mat')
    
    %%% clip recording to relavent time
    interval = cumsum(info.blockduration(1 : 2));
    interval = [1, interval(2)];
    
    % 19 = PYR, 17 = INT
    
    %%% load spikes
    spikes = getSpikes('basepath', basepath, 'saveMat', false, 'noPrompts', true);
    lfp = getLFP('basepath', basepath, 'chans', [1 : 16], 'chavg', {}, 'fs', 1250,...
        'interval', [0 inf], 'savevar', false);
    
    %%% burst suppression
    % idx = 50 * 60 * lfp.fs : 60 * 60 * lfp.fs;
    % sig = lfp.data(idx, 16);
    % getBS('sig', sig, 'fs', lfp.fs)
    
    % Isolation distance. When including all cells, only 7 clusters were SU.
    % Removing pre-determined MU increased SU count to 10.
    
    %     mu = [3, 4, 6, 10, 13, 16, 18, 19, 20, 23, 25, 26];
    %     spikes = cluVal(spikes, 'basepath', basepath, 'saveVar', false,...
    %         'saveFig', true, 'force', true, 'mu', mu, 'graphics', false);
    %
    
       
    %%% firing rate
    fr = FR(spikes.times, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
        'binsize', 300, 'saveVar', true, 'smet', 'MA', 'winBL', [1, 30 * 60], 'select', {'thr'});
    [nunits, nbins] = size(fr.strd);
    tFR = ([1 : nbins] / (60 / fr.binsize) / 60);
    
    idxControl = [1 : find(tFR * 60 > 30, 1)];
    idxPsam = [find(tFR * 60 > 45, 1) : find(tFR * 60 > 75, 1)];
    idxReturn = [find(tFR * 60 > 240, 1) : find(tFR * 60 > 270, 1)];
    
    select = [2];
    lns = cumsum(info.blockduration / 60 / 60);
    lns = [lns([select - 1])];
    lbs = {'uPSEM'};
    
    %%% INT vs PYR
    % units = [1 : 5, 9 : 13, 15 : 17, 19 : 22, 24 : 27];
    CellClass = cellClass('waves', cat(1, spikes.rawWaveform{:})', 'u', u,...
        'man', false, 'fs', spikes.samplingRate, 'mfr', fr.mfr(:),...
        'saveVar', false, 'graphics', false);
    
    % high res
binSize = 0.0002; % [s]
dur = 0.01;
[ccg2, t2] = CCG({spikes.times{u}}, [], 'duration', dur, 'binSize', binSize);

% low res
binSize = 0.001; % [s]
dur = 0.06;
[ccg1, t1] = CCG({spikes.times{u}}, [], 'duration', dur, 'binSize', binSize);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
f = figure;

c = {'b', 'r'};

% exp. design
subplot(4, 4, 1)
title('Exp. Design')
box off
yticks([])
xticks([])

% histology
subplot(4, 4, 5)
title('Histology')
box off
yticks([])
xticks([])

% waveform u1
subplot(4, 4, 2)
plotWaveform('avgwv', spikes.avgWaveform{u(1)}, 'stdwv', spikes.stdWaveform{u(1)},...
    'orient', 'vert', 'fs', spikes.samplingRate, 'c', c{1});
title(sprintf('u%d', u(1)))

% waveform u2
subplot(4, 4, 6)
plotWaveform('avgwv', spikes.avgWaveform{u(2)}, 'stdwv', spikes.stdWaveform{u(2)},...
    'orient', 'vert', 'fs', spikes.samplingRate, 'c', c{2});
title(sprintf('u%d', u(2)))

% ACC u1
subplot(4, 4, 3)
plotCCG('ccg', ccg1(:, 1, 1), 't', t1, 'basepath', basepath,...
    'saveFig', false, 'c', {c{1}});
xlabel('')
% ylabel('')
% yticks([])
xticks([])

% ACC u2
subplot(4, 4, 7)
plotCCG('ccg', ccg1(:, 2, 2), 't', t1, 'basepath', basepath,...
    'saveFig', false, 'c', {c{2}});
xticks([-30 0 30])

% ccg
subplot(4, 4, 4)
plotCCG('ccg', ccg2(:, 1, 2), 't', t2, 'basepath', basepath,...
    'saveFig', false, 'c', {'k'});

% INT vs. PYR
subplot(4, 4, 8)
mfr = rescale(fr.mfr, 10, 50);
s = scatter(CellClass.tp, CellClass.spkw, mfr, 'filled', 'k');
hold on
scatter(CellClass.tp(u(1)), CellClass.spkw(u(1)), mfr(u(1)), 'filled', c{1})
scatter(CellClass.tp(u(2)), CellClass.spkw(u(2)), mfr(u(2)), 'filled', c{2})
xlabel('trough-to-peak [ms]')
ylabel('spike width [ms]')
set(gca, 'TickLength', [0 0])

% FR units
subplot(4, 4, 9 : 11)
hold on
for i = 1 : nunits
    p = plot(tFR, log10(fr.strd(i, :)));
    if i == u(1)
        p.LineWidth = 3;
        p.Color = c{1};
    elseif i == u(2)
        p.LineWidth = 3;
        p.Color = c{2};
    else
        p.Color = 'k';
    end
end
axis tight
ylabel('Firing Rate [log(Hz)]')
if ~isempty(lns)
    nlines = length(lns);
    x = repmat(lns, 2, 1);
    y = repmat(ylim, nlines, 1);
    plot(x, y' ,'--k');
    if ~isempty(lbs)
        text(x(1, :), y(1, 2)*ones(nlines, 1), lbs)
    end
end
hold on
y = ylim;
z = idxControl * fr.binsize / 60 / 60;
p = patch([z(1) z(end) z(end) z(1)], [y(1) y(1) y(2) y(2)], 'g');
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';
z = idxPsam * fr.binsize / 60 / 60;
p = patch([z(1) z(end) z(end) z(1)], [y(1) y(1) y(2) y(2)], 'm');
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';
z = idxReturn * fr.binsize / 60 / 60;
p = patch([z(1) z(end) z(end) z(1)], [y(1) y(1) y(2) y(2)], 'y');
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';
xlim([0 5])
xticks([0 : 5])
set(gca, 'TickLength', [0 0])

% FR avg
subplot(4, 4, 13 : 15)
plotFRtime('fr', fr, 'units', false, 'lns', lns, 'lbs', lbs,...
    'avg', true, 'raster', false, 'saveFig', false);
title('')
hold on
y = ylim;
z = idxControl * fr.binsize / 60 / 60;
p = patch([z(1) z(end) z(end) z(1)], [y(1) y(1) y(2) y(2)], 'g');
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';
z = idxPsam * fr.binsize / 60 / 60;
p = patch([z(1) z(end) z(end) z(1)], [y(1) y(1) y(2) y(2)], 'm');
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';
z = idxReturn * fr.binsize / 60 / 60;
p = patch([z(1) z(end) z(end) z(1)], [y(1) y(1) y(2) y(2)], 'y');
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';
xlim([0 5])
xticks([0 : 5])
set(gca, 'TickLength', [0 0])

% change in norm. fr
subplot(4, 4, 16)
FRcontrol = mean(fr.norm(:, idxControl), 2);
FRpsam = mean(fr.norm(:, idxPsam), 2);
FRreturn = mean(fr.norm(:, idxReturn), 2);
mat = [FRcontrol FRpsam FRreturn];
xpoints = 1 : 3;
line(xpoints, mat, 'Color', [0.5 0.5 0.5])
hold on
line(1 : 3, mat(15, :), 'Color', c{1}, 'LineWidth', 3)
line(1 : 3, mat(22, :), 'Color', c{2}, 'LineWidth', 3)
line(xpoints, [mean(mat)], 'Color', 'k', 'LineWidth', 5)
axis tight
xlim([xpoints(1) - 0.2, xpoints(end) + 0.2])
ylabel('Norm. Firing Rate')
ax = gca;
ax.XTick = xpoints;
ax.XTickLabel = {'BL', 'uPSEM', 'WO'};
xlabel('')
% ylim([0 2])

% FR distribution (Hz, not norm)
subplot(4, 4, 12)
FRcontrol = mean(fr.strd(:, idxControl), 2);
FRpsam = mean(fr.strd(:, idxPsam), 2);
FRreturn = mean(fr.strd(:, idxReturn), 2);
h = histogram(log10(FRcontrol), 5);
h.EdgeColor = 'none';
h.FaceColor = 'g';
h.FaceAlpha = 0.3;
hold on
h = histogram(log10(FRpsam), 5);
h.EdgeColor = 'none';
h.FaceColor = 'm';
h.FaceAlpha = 0.3;
h = histogram(log10(FRreturn), 5);
h.EdgeColor = 'none';
h.FaceColor = 'y';
h.FaceAlpha = 0.3;
box off
axis tight
xlabel('Firing Rate [log(Hz)]')
ylabel('Number of Units')
set(gca, 'TickLength', [0 0])




%%% save
% filename = 'lh24';
% savePdf(filename, basepath, f)



