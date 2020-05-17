path{1} = 'H:\Data\Dat\lh44\lh44_200208';
path{2} = 'H:\Data\Dat\lh46\lh46_200225a';
path{3} = 'H:\Data\Dat\lh46\lh46_200227a';
path{4} = 'H:\Data\Dat\lh49\lh49_200319';
path{5} = 'H:\Data\Dat\lh46\lh46_200225b';
path{6} = 'H:\Data\Dat\lh45\lh45_200209';
path{7} = 'H:\Data\Dat\Refaela\120520_Bsl';

force = true;
nsessions = length(path);

if force
    waves = [];
    mfr = [];
    spktimes = {};
    for i = 1 : nsessions
        basepath = path{i};
        cd(basepath)
        bname{i} = bz_BasenameFromBasepath(basepath);
       
        s{i} = getSpikes('basepath', basepath, 'saveMat', true,...
            'noPrompts', true, 'forceL', false);
        spktimes = [spktimes{:}, s{i}.times];
        nu(i) = length(s{i}.UID);
        
        waves = [waves, cat(1, s{i}.maxwv)'];
        
        fr{i} = FR(s{i}.times, 'basepath', basepath, 'graphics', false, 'saveFig', false,...
            'binsize', 60, 'saveVar', false, 'smet', 'MA');
        mfr = [mfr, mean(fr{i}.strd, 2)'];
    end
end

nunits = size(waves, 2);

% waveform parameters
cc = cellclass('waves', waves, 'saveVar', false,...
    'fs', s{1}.samplingRate, 'man', false); 

% burstiness
[b, bidx] = isbursty('spktimes', spktimes);
b(b == 0) = 2;
c = {'k', 'r'};

% ACG
binSize = 0.001; % [s]
dur = 0.06;
[ccg, t1] = CCG(spktimes, [], 'duration', dur, 'binSize', binSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examine separting matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbins = 70;
fh = figure;
subplot(2, 3, 1)
histogram(cc.tp, nbins)
xlabel('tp')
subplot(2, 3, 2)
histogram(cc.spkw, nbins)
xlabel('spkw')
subplot(2, 3, 3)
histogram(cc.hpk, nbins)
xlabel('hpk')
subplot(2, 3, 4)
histogram(cc.asym, nbins)
xlabel('asym')
subplot(2, 3, 5)
histogram(log10(mfr), nbins)
xlabel('log(mfr)')
subplot(2, 3, 6)
histogram(log10(bidx), nbins)
xlabel('log(burstiness)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% review individual neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smfr = rescale(mfr, 10, 30);

for i = 107 : 120
    fh = figure;
    suptitle(['clu ' num2str(i)])
    subplot(2, 2, 1)
    plot(waves(:, i));
    
    subplot(2, 2, 2)
    scatter(cc.tp, cc.spkw, smfr, 'filled', 'k')
    hold on
    scatter(cc.tp(i), cc.spkw(i), smfr(i), 'filled', 'r')
    xlabel('trough-to-peak [ms]')
    ylabel('spike width [ms]')
    set(gca, 'TickLength', [0 0])

    subplot(2, 2, 3)
    plotCCG('ccg', ccg(:, i, i), 't', t1, 'basepath', basepath,...
        'saveFig', false, 'c', {c{b(i)}});
    xlabel('')
    % ylabel('')
    % yticks([])
    xticks([-30 0 30])
    
    subplot(2, 2, 4)
    isi = diff(spktimes{i});
    bins = [0 : 0.001 : 0.05];
    histogram(isi, bins)
end