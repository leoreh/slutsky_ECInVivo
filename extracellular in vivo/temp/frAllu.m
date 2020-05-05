close all

force = false;

path{1} = 'H:\Data\Dat\lh50\lh50_200426\080458_e3r1-1';
path{2} = 'H:\Data\Dat\lh50\lh50_200426\110419_e8r1-1';
path{3} = 'H:\Data\Dat\lh50\lh50_200426\120413_e10r1-1';
path{4} = 'H:\Data\Dat\lh50\lh50_200427\080449_e1r1-1';
path{5} = 'H:\Data\Dat\lh50\lh50_200427\090420_e3r1-1';
path{6} = 'H:\Data\Dat\lh50\lh50_200428\080422_e1r1-1';
path{7} = 'H:\Data\Dat\lh50\lh50_200428\090410_e3r1-1';

if force
    for i = 1 : length(path)
        basepath = path{i};
        cd(basepath)
        basename = bz_BasenameFromBasepath(basepath);
        
        accname = [basename '.acceleration.mat'];
        load(accname)
        
        spikes = getSpikes('basepath', basepath, 'saveMat', true,...
            'noPrompts', true, 'forceL', false);
        
        dur(i) = acc.tstamps(end);
        y{i} = {sort(vertcat(spikes.times{:}))};
        
        for j = 1 : length(unique(spikes.shankID))
            x{j} = sort(vertcat(spikes.times{spikes.shankID == j}));
        end
        
        [fr.strd, ~, fr.tstamps] = calcFR(x, 'binsize', 60,...
            'winCalc', [1 Inf], 'smet', 'MA');
        
        f{i} = fr;
        a{i} = acc;
        b{i} = basename;
        
    end
    
%     y = [y{:}];
%     for i = 1 : length(y)
%     nspk(i) = length(y{i});
%     end
%     nspk ./ dur
    
    path{1} = 'H:\Data\Dat\lh50\lh50_200426\100442_e6r1-4';
    path{2} = 'H:\Data\Dat\lh50\lh50_200426\120451_e9r1-5';
    path{3} = 'H:\Data\Dat\lh50\lh50_200426\130425_e11r1-5';
    path{4} = 'H:\Data\Dat\lh50\lh50_200427\090428_e2r1-5';
    path{5} = 'H:\Data\Dat\lh50\lh50_200427\100440_e4r1-5';
    path{6} = 'H:\Data\Dat\lh50\lh50_200428\090403_e2r1-5';
    path{7} = 'H:\Data\Dat\lh50\lh50_200428\100405_e4r1-5';
    
    for i = 1 : length(path)
        basepath = path{i};
        cd(basepath)
        basename = bz_BasenameFromBasepath(basepath);
        
        fepspname = [basename '.fepsp.mat'];
        load(fepspname)
        
        fe{i} = fepsp;
    end
    path{1} = 'H:\Data\Dat\lh50\lh50_200426\100442_e6r1-4';
    path{2} = 'H:\Data\Dat\lh50\lh50_200426\120451_e9r1-5';
    path{3} = 'H:\Data\Dat\lh50\lh50_200426\130425_e11r1-5';
    path{4} = 'H:\Data\Dat\lh50\lh50_200427\090428_e2r1-5';
    path{5} = 'H:\Data\Dat\lh50\lh50_200427\100440_e4r1-5';
    path{6} = 'H:\Data\Dat\lh50\lh50_200428\090403_e2r1-5';
    path{7} = 'H:\Data\Dat\lh50\lh50_200428\100405_e4r1-5';
    
    for i = 1 : length(path)
        basepath = path{i};
        cd(basepath)
        basename = bz_BasenameFromBasepath(basepath);
        
        fepspname = [basename '.fepsp.mat'];
        load(fepspname)
        
        fe{i} = fepsp;
        be{i} = basename;
    end
end

figure
for i = 1 : length(path)
    % field
    subplot(3, length(path), i)
    plot(fe{i}.amp')
    if i == 1
        ylabel('Amplitude [mV]')
        set(gca, 'TickLength', [0 0], 'XTickLabel', [],...
            'Color', 'none', 'XColor', 'none')
        legend(split(num2str(unique(spikes.shankID))))
    else
        set(gca, 'TickLength', [0 0], 'XTickLabel', [], 'YTickLabel', [],...
            'Color', 'none', 'XColor', 'none')
    end
    box off
    axis tight
    ylim([0 7])
    title(be{i}, 'Interpreter', 'none')
    
    % firing rate
    subplot(3, length(path), length(path) + i)
    plot(f{i}.tstamps / 60, f{i}.strd)
    if i == 1
        ylabel('Total Spike Rate')
        set(gca, 'TickLength', [0 0], 'XTickLabel', [],...
            'Color', 'none', 'XColor', 'none')
    else
        set(gca, 'TickLength', [0 0], 'XTickLabel', [], 'YTickLabel', [],...
            'Color', 'none', 'XColor', 'none')
    end
    box off
    axis tight
    ylim([0 300])
    title(b{i}, 'Interpreter', 'none')
    
    % acceleration
    subplot(3, length(path), length(path) * 2 + i)
    plot(a{i}.tband / 60, a{i}.pband)
    %     plot(a{i}.tstamps / a{i}.fs_orig / 60, a{i}.mag)
    if i == 1
        ylabel('Acceleration Power')
        set(gca, 'TickLength', [0 0],...
            'Color', 'none')
    else
        set(gca, 'TickLength', [0 0], 'YTickLabel', [],...
            'Color', 'none')
    end
    box off
    hold on
    axis tight
    ylim([600 1800])
    Y = ylim;
    fill([a{i}.sleep fliplr(a{i}.sleep)]' / 60, [Y(1) Y(1) Y(2) Y(2)],...
        'k', 'FaceAlpha', 0.25,  'EdgeAlpha', 0);
    xlabel('Time [m]')
    
    figname = 'LTP experiment';
    export_fig(figname, '-tif', '-transparent')
end


