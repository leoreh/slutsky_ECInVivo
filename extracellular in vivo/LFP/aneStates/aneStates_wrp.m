
% this is a wrapper for the analysis of anesthesia states across mice and
% groups. directs the selected files to aneStates.m and arranges data in
% the struct as. within as the fields bsr, iisRate, and dband are cells of
% mats (one cell per grp, one row per mouse). spks and duration data are
% arranged in mats where each column is a grp and each row is a mouse.

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% group number refers to the order in the grppath cell.
grp = [1 : 4];
% if mouse empty will load all mice in group
mouse = [1 : 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path to data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grppath{1} = 'E:\Data\Field\IIS\WT';
grppath{2} = 'E:\Data\Field\IIS\APPPS1';
grppath{3} = 'E:\Data\Field\IIS\APPKi';
grppath{4} = 'E:\Data\Field\IIS\FADx5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forceA = false;
forceL = false;
saveFig = false;
graphics = false;
saveVar = false;
ch = 1;
smf = 7;        % smooth factor
fs = 1250;      % sampling frequency
binsize = 30;   % for BSR and delta [s]
maxmice = [];   % initialize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% send selected files to anestats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go over each grp
for j = grp
    
    basepath = grppath{j};
    cd(basepath)
    [~, grpname] = fileparts(basepath);
    as.grpname{j} = grpname;
    filename = dir('*.lfp.*');
    files = natsort({filename.name});
    nfiles = 1 : length(files);
    maxmice = max([maxmice, nfiles]);
    
    if isempty(mouse)
        mouse = 1 : nfiles;
    end
    
    % go over each mouse
    for i = mouse
        
        [~, basename] = fileparts(files{nfiles(i)});
        [~, basename] = fileparts(basename);
        
        % send to anestats
        [bs, iis, ep] = aneStates('ch', ch, 'basepath', basepath,...
            'basename', basename, 'graphics', graphics,...
            'saveVar', saveVar, 'saveFig', saveFig, 'forceA', forceA,...
            'binsize', binsize, 'smf', smf, 'thrMet', 1);
        
        % arrange group data
        as.mouse{j}{i} = basename;
        recDur{j}(i) = ep.recDur;
        deepDur{j}(i) = ep.deepDur;
        surDur{j}(i) = ep.surDur;
        nspks{j}(i) = ep.nspks;
        deep_nspks{j}(i) = ep.deep_nspks;
        sur_nspks{j}(i) = ep.sur_nspks;
        thr{j}(i) = iis.thr(2);
        sur_delta{j}(i) = ep.sur_delta;
        deep_delta{j}(i) = ep.deep_delta;
        iisRate{i} = iis.rate;
        bsr{i} = bs.bsr;
        dband{i} = ep.dband;
        t{i} = bs.cents / fs / 60;
        
        % after last recording is loaded, arrange in mat and save
        if i == mouse(end)
            
            [~, idx] = max(recDur{j});      % index to longest recording
            maxdur = length(bsr{idx});
            
            % bsr, iis and delta
            mat = cellfun(@(x)[x(:); NaN(maxdur-length(x), 1)], bsr,...
                'UniformOutput', false);
            as.bsr{j} = cell2mat(mat)';
            mat = cellfun(@(x)[x(:); NaN(maxdur-length(x), 1)], iisRate,...
                'UniformOutput', false);
            as.iisRate{j} = cell2mat(mat)';
            mat = cellfun(@(x)[x(:); NaN(maxdur-length(x), 1)], dband,...
                'UniformOutput', false);
            as.dband{j} = cell2mat(mat)';
            % timestamps
            mat = cellfun(@(x)[x(:); NaN(maxdur-length(x), 1)], t,...
                'UniformOutput', false);
            as.t{j} = cell2mat(mat)';
            as.t{j} = as.t{j}(idx, :);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange struct for prism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mat = cellfun(@(x)[x(:); NaN(maxmice-length(x), 1)], nspks,...
    'UniformOutput', false);
as.nspks = cell2mat(mat);

mat = cellfun(@(x)[x(:); NaN(maxmice-length(x), 1)], deep_nspks,...
    'UniformOutput', false);
as.deep_nspks = cell2mat(mat);

mat = cellfun(@(x)[x(:); NaN(maxmice-length(x), 1)], sur_nspks,...
    'UniformOutput', false);
as.sur_nspks = cell2mat(mat);

mat = cellfun(@(x)[x(:); NaN(maxmice-length(x), 1)], thr,...
    'UniformOutput', false);
as.thr = cell2mat(mat);

mat = cellfun(@(x)[x(:); NaN(maxmice-length(x), 1)], recDur,...
    'UniformOutput', false);
recDur = cell2mat(mat);
as.recDur = recDur / fs / 60;

mat = cellfun(@(x)[x(:); NaN(maxmice-length(x), 1)], deepDur,...
    'UniformOutput', false);
deepDur = cell2mat(mat);
as.deepDur = deepDur / fs / 60;

mat = cellfun(@(x)[x(:); NaN(maxmice-length(x), 1)], surDur,...
    'UniformOutput', false);
surDur = cell2mat(mat);
as.surDur = surDur / fs / 60;

as.surFraction = surDur ./ recDur;
as.deepFraction = deepDur ./ recDur;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if graphics
    
    recmean = cellfun(@mean, recDur);
    bsmean = cellfun(@mean, bsDur);
    bmean = cellfun(@mean, bDur);
    
    c = [1 0 0; 1 0 1; 0 0 1; 0 1 1];
    c2 = 'rmbc';
    gidx = [ones(1, length(bsDur{1})), ones(1, length(bsDur{2})) * 2,...
        ones(1, length(bsDur{3})) * 3, ones(1, length(bsDur{4})) * 4];
    
    figure
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    % short recordings
    subplot(3, 4, 1 : 2)
    hold on
    for i = [1, 2]
        stdshade(as{i}.iis, 0.1, c2(i), as{i}.t / fs / 60)
    end
    axis tight
    ylabel('IIS Rate [spikes / bin]')
    xlim([1 50])
    box off
    set(gca, 'TickLength', [0 0])
    title('Short Recordings')
    
    % long recordings
    subplot(3, 4, 5 : 6)
    hold on
    for i = [3, 4]
        stdshade(as{i}.iis, 0.1, c2(i), as{i}.t / fs / 60)
    end
    axis tight
    ylabel('IIS Rate [spikes / bin]')
    xlabel('Time [m]')
    xlim([1 100])
    box off
    set(gca, 'TickLength', [0 0])
    title('Long Recordings')
    
    % IIS total
    subplot(3, 4, 3)
    boxplot([nspks{:}] ./ ([recDur{:}] / fs / 60), gidx,...
        'BoxStyle', 'outline', 'Color', c2, 'notch', 'off')
    hold on
    gscatter(gidx, [nspks{:}] ./ ([recDur{:}] / fs / 60), gidx, c2)
    legend off
    xlabel('')
    ylabel('IIS Rate [spikes / bin]')
    xticklabels(as)
    title('IIS total')
    box off
    set(gca, 'TickLength', [0 0])
    
    % IIS in BS
    subplot(3, 4, 7)
    boxplot([nspksBS{:}] ./ ([bsDur{:}] / fs / 60), gidx,...
        'BoxStyle', 'outline', 'Color', c2, 'notch', 'off')
    hold on
    gscatter(gidx, [nspksBS{:}] ./ ([bsDur{:}] / fs / 60), gidx, c2)
    legend off
    xlabel('')
    ylabel('IIS Rate [spikes / bin]')
    xticklabels(as)
    title('IIS in BS')
    box off
    set(gca, 'TickLength', [0 0])
    
    % IIS in B
    subplot(3, 4, 11)
    boxplot([nspksB{:}] ./ ([bDur{:}] / fs / 60), gidx,...
        'BoxStyle', 'outline', 'Color', c2, 'notch', 'off')
    hold on
    gscatter(gidx, [nspksB{:}] ./ ([bDur{:}] / fs / 60), gidx, c2)
    legend off
    xlabel('')
    ylabel('IIS Rate [spikes / bin]')
    xticklabels(as)
    title('IIS in B')
    box off
    set(gca, 'TickLength', [0 0])
    
    % Duration total
    subplot(3, 4, 4)
    boxplot([recDur{:}] / fs / 60, gidx,...
        'BoxStyle', 'outline', 'Color', c2, 'notch', 'off')
    hold on
    gscatter(gidx, [recDur{:}] / fs / 60, gidx, c2)
    legend off
    xlabel('')
    ylabel('Duration [m]')
    xticklabels(as)
    title('Duration Total')
    box off
    set(gca, 'TickLength', [0 0])
    
    % Duration BS
    subplot(3, 4, 8)
    boxplot(([bsDur{:}] / fs / 60) ./ ([recDur{:}] / fs / 60), gidx,...
        'BoxStyle', 'outline', 'Color', c2, 'notch', 'off')
    hold on
    gscatter(gidx, ([bsDur{:}] / fs / 60) ./ ([recDur{:}] / fs / 60), gidx, c2)
    legend off
    xlabel('')
    ylabel('BS / total')
    xticklabels(as)
    title('Duration of BS')
    box off
    set(gca, 'TickLength', [0 0])
    
    % Duration B
    subplot(3, 4, 12)
    boxplot(([bDur{:}] / fs / 60) ./ ([recDur{:}] / fs / 60), gidx,...
        'BoxStyle', 'outline', 'Color', c2, 'notch', 'off')
    hold on
    gscatter(gidx, ([bDur{:}] / fs / 60) ./ ([recDur{:}] / fs / 60), gidx, c2)
    legend off
    xlabel('')
    ylabel('B / total')
    xticklabels(as)
    title('Duration of B')
    box off
    set(gca, 'TickLength', [0 0])
    
    % duration bar
    subplot(3, 4, 9)
    b = bar(([bsmean; bmean; recmean]' / fs / 60), 'stacked',...
        'FaceColor', 'flat');
    for i = 1 : 4
        b(1).CData(i, :) = c(i, :) + 0.8;
        b(2).CData(i, :) = c(i, :) + 0.4;
        b(3).CData(i, :) = c(i, :);
    end
    axis tight
    xticklabels(as)
    set(gca,'TickLabelInterpreter','none')
    legend({'0.3 < BSR < 0.8', 'BSR < 0.3', 'BSR > 0.8'})
    ylabel('Duration [m]')
    title('Duration')
    box off
    set(gca, 'TickLength', [0 0])
    
    % threshold bar
    subplot(3, 4, 10)
    boxplot([thr{:}], gidx, 'PlotStyle', 'traditional',...
        'BoxStyle', 'outline', 'Color', c2, 'notch', 'off')
    hold on
    gscatter(gidx, [thr{:}], gidx, c2)
    legend off
    xlabel('')
    xticklabels(as)
    ylabel('Threshold [mV]')
    title('Threshold')
    box off
    set(gca, 'TickLength', [0 0])
    
    if saveFig
        figname = ['summary'];
        export_fig(figname, '-tif', '-transparent')
        % savePdf(figname, basepath, ff)
    end
end



