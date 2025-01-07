
basepath = 'E:\Data\Others\DZ\Field\Acute recordings\2h-3h\WT';
cd(basepath)
filename = dir('*.abf');
files = {filename.name};
nfiles = 1 : length(files);     % address specific files

forceA = true;
forceload = false;
saveFig = false;
tetrodes = false;
popGraphics = false;
mouseGraphics = false;
saveVar = false;

maxidx = 10;         % longest recording
% pop.bsrmat = nan(length(files), length(bs.bsr));
% pop.dmat = nan(length(files), length(dband));
% pop.iismat = nan(length(files), length(iis.rate));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = nfiles
    
    if tetrodes
        ch = 5;
        [~, basename] = fileparts(basepath);
    else
        ch = 1;
        [~, basename] = fileparts(files{i});
    end
    
    lfp = getLFP('basepath', basepath, 'chans', ch, 'chavg', {},...
        'fs', 1250, 'interval', [0 inf], 'extension', 'abf', 'pli', true,...
        'savevar', true, 'force', forceload, 'basename', basename);
    
    sig = double(lfp.data(:, ch));
    fs = lfp.fs;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % burst suppression
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vars = {'std', 'max', 'sum'};
    bs = getBS('sig', sig, 'fs', fs, 'basepath', basepath,...
        'graphics', true, 'saveVar', true, 'binsize', 0.5,...
        'clustmet', 'gmm', 'vars', vars, 'basename', basename,...
        'saveFig', false, 'forceA', forceA);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % thr during suppression
    thr(1) = 0;
    thr(2) = mean(sig(~bs.binary)) + 15 * std(sig(~bs.binary));
    
    iis = getIIS('sig', sig, 'fs', fs, 'basepath', basepath,...
        'graphics', true, 'saveVar', saveVar, 'binsize', 300,...
        'marg', [], 'basename', basename, 'thr', thr,...
        'saveFig', false, 'forceA', forceA);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % burst suppression after IIS filteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(iis, 'filtered') && size(iis.wv, 1) > 50
        bs = getBS('sig', iis.filtered, 'fs', fs, 'basepath', basepath,...
            'graphics', false, 'saveVar', true, 'binsize', 1,...
            'clustmet', 'gmm', 'vars', vars, 'basename', basename,...
            'saveFig', false, 'forceA', forceA);    
    end  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % anesthesia state bouts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calc dband
    [dband, tband] = specBand('sig', sig, 'graphics', false, 'band', [1 4]);
    
    % find bouts
    limbs = [0 0.5];
    limd = [0.5 1];
    vec = bs.bsr' > limbs(1) & bs.bsr' < limbs(2) &...
        dband' > limd(1) & dband' < limd(2);
    
    % binsize of bsr and delta power is 10 s. minimum duration is 1 min [6
    % bins] and minimum time between bouts of deep anesthesia is 1 minutes
    % [6 bins]
    bouts = binary2bouts('vec', vec, 'minDur', 6, 'interDur', 6);
    
    % bout bins to sample idx
    bouts = bs.cents(bouts);
    
    % find IIS within anesthesia bouts
    idx = [];
    for j = 1 : size(bouts, 1)
        idx = [idx; find(iis.peakPos > bouts(j, 1) &...
            iis.peakPos < bouts(j, 2))];
    end
    wv = iis.wv(idx, :);
    marg = round(0.1 * fs);
    wvstamps = -marg / fs * 1000 : 1 / fs * 1000 : marg / fs * 1000;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % arrange population data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if exist('pop.mat')
        load('pop.mat')
    else
        pop.recDur(i) = length(sig);
        pop.iis{i} = iis.rate;
        
        % arrange states and IIS within states
        if ~isempty(bouts)
            pop.stateDur(i) = sum(bouts(:, 2) - bouts(:, 1)) / fs / 60;
        end
        pop.IISinState = size(wv, 1);
        
        % equalize sampling rate of bsr and iis
        btimes = (find(~bs.binary));
        bsrbinsize = 300 * fs;      % fit spectrogram w/ 10 s window and 0% overlap
        [pop.bsr{i}, ~, ~] = calcFR(btimes, 'winCalc', [1, length(bs.binary)],...
            'binsize', bsrbinsize, 'smet', 'none');
        
        % equalize sampling rate of delta and iis
        pop.dband{i} = interp1(tband, dband, 1 : 300 : pop.recDur(i) / fs);
        pop.dband{i}(1) = [];
               
        % arrange in matrices after last recording is loaded and save
        if i == length(files)
            
            %         xx = cellfun(@(x)[x(1:end), NaN(1, 90-length(x))], pop.dband,...
            %             'UniformOutput', false);
            %
            %         xx = reshape(cell2mat(xx), 8, 90);
            
            pop.bsrmat(i, 1 : length(bs.bsr)) = bs.bsr;
            pop.dmat(i, 1 : length(dband)) = dband;
            pop.iismat(i, 1 : length(iis.rate)) = iis.rate;
            if i == maxidx
                pop.td = tband;
                pop.tb = bs.cents;
                pop.ti = iis.cents;
            end
            
            if saveVar
                save('pop', 'pop')
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % graphics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mouseGraphics
        
        f = figure;
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        suptitle(basename)
        
        % raw
        sb1 = subplot(3, 2, 1 : 2);
        plot((1 : length(sig)) / fs / 60, sig, 'k')
        axis tight
        ylim([-1 2])
        hold on
        plot(xlim, [iis.thr(2) iis.thr(2)], '--r')
        ylabel('Voltage [mV]')
        yyaxis right
        p1 = plot(iis.cents / fs / 60, iis.rate, 'r', 'LineWidth', 3);
        if ~isempty(p1)
            p1.Color(4) = 0.3;
        end
        ylim([0 1])
        ylabel('Rate [IIS / min]')
        legend({'Raw', 'IIS thr', 'IIS rate'}, 'Location', 'northwest')
        axis tight
        set(gca, 'TickLength', [0 0])
        box off
        
        % bsr and delta
        sb2 = subplot(3, 2, 3 : 4);
        plot(bs.cents / fs / 60, bs.bsr, 'k', 'LineWidth', 2)
        hold on
        plot(tband / 60, dband, 'b', 'LineWidth', 2)
        legend({'BSR', 'Delta'})
        Y = ylim;
        if ~isempty(bouts)
            fill([bouts fliplr(bouts)]' / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
                'b', 'FaceAlpha', 0.2,  'EdgeAlpha', 0, 'HandleVisibility', 'off');
        end
        xlabel('Time [m]')
        ylabel('[a.u.]')
        axis tight
        set(gca, 'TickLength', [0 0])
        box off
        title('Surgical Anesthesia')
        
        % bs
        subplot(3, 2, 5);
        midsig = round(length(sig) / 2 / fs / 60);
        ylimsig = sig(midsig * fs * 60 : (midsig + 2) * fs * 60);
        plot((1 : length(sig)) / fs / 60, sig)
        hold on
        axis tight
        ylim([min(ylimsig) max(ylimsig)])
        Y = ylim;
        fill([bs.stamps fliplr(bs.stamps)] / fs / 60, [Y(1) Y(1) Y(2) Y(2)],...
            'k', 'FaceAlpha', 0.4,  'EdgeAlpha', 0);
        xlim([midsig midsig + 2])
        ylabel('Voltage [mV]')
        set(gca, 'TickLength', [0 0])
        box off
        title('BS')
        
        % iis waveforms
        if ~isempty(iis.wv)
            subplot(3, 2, 6)
            plot(wvstamps, iis.wv)
            ylabel('Voltage [mV]')
            xlabel('Time [ms]')
            axis tight
            xticks([-marg, 0, marg] / fs * 1000);
            set(gca, 'TickLength', [0 0])
            box off
            title('IIS waveform')
            
            
            % mean + std waveform
            axes('Position',[.571 .11 .15 .1])
            box on
            stdshade(iis.wv, 0.5, 'k', wvstamps)
            axis tight
            set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
                'XColor', 'none', 'YColor', 'none', 'Color', 'none')
            title(sprintf('n = %d', size(iis.wv, 1)));
            box off
        end
        
%         linkaxes([sb1, sb2], 'x');
        
        if saveFig
            figname = [basename];
            export_fig(figname, '-tif', '-transparent')
            % savePdf(figname, basepath, ff)
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % population summary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if popGraphics
        
        % graphics
        figure;
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        %     iis
        subplot(3, 3, [1 : 2])
        plot(pop.ti / fs / 60, pop.iismat, 'LineWidth', 1)
        hold on
        stdshade(pop.iismat, 0.5, 'k', pop.ti / fs / 60)
        axis tight
        ylim([0 10])
        ylabel('Rate [IIS / min]')
        set(gca, 'TickLength', [0 0])
        box off
        title('IIS rate [5 min win]')
        
        %     bs
        subplot(3, 3, [4 : 5])
        plot(pop.tb / fs / 60, pop.bsrmat, 'LineWidth', 1)
        hold on
        stdshade(pop.bsrmat, 0.5, 'k', pop.tb / fs / 60)
        axis tight
        ylim([0 1])
        ylabel('BSR')
        set(gca, 'TickLength', [0 0])
        box off
        title('BSR [10 s win]')
        
        %     delta
        subplot(3, 3, [7 : 8])
        plot(pop.td / 60, pop.dmat, 'LineWidth', 1);
        hold on
        stdshade(pop.dmat, 0.5, 'k', pop.td / 60)
        axis tight
        ylim([0 1])
        ylabel('norm. delta power')
        set(gca, 'TickLength', [0 0])
        box off
        title('Delta [10 s win]')
        
        %     iis vs bsr
        subplot(3, 3, 3)
        for i = 1 : length(files)
            if ~isempty(pop.iis{i})
                scatter(pop.bsr{i}, pop.iis{i})
                hold on
            end
        end
        xlabel('BSR')
        ylabel('IIS rate')
        set(gca, 'TickLength', [0 0])
        box off
        title('IIS vs. BSR')
        
        %     iis vs delta
        subplot(3, 3, 6)
        for i = 1 : length(files)
            if ~isempty(pop.iis{i})
                scatter(pop.dband{i}, pop.iis{i})
                hold on
            end
        end
        xlabel('delta power')
        ylabel('IIS rate')
        set(gca, 'TickLength', [0 0])
        box off
        title('IIS vs. Delta')
        
        %     bsr vs delta
        subplot(3, 3, 9)
        for i = 1 : length(files)
            if ~isempty(pop.iis{i})
                scatter(pop.dband{i}, pop.bsr{i})
                hold on
            end
        end
        xlabel('delta power')
        ylabel('BSR')
        set(gca, 'TickLength', [0 0])
        box off
        title('BSR vs. Delta')
        
        if saveFig
            figname = ['param_summary'];
            export_fig(figname, '-tif', '-transparent')
        end
        break
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% population stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% basepath{1} = 'E:\Data\Others\DZ\Field\Acute recordings\Long recordings\WT';
% basepath{2} = 'E:\Data\Others\DZ\Field\Acute recordings\2h-3h\WT';
% basepath{3} = 'E:\Data\Others\DZ\Field\Acute recordings\Long recordings\APPPS1';
% basepath{4} = 'E:\Data\Others\DZ\Field\Acute recordings\2h-3h\APPPS1';
%
% for i = 1 : length(basepath)
%     load([basepath{i} '\anestats.mat']);
%     stats{i} = anestats;
%
%     x(i, :) = nanmean(anestats);
% end
%
%
% figure
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% subplot(1, 2, 1)
% bar(([x(:, 4), x(:, 5) - x(:, 4)]), 'stacked')
% hold on
% for i = 1 : length(basepath)
% scatter(ones(size(stats{i}, 1), 1) * i, (stats{i}(:, 5)), 500, '.k')
% end
% axis tight
% xticklabels({'WT long', 'WT short', 'APP long', 'APP short'})
% legend({'IIS in state', 'IIS out state', 'IIS total'})
% ylabel('Number of IIS')
% title('IIS')
% box off
% set(gca, 'YScale', 'log')
% set(gca, 'TickLength', [0 0])
%
% subplot(1, 2, 2)
% bar(([x(:, 1), x(:, 2) - x(:, 1)]), 'stacked')
% hold on
% for i = 1 : length(basepath)
% scatter(ones(size(stats{i}, 1), 1) * i, (stats{i}(:, 2)), 500, '.k')
% end
% axis tight
% xticklabels({'WT long', 'WT short', 'APP long', 'APP short'})
% legend({'state duration', 'non-state duration', 'recording duration'})
% ylabel('Duration [m]')
% title('State duration')
% box off
% set(gca, 'TickLength', [0 0])
%
% if saveFig
%     figname = ['summary'];
%     export_fig(figname, '-tif', '-transparent')
%     % savePdf(figname, basepath, ff)
% end

