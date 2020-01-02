
basepath = 'E:\Data\Others\DZ\Field\Acute recordings\2h-3h\WT';
cd(basepath)
filename = dir('*.abf');
files = {filename.name};
nfiles = 1 : length(files);     % address specific files

forceLoad = true;
analyze = true;
saveFig = true;
tetrodes = false;

for i = nfiles
    
    if forceLoad
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % tetrodes
        if tetrodes
            ch = 5;
            [~, basename] = fileparts(basepath);
            load([basename '.lfp.mat'])
            fs = lfp.fs;
            sig = double(lfp.data(:, ch));
            tstamps = lfp.timestamps;
            
            % field
        else
            [~, basename] = fileparts(files{i});
            if exist([basename '_lfp.mat'])
                load([basename '_lfp.mat'])
                load([basename '_info.mat'])
            else
                filename = dir('*abf');
                filename = filename.name;
                [lfp.data, info] = abf2load(filename);
                
                save([basename '_lfp.mat'], 'data')
                save([basename '_info.mat'], 'info')
            end
            fs_orig = info.fADCSequenceInterval;
            fs_orig = 1 / (fs_orig / 1000000);
            fs = 1250;
        end
    end
    
    
    if analyze
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % prepare signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % resmaple
        if ~tetrodes
            sig = resample(double(lfp.data), fs, round(fs_orig));
            sig(end : -1 : end - 60 * fs) = [];
            tstamps = [1 : length(sig)] / fs;
        end

        
        % remove line
        
        % linet = lineDetect('x', x, 'fs', fs, 'graphics', false);
        % [x] = lineRemove(x, linet, [], [], 0, 1);
%         x = filterLFP(sig, 'fs', fs, 'stopband', [45 55], 'order', 6,...
%             'type', 'butter', 'dataOnly', true, 'graphics', false,...
%             'saveVar', false);
%         
%            x = filterLFP(sig, 'fs', fs, 'passband', [0 10], 'order', 6,...
%             'type', 'butter', 'dataOnly', true, 'graphics', true,...
%             'saveVar', false);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % delta power
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % broad-band spectrogram
        freq = logspace(0, 2, 100);
        winsize = 10;       % win length [s]
        win = hann(2 ^ nextpow2(winsize * fs));
        
        [s, f, t, p] = spectrogram(sig, win, round(length(win) / 10), freq, fs,...
            'yaxis', 'psd');
        
        % z-score. great way of comparing changes within a signal
        z = zscore(abs(log10(s)));
        smf = 15;
        specdt = mode(diff(t));
        % sw = smooth(z, smf ./ specdt);
        % z = bz_NormToRange(z, [0 1]);
        
        % integrate power over the delta band
        deltaf = [0.5 4];
        [~, deltaidx] = min(abs(f - deltaf));
        zdelta = mean(z(deltaidx(1) : deltaidx(2), :), 1);
        zdelta = bz_NormToRange(zdelta, [0 1]);
        zdelta = smooth(zdelta, smf ./ specdt);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % burst suppression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        vars = {'std', 'sum', 'max'};
        bs = getBS('sig', sig, 'fs', fs, 'basepath', basepath,...
            'graphics', true, 'saveVar', false, 'binsize', 1,...
            'clustmet', 'gmm', 'vars', vars, 'basename', basename,...
            'saveFig', saveFig, 'forceAnalyze', true);

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % graphics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ff = gcf;
    
    % spectrogram
    s1 = subplot(3, 4, [9 : 10]);
    surf(t / 60, f, 10*log10(abs(p)), 'EdgeColor', 'none');
    axis xy;
    axis tight;
    view(0,90);
    origSize = get(gca, 'Position');
    colormap(jet);
    colorbar;
    ylabel('Frequency [Hz]');
    set(gca, 'YScale', 'log')
    set(s1, 'Position', origSize);
    set(gca, 'TickLength', [0 0])
    box off
    title('Wideband spectrogram')
    
    % delta power
    splot = subplot(3, 4, [5 : 6]);
    subplot(splot)
    yyaxis right
    plot(t / 60, zdelta)
    xlabel('Time [min]');
    ylabel('Norm. z-score')
    axis tight
    set(gca, 'TickLength', [0 0])
    box off
    title('Delta power and BSR')
    ylim([0 1])

    
%     if saveFig
%         figname = [basename '_anesthesia'];
%         export_fig(figname, '-tif', '-transparent')
%         % avePdf(figname, basepath, ff)
%     end
%     
    
    
end