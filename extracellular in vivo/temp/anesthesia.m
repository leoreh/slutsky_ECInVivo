
basepath = 'E:\Data\Dat\lh39';
cd(basepath)
filename = dir('*.abf');
files = {filename.name};
nfiles = 1 : length(files);     % address specific files

forceLoad = true;
saveFig = true;
tetrodes = true;

for i = 1
    
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
                filename = [basename '.abf'];
                [lfp.data, info] = abf2load(filename);
                
                save([basename '_lfp.mat'], 'lfp')
                save([basename '_info.mat'], 'info')
            end
            fs_orig = info.fADCSequenceInterval;
            fs_orig = 1 / (fs_orig / 1000000);
            fs = 1250;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prepare signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % resmaple
    if ~tetrodes
        sig = resample(double(lfp.data), fs, round(fs_orig));
        sig(end : -1 : end - 60 * fs) = [];
        tstamps = [1 : length(sig)] / fs;
    end
    
    % filter
    linet = lineDetect('x', sig, 'fs', fs, 'graphics', false);
    sig = lineRemove(sig, linet, [], [], 0, 1);
    
    %         x = filterLFP(sig, 'fs', fs, 'stopband', [45 55], 'order', 6,...
    %             'type', 'butter', 'dataOnly', true, 'graphics', false,...
    %             'saveVar', false);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delta power
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % broad-band spectrogram
    freq = logspace(0, 2, 100);
    winsize = 1;       % win length [s]
    win = hann(2 ^ nextpow2(winsize * fs));
    
    [s, f, t, p] = spectrogram(sig, win, round(length(win) / 10), freq, fs,...
        'yaxis', 'psd');
    
    % z-score. great way of comparing changes within a signal
    z = zscore(10 * log10(abs(p)));
    
    % integrate power over delta and sigma band
    deltaf = [1 4];
    [~, deltaidx] = min(abs(f - deltaf));
    zdelta = sum(z(deltaidx(1) : deltaidx(2), :), 1);
    sigmaf = [9 25];
    [~, sigmaidx] = min(abs(f - sigmaf));
    zsigma = sum(z(sigmaidx(1) : sigmaidx(2), :), 1);
    
    smf = round(15 / mode(diff(t)));
    zdelta = bz_NormToRange(zdelta, [0 1]);
    zdelta = smooth(zdelta, smf);
    zsigma = bz_NormToRange(zsigma, [0 1]);
    zsigma = smooth(zsigma, smf);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % burst suppression
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    vars = {'std', 'sum', 'max'};
    bs = getBS('sig', sig, 'fs', fs, 'basepath', basepath,...
        'graphics', true, 'saveVar', false, 'binsize', 1,...
        'clustmet', 'gmm', 'vars', vars, 'basename', basename,...
        'saveFig', false, 'forceA', true);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    iis = getIIS('sig', sig, 'fs', fs, 'basepath', basepath,...
        'graphics', true, 'saveVar', false, 'binsize', 600,...
        'marg', [], 'basename', basename, 'thr', 50,...
        'saveFig', false, 'forceA', true);
    
    
    
    
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
    hold on
    plot(t / 60, zdelta, 'r')
    plot(t / 60, zsigma, 'b')
    legend({'BSR', '[1-4 Hz]', '[9-25 Hz]'})
    xlabel('Time [min]');
    ylabel('[a.u.]')
    axis tight
    set(gca, 'TickLength', [0 0])
    box off
    title('Delta power and BSR')
    ylim([0 1])
    
    
    if saveFig
        figname = [basename '_anesthesia'];
        export_fig(figname, '-tif', '-transparent')
        % savePdf(figname, basepath, ff)
    end
    
    
    
end