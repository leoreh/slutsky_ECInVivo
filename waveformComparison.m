% spike width comparison

basepath{1} = 'E:\Data\Styr2019\TERI';
basepath{2} = 'E:\Data\Styr2019\DMSO';
basepath{3} = 'E:\Data\Others\Buzsaki\HC-2\ec013.527';
basepath{4} = 'E:\Data\Others\Buzsaki\Mouse17-130129';
basepath{5} = 'E:\Data\Others\Refaela\Bruce\bruce_120519';
basepath{6} = 'E:\Data\Others\Refaela\Bruce\bruce_140519';
basepath{7} = 'E:\Data\Others\Refaela\Bane\310319HC';
basepath{8} = 'E:\Data\Others\Refaela\alfred_150519';
basepath{9} = 'E:\Data\Others\DerdikmanSlutsky\Data_with_clusters_after_pre_process';

upsamp = 4;


for j = 1 : length(basepath)
    %     spikes = getSpikes('basepath', basepath{j}, 'saveMat', false, 'noPrompts', true);
    
    [~, basename] = fileparts(basepath{j});
    filename = fullfile(basepath{j}, [basename '.spikes.mat']);
    load(filename)
       
    nunits = length(spikes.UID);
    waves = cat(1, spikes.rawWaveform{:})';
    fs = spikes.samplingRate(1);
    clear tp spkw
    
    for i = 1 : nunits
        w = waves(:, i);
        w = fft_upsample(waves(:, i), upsamp);
        [minval, minpos] = min(w);
        [maxval, maxpos] = max(w(1 : minpos - 1));
        [maxvalpost, maxpost] = max(w(minpos : end));
        if ~isempty(maxpost)
            % trough-to-peak - Stark et al., 2013; Bartho et al., 2004
            tp(i) = maxpost;
            if ~isempty(maxval)
                % asymmetry - Sirota et al., 2008
                asym(i) = (maxvalpost - maxval) / (maxvalpost + maxval);
            end
        else
            warning('waveform may be corrupted')
            tp(i) = NaN;
            asym(i) = NaN;
        end
        
        %%% spike width
        [~, halfpos1] = min(abs((minval / 2) - w(1 : length(w) / 2)));
        [~, halfpos2] = min(abs((minval / 2) - w(length(w) / 2 + 1 : end)));
        spkw(i) = (halfpos2 - halfpos1) / fs / upsamp * 1000;

    end
    % samples to ms
    tp = tp / fs * 1000 / upsamp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spike width by inverse of max frequency in spectrum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for i = 1 : nunits
%         w = waves(:, i);
%         
%         % w = fft_upsample(waves(:, i), upsamp);
%         w = [w(1) * ones(1000, 1); w; w(end) * ones(1000, 1)];
%         [wave f t] = getWavelet(w, fs, 500, 3000, 128);
%         % wt = cwt(w, fs, 'amor', 'FrequencyLimits', [500 3000]);
%         wave = wave(:, int16(length(t) / 4) : 3 * int16(length(t) / 4));
%         
%         % find maximum
%         [maxPow, ix] = max(wave);
%         [~, mix] = max(maxPow);
%         ix = ix(mix);
%         spkw(i) = 1000 / f(ix);
%     end
    
    spikewidth{j} = spkw;
    troughpeak{j} = tp;
    
end

% graphics
names = {'TERI', 'DMSO', 'Buzsaki', 'Buzsaki', 'bruce tdt', 'bruce intan', 'bane', 'alfred', 'Dori'};
close all
figure
for j = 1 : length(basepath)
    subplot(5, 2, j)
    plot(troughpeak{j}, spikewidth{j}, '.')
    hold on
%     ylim([0 1])
    if j == 9
%         ylim([0 2.5])
    end
%     xlim([0 1])
    title(names{j})
    xlabel('trough-to-peak [ms]')
    ylabel('spike width [ms]')
end

