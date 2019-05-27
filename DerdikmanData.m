% for i = 1 : 100
%     figure
%     plot(library_of_acceptable_spike_shapes(i, :))
% end

basepath = 'E:\Data\Others\DerdikmanSlutsky\Data_with_clusters_after_pre_process';
files = dir(basepath);
dirFlags = [files.isdir];
sub = files(dirFlags);

% load and arrange units
t = 1;
for k = 3 : length(sub)
        foldername = fullfile(basepath, sub(k).name);
        cd(foldername)
        ntet = length(dir('*.Ntt'));
    for i = 1 : ntet
        load_ntt_from_GUI
        units = unique(tmp_CellNumbers);
        for j = 1 : length(units)
            idx = find(tmp_CellNumbers == units(j));
            spikes.rawWaveform{t} = -mean(tmp_Samples(:, :, idx), 3);
            spikes.UID(t) = units(j);
            spikes.tetrode(t) = i;
            spikes.dir{t} = sub(k).name;
            t = t + 1;
        end
    end
end


spikes.samplingRate = 32000;
nunits = length(spikes.UID);

for i = 1 : nunits
    spikes.rawWaveform{i} = spikes.rawWaveform{i}';
    for j = 1 : size(spikes.rawWaveform{i}, 1)
        [~, pos] = max(peak2peak(spikes.rawWaveform{i}, 2));
        spikes.rawWaveform{i} = spikes.rawWaveform{i}(pos, :);
    end
end


for i = 1 : nunits
    f = figure;
    plotWaveform(spikes.rawWaveform{i}', [], [], 'vert', spikes.samplingRate)
    txt = sprintf('cluster %d tetrode %d', spikes.UID(i), spikes.tetrode(i));
    title(txt);
    saveas(f, txt, 'png')
end

cd(basepath)
save('spikes.mat', 'spikes')
