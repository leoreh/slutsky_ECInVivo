
% demonstrate MFR homeostasis by ploting traces from three sessions
close all
b2uv = 0.195;
dur = 0.3;
saveFig = false;

basepath = 'F:\Data\Processed\lh58';
dirnames = ["lh58_200830_090851";
    "lh58_200831_080808";
    "lh58_200903_080936"];

cd(fullfile(basepath, dirnames(1)))
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
tstamps = 1 / fs : 1 / fs : dur;

% dirOrder = randperm(3);
% dirnames = dirnames(dirOrder);
pathPieces = regexp(dirnames(:), '_', 'split'); % assumes filename structure: animal_date_time
sessionDate = [pathPieces{:}];
sessionDate = sessionDate(2 : 3 : end);
titleLabels = ["Baseline"; "12 Hours CNO"; "96 Hours CNO"];

start(1) = 69 * 60 + 58.8;
start(2) = 58 * 60 + 39.7;
start(3) = 151 * 60 + 45.35;

fh = figure;
k = 1;
for i = 1 : length(dirnames)
    filepath = fullfile(basepath, dirnames{i});
    cd(filepath)
    [~, basename] = fileparts(filepath);
    filename = [char(basename) '.dat'];
    
    %     start = randperm(1e8, 1) / fs;
    sig = bz_LoadBinary(filename, 'start', start(i), 'duration', dur,...
        'nChannels', nchans, 'channels', spkgrp{3});
    sig = double(sig) * b2uv;
    sig = sig - mean(sig);
    
    sb{i} = subplot(3, 1, i);
    yOffset = 500;
    for ii = 1 : size(sig, 2)
        plot(tstamps, sig(:, ii) + (yOffset / 2) * ii, 'k', 'LineWidth', 1)
        hold on
    end
    ylim([-100 1400])
    ylabel('Amplitude [uV]')
    yticks([])
    yticklabels('')
    title(titleLabels{i})
    xticks([])
end
xticks([0 0.3])
xlabel('Time [s]')

if saveFig
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = [figpath '\MFRhomeostasis'];
    export_fig(figname, '-tif', '-r300', '-transparent')
end