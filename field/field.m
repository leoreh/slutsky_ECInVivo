
% this is a wrapper for fEPSP analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basepath = 'E:\Data\Field\lh47_200212';
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
files = [20];
=======
basepath = 'E:\Data\19.02.20';
files = [25 26];
>>>>>>> parent of cd38659... Update field.m
=======
files = [20];
>>>>>>> parent of 6b59bba... field
=======
files = [20];
>>>>>>> parent of 6b59bba... field

cd(basepath)
filename = dir('*.wcp');
filename = natsort({filename.name});
ch = [1, 3];

% typical params for stp
start = [0.01 : 0.02 : 0.09];
nstim = 5;
% typical params for IO / stability
start = 0.03;
nstim = 1;

% go over files, calc response and save data as .mat
for i = files
    raw = import_wcp(fullfile(basepath, filename{i}));
    fprintf('\nloading %s\n', fullfile(basepath, filename{i}));
    [~, basename] = fileparts(filename{i});
    for j = ch
        [data.amp{j}, data.rm{j}] = getFieldAmp('sig', raw.S{j},...
            'fs', raw.fs, 'start', start, 'stop', [], 'inspect', false,...
            'basepath', basepath, 'nstim', nstim,...
            'graphics', false, 'saveVar', false, 'saveFig', true,...
            'filename', basename);
        data.sig{j} = raw.S{j};
    end   
    % arrange and save
    data.t = raw.T;
    data.fs = raw.fs;
    save(basename, 'data')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% concatenate traces from different files
amp = cell(ch);
k = 1;
for i = files
    [~, basename] = fileparts(filename{i});
    load(basename)
    for j = ch
        amp{j} = [amp{j} data.amp{j}];
        trace{k, j} = data.sig{j};
    end
    ntraces(k) = size(amp{j}, 2);
    k = k + 1;
end

% graphics
tstability = [1 : length(amp{1})] * 30 / 60;
lbs = {'STP', 'PBS', 'STP', 'PSEM', 'STP'};

fh = figure;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
subplot(4, 3, 1 : 3)
scatter(tstability, amp{1}, 20, 'k', 'filled')
hold on
ylabel('Amplitude [mV]')
set(gca, 'TickLength', [0 0])
box off
axis tight
title('Right CA1')
ylim([-0.5 1.5])
addLns('lns', tstability(ntraces), 'lbs', lbs)

subplot(4, 3, 7 : 9)
scatter(tstability, amp{3}, 20, 'k', 'filled')
hold on
xlabel('Time [m]')
ylabel('Amplitude [mV]')
set(gca, 'TickLength', [0 0])
box off
axis tight
title('Left CA1')
ylim([-0.5 1.5])
addLns('lns', tstability(ntraces), 'lbs', lbs)

k = 0;
for i = [1, 3, 5]
    subplot(4, 3, 4 + k)
    stdshade(trace{i, 1}', 0.5, 'k', data.t)
    ylabel('Voltage [mV]')
    if i ~= 1
        set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
            'XColor', 'none', 'YColor', 'none', 'Color', 'none')
    else
        xlabel('Time [s]')
        xticks([0 0.2])
        set(gca, 'TickLength', [0 0])
    end
    box off
    axis tight
    k = k + 1;
end
ylim([-2 1.5])

k = 0;
for i = [1, 3, 5]
    subplot(4, 3, 10 + k)
    stdshade(trace{i, 3}', 0.5, 'k', data.t)
    ylabel('Voltage [mV]')
    if i ~= 1
        set(gca, 'TickLength', [0 0], 'YTickLabel', [], 'XTickLabel', [],...
            'XColor', 'none', 'YColor', 'none', 'Color', 'none')
    else
        xlabel('Time [s]')
        xticks([0 0.2])
        set(gca, 'TickLength', [0 0])
    end
    box off
    axis tight
    k = k + 1;
end
ylim([-2 1.5])

savePdf('stability', basepath, fh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

for j = ch
    subplot(1, 4, j)
    plot(data.t, mean(data.sig{j}'), 'k', 'LineWidth', 2)
    xlabel('Time [s]')
    set(gca, 'TickLength', [0 0])
    box off
    axis tight
    ylim([-1.5 1.5])
    xticks([0 0.2])
    if j == 1
        title('Right CA1')
        ylabel('Amplitude [mV]')
    else
        title('Left CA1')
    end
    
    subplot(1, 4, j + 1)
    errorbar([1 : nstim], mean(data.amp{j}'), std(data.amp{j}'),...
        '--*k', 'LineStyle', 'none', 'LineWidth', 2)
    set(gca, 'TickLength', [0 0])
    box off
    axis tight
    xlabel('Stimulus [#]')
    ylim([0 1])
    xlim([0 6])
    xticks(1 : 5)
end

savePdf('io_stimL', basepath, fh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vector of intensities (mA). The order must correspond to the order of
% files.
intensity = [0.2 : 0.1 : 0.4];
files = [2 : 4];

% arrange data
k = 1;
for i = files
    [~, basename] = fileparts(filename{i});
    load(basename)
    for j = ch
        amp{k, j} = (data.amp{j});
        trace{j}(k, :) = mean(data.sig{j}');
    end
    k = k + 1;
end

% graphics

fh = figure;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

for j = ch
    subplot(1, 4, j)
    errorbar(intensity, ampavg{j}, ampstd{j}, '--ok', 'LineStyle', 'none', 'LineWidth', 3)
    set(gca, 'TickLength', [0 0])
    box off
    axis tight
    if j == 1
        title('Right CA1')
        ylabel('Amplitude [mV]')
    else
        title('Left CA1')
    end
    xlabel('Intensity [mA]')
    ylim([-1 5])
    xlim([0.1 0.5])
    xticks(intensity)
    
    subplot(1, 4, j + 1)
    hold on
    for i = 1 : size(trace, 1)
        plot(data.t, mean(trace{i, j}'))
    end
    if j == 3
        legend(strsplit(num2str(intensity)))
    end
    xlabel('Time [s]')
    set(gca, 'TickLength', [0 0])
    box off
    axis tight
    ylim([-4 4])
    xticks([0 0.2])
end

savePdf('io_stimL', basepath, fh)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intensity = [0.3];
clear filename
filename{1} = 'STP2_5AP_50Hz_0.5msduration_0.4mA';
filename{2} = 'STP2_5AP_50Hz_0.5msduration_0.8mA';



pk_stp = stp('filename', {filename{1}}, 'intensity', intensity,...
    'inspect', true, 'basepath', basepath, 'graphics', graphics,...
    'nstim', 5, 'saveFig', saveFig);
