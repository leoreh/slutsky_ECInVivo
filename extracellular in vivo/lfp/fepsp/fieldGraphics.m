% field graphics
% created to clean field.m

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


savePdf('io_stimL', basepath, fh)