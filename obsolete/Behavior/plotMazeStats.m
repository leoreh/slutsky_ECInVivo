function plotMazeStats(correct, dur, weight, idx)

%%% success and duration, and weight across sessions %%%
sticks = [1, 2, 4, 6, 7, 9 : 2 : 24];
sidx = sticks(11);

f = figure;

subplot(2, 1, 1)
plot(correct(idx, :))
hold on
errorbar(nanmean(correct(idx, :), 2), std(correct(idx, :), [], 2, 'omitnan'), 'k', 'LineWidth', 2)
ylim([0 1])
y = ylim;
plot([sidx + 0.5, sidx + 0.5], y, '--k', 'LineWidth', 2)
yticks([0 0.5 0.75 1])
xticks([])
xlim([1 size(correct(idx, :), 1)])
box off
set(gca, 'TickLength', [0 0])
ylabel('Success Rate [%]')

subplot(2, 1, 2)
plot(dur(idx, :))
hold on
errorbar(nanmean(dur(idx, :), 2), std(dur(idx, :), [], 2, 'omitnan'), 'k', 'LineWidth', 2)

% {'30.06', '01.07', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'})
box off
set(gca, 'TickLength', [0 0])
ylabel('Running Duration [s]')
xlabel('Day [#]')

% if ~isempty(weight)
%     subplot(3, 1, 3)
%     plot(weight(idx, :))
%     hold on
%     errorbar(nanmean(weight(idx, :), 2), std(weight(idx, :), [], 2, 'omitnan'), 'k', 'LineWidth', 2)
%     
%     xlim([1 size(weight(idx, :), 1)])
%     box off
%     set(gca, 'TickLength', [0 0])
%     ylabel('Weight [g]')
% end

% xticks([1 : size(correct, 1)])
axis tight
xlim([1 size(dur(idx, :), 1)])
y = ylim;
plot([sidx + 0.5, sidx + 0.5], y, '--k', 'LineWidth', 2)
xticks(sticks)
xticklabels(num2str([1 : 13]'))

end