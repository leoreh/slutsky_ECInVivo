
force = false;

if force
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename = 'C:\Users\heiml\AppData\Local\Temp\OneNote\16.0\Exported\{46DB0164-023F-4175-B436-E0C76B451EC2}\NT\2\Cohort2_s13.xlsx';
    ndays = 12;
    
    % initialize
    tmaze.duration = [];
    tmaze.correct = [];
    tmaze.trials = [];
    tmaze.day = [];
    tmaze.delay = [];
    for i = 1 : ndays
        sheet = sprintf('Sheet%d', i);
        [num, ~, raw] = xlsread(filename, sheet);
        
        tmaze.date{i} = raw{1, 1};
        tmaze.mice{i} = raw(3, ~isnan(num(1, :)));
        tmaze.weight(i, :) =  num(2, 3 : 2 : end);
        tmaze.duration = [tmaze.duration; num(4 : end, 3 : 2 : end)];
        tmaze.delay = [tmaze.delay; num(4 : end, 2)];
        tmaze.correct = [tmaze.correct; num(4 : end, 4 : 2 : end)];
        tmaze.success(i, :) = nanmean(num(4 : end, 4 : 2 : end));
        tmaze.trials = [tmaze.trials; [1 : length(num(4 : end, 1))]'];
        tmaze.day = [tmaze.day; ones(length(num(4 : end, 1)), 1) * i];
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
f = figure;

% mean success across sessions
subplot(1, 4, 3)
mice = 5 : 9;
days = 1 : 10;
ndays = length(days);
s = tmaze.success(days, mice);
x = 1 : ndays;
savg = mean(s', 1);
sstd = nanstd(s', 0, 1);
errbounds = [(savg) + (sstd);...
    (savg) - (sstd)];
p = patch([x, x(end : -1 : 1)], [errbounds(1 ,:), errbounds(2, end : -1 : 1)], [.5 .5 .5]);
p.EdgeColor = 'none';
p.FaceAlpha = 0.5;
hold on
plot(savg, 'lineWidth', 3, 'Color', 'k')    
axis tight
x = xlim;
plot(x, [0.5 0.5], '--k')
xticks(1 : ndays)
xticklabels(num2cell([1 : ndays]));
ylim([0.5 1])
yticks([0 0.5 1])
yticklabels(num2cell([0 0.5 1]));
box off
set(gca, 'TickLength', [0 0])
ylabel('Success Rate [%]')
xlabel('Time [Days]')

% perfrmance vs delay
subplot(1, 4, [4])
plot(nanmean(delay.avg', 1), 'LineWidth', 3, 'Color', 'k')
hold on
savg = nanmean(delay.avg', 1);
sstd = std(delay.avg', [], 1, 'omitnan');
errbounds = [(savg) + (sstd);...
    (savg) - (sstd)];
x = 1 : 5;
p = patch([x, x(end : -1 : 1)], [errbounds(1 ,:), errbounds(2, end : -1 : 1)], [.5 .5 .5]);
p.EdgeColor = 'none';
p.FaceAlpha = 0.5;
xticks([1 : length(delay.times)])
xticklabels(num2cell(delay.times))
axis tight
x = xlim;
plot(x, [0.5 0.5], '--k')
ylim([0.5 1])
yticks([0 0.5 1])
yticklabels(num2cell([0 0.5 1]));
box off
set(gca, 'TickLength', [0 0])
xlabel('Delay [s]')

% norm weight
% weight in gram and normalized
[num, txt, raw] = xlsread(filename, 'Weight');
weightMat = num([2, 5 : end], 2 : end);
subplot(1, 4, 2)
mat = weightMat ./ weightMat(1, :);
savg = nanmean(mat', 1);
sstd = std(mat', [], 1, 'omitnan');
errbounds = [(savg) + (sstd);...
    (savg) - (sstd)];
x = 1 : 15;
p = patch([x, x(end : -1 : 1)], [errbounds(1 ,:), errbounds(2, end : -1 : 1)], [.5 .5 .5]);
p.EdgeColor = 'none';
p.FaceAlpha = 0.5;
hold on
plot(savg, 'LineWidth', 3, 'Color', 'k')
box off
set(gca, 'TickLength', [0 0])
ylabel('Norm. Weight')
xlabel('Time [Days]')
ylim([0.7 1])
yticks([0.7 0.8 0.85 1])
yticklabels(num2cell([0.7 0.8 0.85 1]))
xlim([1 ndays])
x = xlim;
plot(x, [0.85 0.85], '--k')
plot(x, [0.8 0.8], '--r', 'LineWidth', 2)

% schematics
subplot(1, 4, 1)
box off
yticks([])
xticks([])


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % learning - trial duration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % trial duration vs. trials. correct blue, wrong red, nan black.
% mice = 1 : 9;
% f = figure;
% clip = [0 90];
% k = 1;
% for j = mice
%     subplot(length(mice), 1, k)
%     x = min(tmaze.duration(:, j), clip(2), 'includenan');
%     plot(x, '.k', 'MarkerSize', 12)
%     hold on
%     idx = find(tmaze.correct(:, j) == 1);
%     plot(idx, x(idx), '.b', 'MarkerSize', 12)
%     idx = find(tmaze.correct(:, j) == 0);
%     plot(idx, x(idx), '.r', 'MarkerSize', 12)
%     idx = find(tmaze.trials == 1);
%     xticks(idx)
%     plot(repmat(idx, 2), clip, '--k')
%     xlim([1 length(x)])
%     ylim(clip)
%     yticks([])
%     box off
%     set(gca, 'TickLength', [0 0])
%     title(tmaze.mice{1}{j})
%     k = k + 1;
%     if j == mice(end)
%         xlabel('Trial [#]')
%         ylabel('Trial Duration [s]')
%         yticks(clip)
%     end
% end
% legend({'NaN', 'Correct', 'Wrong'})
% txt = 'Learning to Run';
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % learning - insight
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mice = 1 : 4;
% days = 7 : 9;
% ndays = length(days);
% f = figure;
% x = tmaze.success(days, mice);
% plot(x, 'LineWidth', 2)
% hold on
% errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% x = xlim;
% plot(x, [0.5 0.5], '--k')
% xticks(1 : ndays)
% xticklabels(num2cell([1 : ndays]));
% ylim([0.5 1])
% yticks([0.5 0.75 1])
% yticklabels(num2cell([0.5 0.75 1]));
% box off
% set(gca, 'TickLength', [0 0])
% ylabel('Success Rate [%]')
% xlabel('Time [Days]')
% legend(tmaze.mice{10}(mice), 'Location', 'South East')
% txt = 'insight';
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % learning - success
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % mean success across sessions
% mice = 5 : 9;
% days = 1 : 10;
% ndays = length(days);
% f = figure;
% subplot(length(mice) * 2, 1, [1 : length(mice) - 1])
% x = tmaze.success(days, mice);
% plot(x, 'LineWidth', 2)
% hold on
% errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% x = xlim;
% plot(x, [0.5 0.5], '--k')
% xticks(1 : ndays)
% xticklabels(num2cell([1 : ndays]));
% ylim([0.5 1])
% yticks([0.5 0.75 1])
% yticklabels(num2cell([0.5 0.75 1]));
% box off
% set(gca, 'TickLength', [0 0])
% ylabel('Success Rate [%]')
% xlabel('Time [Days]')
% legend(tmaze.mice{10}(mice), 'Location', 'South East')
% 
% % binned success across sessions
% binsize = 4;
% bin = cell(length(mice), 1);
% idx = diff(find(tmaze.trials == 1));
% nbins = ceil(idx / binsize);
% binidx = cumsum(nbins);
% % ceil(find(tmaze.trials == 1) / binsize);
% k = 1;
% for j = mice
%     for i = 1 : ndays
%         y = nan(nbins(i), 1)';
%         x = tmaze.correct(tmaze.day == i, j);
%         x = x(~isnan(x));
%         
%         % divide to bins of n values
%         m = numel(x);
%         y(1 : ceil(m / binsize)) = nanmean(reshape([x(:); nan(mod(-m, binsize), 1)], binsize, []));
%         bin{k} = [bin{k}, y];  % calculate the mean over the 1st dim
%     end
%     
%     % plot
%     subplot(length(mice) * 2, 1, length(mice) + k)
%     plot(bin{k})
%     hold on
%     xlim([1 sum(nbins)])
%     xticks([1; binidx + 1])
%     xticklabels(num2cell([1; binidx(2 : end) * binsize]))
%     plot(repmat(binidx + 1, 2), [0 1], '--k')
%     box off
%     set(gca, 'TickLength', [0 0])
%     title(tmaze.mice{1}{j})
%     
%     k = k + 1;
%     
%     if j == mice(end)
%         xlabel('Trial [#]')
%         ylabel('Success Rate [%]')
%     end
% end
% txt = 'Learning to Alternate';
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % weight
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % weight in gram and normalized
% [num, txt, raw] = xlsread(filename, 'Weight');
% weightMat = num([2, 5 : end], 2 : end);
% f = figure;
% subplot(2, 2, 1)
% plot(weightMat)
% xlim([1 ndays])
% box off
% set(gca, 'TickLength', [0 0])
% ylabel('Weight [Kg]')
% 
% subplot(2, 2, 3)
% plot(weightMat ./ weightMat(1, :))
% hold on
% box off
% set(gca, 'TickLength', [0 0])
% ylabel('Norm. Weight')
% xlabel('Time [Days]')
% ylim([0.7 1])
% yticks([0.7 0.8 0.85 1])
% yticklabels(num2cell([0.7 0.8 0.85 1]))
% xlim([1 ndays])
% x = xlim;
% plot(x, [0.85 0.85], '--k')
% plot(x, [0.8 0.8], '--r', 'LineWidth', 2)
% 
% % performance vs. weight
% ssdays = 5 : 10;
% % success rate
% subplot(2, 2, 2)
% plot(weightMat(ssdays + 3, :) ./ weightMat(1, :), tmaze.success(ssdays, :), '.', 'MarkerSize', 8)
% box off
% set(gca, 'TickLength', [0 0])
% xlim([0.7 1])
% ylim([0.5 1])
% yticks([0.5 1])
% xticks([0.7 1])
% ylabel('Success Rate [%]')
% txt = 'Performance vs. Weight';
% % title(txt)
% % running time
% k = 1;
% clear x
% for i = ssdays
%     x(k, :) = nanmean(tmaze.duration(tmaze.day == i, :));
%     k = k + 1;
% end
% subplot(2, 2, 4)
% plot(weightMat(ssdays + 3, :) ./ weightMat(1, :), x, '.', 'MarkerSize', 8)
% box off
% set(gca, 'TickLength', [0 0])
% xlabel('Norm. Weight')
% xlim([0.7 1])
% xticks([0.7 1])
% ylabel('Trial Duration [s]')
% legend(tmaze.mice{10}(mice), 'Location', 'North East')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % changing environment
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ssdays = 5 : 10;
% mice = 5 : 9;
% idx = [find(tmaze.day == ssdays(1), 1) : find(tmaze.day == ssdays(end), 1, 'last')];
% x = zeros(2, length(mice));
% x(1, :) = nanmean(tmaze.correct(idx, mice));
% x(2, :) = nanmean(tmaze.correct(tmaze.day == 11, mice));
% f = figure;
% plot(x, 'LineWidth', 2)
% hold on
% errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% xlim([0.8 2.2])
% x = xlim;
% plot(x, [0.5 0.5], '--k')
% xticks([1 2])
% xticklabels({'familiar', 'novel'})
% ylim([0.25 1])
% yticks([0.25 0.5 0.75 1])
% yticklabels(num2cell([0.25 0.5 0.75 1]));
% box off
% set(gca, 'TickLength', [0 0])
% ylabel('Success Rate [%]')
% xlabel('Environment')
% txt = 'Behavioral Perturbation';
% legend(tmaze.mice{10}(mice), 'Location', 'South East')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % performance vs. delay
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mice = 1 : 4;
% x = tmaze.correct(tmaze.day == 10, mice);
% delay1.times = sort([30 10 50 20]);
% delay1.avg = [nanmean(x(10 : 21, :)); nanmean(x(31 : end, :));...
%     nanmean(x(1 : 9, :)); nanmean(x(22 : 30, :));];
% 
% f = figure;
% subplot(2, 2, 1)
% plot(delay1.avg, 'LineWidth', 2)
% hold on
% errorbar(nanmean(delay1.avg, 2), std(delay1.avg, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% xticks([1 : length(delay1.times)])
% xticklabels(num2cell(delay1.times))
% x = xlim;
% plot(x, [0.5 0.5], '--k')
% ylim([0 1])
% yticks([0 0.5 0.75 1])
% yticklabels(num2cell([0 0.5 0.75 1]));
% box off
% set(gca, 'TickLength', [0 0])
% ylabel('Success Rate [%]')
% xlabel('Delay [s]')
% txt = 'Performance vs. Delay';
% title('Day 1')
% 
% x = tmaze.correct(tmaze.day == 11, mice);
% delay2.times = sort([30 10 50]);
% delay2.avg = zeros(3, 4);
% delay2.avg(1, :) = [nanmean(x(21 : end, 1)); nanmean(x(1 : 9, 2));...
%     nanmean(x(10 : 20, 3)); nanmean(x(21 : end, 4))];
% delay2.avg(2, :) = [nanmean(x(10 : 20, 1)); nanmean(x(10 : 20, 2));...
%     nanmean(x(1 : 9, 3)); nanmean(x(1 : 9, 4))];
% delay2.avg(3, :) = [nanmean(x(1 : 9, 1)); nanmean(x(21 : end, 2));...
%     nanmean(x(21 : end, 3)); nanmean(x(10 : 20, 4))];
% 
% subplot(2, 2, 2)
% plot(delay2.avg, 'LineWidth', 2)
% hold on
% errorbar(nanmean(delay2.avg, 2), std(delay2.avg, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% xticks([1 : length(delay2.times)])
% xticklabels(num2cell(delay2.times))
% x = xlim;
% plot(x, [0.5 0.5], '--k')
% ylim([0 1])
% yticks([0 0.5 0.75 1])
% yticklabels(num2cell([0 0.5 0.75 1]));
% box off
% set(gca, 'TickLength', [0 0])
% xlabel('Delay [s]')
% title('Day 2')
% 
% x = tmaze.correct(tmaze.day == 12, mice);
% delay3.times = sort([10 120]);
% delay3.avg = zeros(2, 4);
% delay3.avg(1, :) = [nanmean(x(11 : end, 1)); nanmean(x(1 : 9, 2));...
%     nanmean(x(11 : end, 3)); nanmean(x(1 : 9, 4))];
% delay3.avg(2, :) = [nanmean(x(1 : 9, 1)); nanmean(x(11 : end, 2));...
%     nanmean(x(1 : 9, 3)); nanmean(x(11 : end, 4))];
% 
% subplot(2, 2, 3)
% plot(delay3.avg, 'LineWidth', 2)
% hold on
% errorbar(nanmean(delay3.avg, 2), std(delay3.avg, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% xticks([1 : length(delay3.times)])
% xticklabels(num2cell(delay3.times))
% x = xlim;
% plot(x, [0.5 0.5], '--k')
% ylim([0 1])
% yticks([0 0.5 0.75 1])
% yticklabels(num2cell([0 0.5 0.75 1]));
% box off
% set(gca, 'TickLength', [0 0])
% xlabel('Delay [s]')
% title('Day 3')
% legend(tmaze.mice{10}(mice), 'Location', 'South East')
% 
% delay.times = sort([10 20 30 50 120]);
% delay.avg = nan(5, 12);
% delay.avg(1, :) = [delay1.avg(1, :), delay2.avg(1, :),...
%     delay3.avg(1, :)];
% delay.avg(2, 1 : 4) = delay1.avg(2, :);
% delay.avg(3, 1 : 8) = [delay1.avg(3, :), delay2.avg(2, :)];
% delay.avg(4, 1 : 8) = [delay1.avg(4, :), delay2.avg(3, :)];
% delay.avg(5, 1 : 4) = delay3.avg(2, :);
% 
% 
% subplot(2, 2, 4)
% plot(delay.avg, '*b', 'LineWidth', 2)
% hold on
% errorbar(nanmean(delay.avg, 2), std(delay.avg, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% xticks([1 : length(delay.times)])
% xticklabels(num2cell(delay.times))
% x = xlim;
% plot(x, [0.5 0.5], '--k')
% ylim([0 1])
% yticks([0 0.5 0.75 1])
% yticklabels(num2cell([0 0.5 0.75 1]));
% box off
% set(gca, 'TickLength', [0 0])
% xlabel('Delay [s]')
% title('Combined')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performance vs. delay % change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mice = 1 : 4;
% x = tmaze.correct(tmaze.day == 10, mice);
% delay1.times = sort([30 50 20]);
% delay1.avg = [nanmean(x(31 : end, :));...
%     nanmean(x(1 : 9, :)); nanmean(x(22 : 30, :))] ./ nanmean(x(10 : 21, :));
% 
% f = figure;
% subplot(2, 2, 1)
% plot(delay1.avg, 'LineWidth', 2)
% hold on
% errorbar(nanmean(delay1.avg, 2), std(delay1.avg, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% xticks([1 : length(delay1.times)])
% xticklabels(num2cell(delay1.times))
% x = xlim;
% ylim([0 1])
% yticks([0 1])
% yticklabels(num2cell([0 1]));
% box off
% set(gca, 'TickLength', [0 0])
% ylabel('Change from 10-s delay [%]')
% xlabel('Delay [s]')
% txt = 'Performance vs. Delay, change';
% title('Day 1')
% legend(tmaze.mice{10}(mice), 'Location', 'South East')
% 
% x = tmaze.correct(tmaze.day == 11, mice);
% delay2.times = sort([30 50]);
% delay2.avg = zeros(3, 4);
% delay2.avg(1, :) = [nanmean(x(21 : end, 1)); nanmean(x(1 : 9, 2));...
%     nanmean(x(10 : 20, 3)); nanmean(x(21 : end, 4))];
% delay2.avg(2, :) = [nanmean(x(10 : 20, 1)); nanmean(x(10 : 20, 2));...
%     nanmean(x(1 : 9, 3)); nanmean(x(1 : 9, 4))];
% delay2.avg(3, :) = [nanmean(x(1 : 9, 1)); nanmean(x(21 : end, 2));...
%     nanmean(x(21 : end, 3)); nanmean(x(10 : 20, 4))];
% delay2.avg = delay2.avg ./ delay2.avg(1, :);
% delay2.avg(1, : ) = [];
% 
% subplot(2, 2, 2)
% plot(delay2.avg, 'LineWidth', 2)
% hold on
% errorbar(nanmean(delay2.avg, 2), std(delay2.avg, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% xticks([1 : length(delay2.times)])
% xticklabels(num2cell(delay2.times))
% x = xlim;
% ylim([0 1])
% yticks([0 1])
% yticklabels(num2cell([0 1]));
% box off
% set(gca, 'TickLength', [0 0])
% xlabel('Delay [s]')
% title('Day 2')
% 
% x = tmaze.correct(tmaze.day == 12, mice);
% delay3.times = sort([120]);
% delay3.avg = zeros(1, 4);
% delay3.avg(1, :) = [nanmean(x(11 : end, 1)); nanmean(x(1 : 9, 2));...
%     nanmean(x(11 : end, 3)); nanmean(x(1 : 9, 4))];
% delay3.avg(2, :) = [nanmean(x(1 : 9, 1)); nanmean(x(11 : end, 2));...
%     nanmean(x(1 : 9, 3)); nanmean(x(11 : end, 4))];
% delay3.avg = delay3.avg ./ delay3.avg(1, :);
% delay3.avg(1, : ) = [];
% 
% subplot(2, 2, 3)
% plot([1 1 1 1], delay3.avg, '*', 'LineWidth', 2)
% hold on
% errorbar(nanmean(delay3.avg), std(delay3.avg, 'omitnan'), 'k', 'LineWidth', 4)
% xticks([1 : length(delay3.times)])
% xticklabels(num2cell(delay3.times))
% x = xlim;
% ylim([0 1])
% yticks([0 1])
% yticklabels(num2cell([0 1]));
% box off
% set(gca, 'TickLength', [0 0])
% xlabel('Delay [s]')
% title('Day 3')
% 
% delay.times = sort([20 30 50 120]);
% delay.avg = nan(5, 12);
% delay.avg(2, 1 : 4) = delay1.avg(1, :);
% delay.avg(3, 1 : 8) = [delay1.avg(2, :), delay2.avg(1, :)];
% delay.avg(4, 1 : 8) = [delay1.avg(3, :), delay2.avg(2, :)];
% delay.avg(5, 1 : 4) = delay3.avg(1, :);
% 
% 
% subplot(2, 2, 4)
% plot(delay.avg, '*b', 'LineWidth', 2)
% hold on
% errorbar(nanmean(delay.avg, 2), std(delay.avg, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% xticks([1 : length(delay.times)])
% xticklabels(num2cell(delay.times))
% x = xlim;
% ylim([0 1])
% yticks([0 1])
% yticklabels(num2cell([0 1]));
% box off
% set(gca, 'TickLength', [0 0])
% xlabel('Delay [s]')
% title('Combined')
% 
% %%% prelimanary results 19-aug-19 %%%
% f = figure;
% % mean success across sessions
% mice = 5 : 9;
% days = 1 : 10;
% ndays = length(days);
% f = figure;
% subplot(1, 3, 1:2)
% x = tmaze.success(days, mice);
% plot(x, 'LineWidth', 2)
% hold on
% errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% xlim([1, ndays])
% x = xlim;
% plot(x, [0.5 0.5], '--k')
% xticks(1 : ndays)
% xticklabels(num2cell([1 : ndays]));
% ylim([0 1])
% yticks([0 0.5 0.75 1])
% yticklabels(num2cell([0 0.5 0.75 1]));
% box off
% set(gca, 'TickLength', [0 0])
% ylabel('Success Rate [%]')
% xlabel('Time [Days]')
% legend(tmaze.mice{10}(mice), 'Location', 'South East')
% 
% subplot(1, 3, 3)
% plot(delay.avg, '*b', 'LineWidth', 2)
% hold on
% errorbar(nanmean(delay.avg, 2), std(delay.avg, [], 2, 'omitnan'), 'k', 'LineWidth', 4)
% xticks([1 : length(delay.times)])
% xticklabels(num2cell(delay.times))
% xlim([1, 5])
% x = xlim;
% plot(x, [0.5 0.5], '--k')
% ylim([0 1])
% yticks([0 0.5 0.75 1])
% yticklabels(num2cell([0 0.5 0.75 1]));
% box off
% set(gca, 'TickLength', [0 0])
% xlabel('Delay [s]')
