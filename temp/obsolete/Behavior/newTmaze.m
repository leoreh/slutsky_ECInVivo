
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = 'E:\Leore';
filename = 'b3.xlsx';
filename = fullfile(basepath, filename);

nmice = 6;
mnames = [60 61 62 63 66 67];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and arrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : nmice
    
%     t = readtable(filename, 'Sheet', 2);
    
%     t.Properties
    
% m{i}.date = t{1, 2 : 2 : end}
% m{i}.time = t{2, 2 : 2 : end};
% m{i}.weight = t{3, 2 : 2 : end};
% m{i}.dur = t{4 : end, 2 : 2 : end};

    micenames{i} = ['m', num2str(mnames(i))];
    [num, ~, raw] = xlsread(filename, micenames{i});
    m{i}.date = num(1, 2 : 2 : end);
    m{i}.time = num(2, 2 : 2 : end);
    m{i}.weight = num(3, 2 : 2 : end);
    m{i}.delay = num(4, 2 : 2 : end);
    m{i}.dur = num(5 : end, 2 : 2 : end);
    m{i}.dur = m{i}.dur(~isnan(m{i}.dur(:, 1)), :);
    
    % arrange direction
    [maxtrials, nsessions] = size(m{i}.dur);
    m{i}.dir = char(zeros(size(m{i}.dur)));
    k = 3;
    for j = 1 : nsessions
        m{i}.dir(:, j) = char(raw{6 : 5 + maxtrials, k});
        k = k + 2;
    end
    
    % arrange correct trials
    for j = 1 : nsessions
        x = strtrim(m{i}.dir(:, j))';
        for k = 1 : maxtrials - 1
            if x(k) == 'R' && x(k + 1) == 'L' || ...
                    x(k) == 'L' && x(k + 1) == 'R'
                m{i}.correct(k, j) = 1;
            elseif x(k) == 'R' && x(k + 1) == 'R' || ...
                    x(k) == 'L' && x(k + 1) == 'L'
                m{i}.correct(k, j) = 0;
            else
                m{i}.correct(k, j) = nan;
%                 warning('wrong trial direction');
            end
        end
    end   
end

[num, ~, raw] = xlsread(filename, 'Weight');
BLweight = num(1, 2 : end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mice = 1 : nmice; 
night = m{1}.time(1:21) > 1600;
c = 'rbgmky';

%%% combine mice to matrix
for i = mice
    correct(:, i) = nanmean(m{i}.correct, 1);
    dur(:, i) = nanmean(m{i}.dur, 1);
%     weight(:, i) = m{i}.weight / BLweight(i);
end

%%% success and duration across sessions %%%
idx = 1 : nsessions;
plotMazeStats(correct, dur, [], idx)
% title('All')
plotMazeStats(correct, dur, [], night)
title('Night')
plotMazeStats(correct, dur, [], ~night)
title('Day')

%%% night vs. day
figure
subplot(3, 1, 1)
x(1, :) = nanmean(correct(~night, :));    % mean of last three days
x(2, :) = nanmean(correct(night, :));   % mean of last three nights
plot(x)
hold on
errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 2)
box off
set(gca, 'TickLength', [0 0])
xlim([0.8 2.2])
ylabel('Success Rate [%]')
xticks([])
[~, p, ~, stats] = ttest(x(1, :), x(2, :));
txt = sprintf('t = %.2g; p = %.2g', stats.tstat, p);
title(txt)

subplot(3, 1, 2)
x(1, :) = nanmean(dur(~night, :));    % mean of last three days
x(2, :) = nanmean(dur(night, :));   % mean of last three nights
plot(x)
hold on
errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 2)
box off
set(gca, 'TickLength', [0 0])
xlim([0.8 2.2])
ylabel('Running Duration [s]')
xticks([])
[~, p, ~, stats] = ttest(x(1, :), x(2, :));
txt = sprintf('t = %.2g; p = %.2g', stats.tstat, p);
title(txt)

subplot(3, 1, 3)
x(1, :) = nanmean(weight(~night, :));    % mean of last three days
x(2, :) = nanmean(weight(night, :));   % mean of last three nights
plot(x)
hold on
errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 2)
box off
set(gca, 'TickLength', [0 0])
xlim([0.8 2.2])
ylabel('Norm. Weight')
xticks([1 2])
xticklabels({'Day', 'Night'})
[~, p, ~, stats] = ttest(x(1, :), x(2, :));
txt = sprintf('t = %.2g; p = %.2g', stats.tstat, p);
title(txt)
sgtitle('Day vs. Night')

%%% correlation running duration and weight %%%
figure
for i = mice
subplot(3, 2, i)
scatter(weight(:, i), dur(:, i), 'filled', c(i))
box off
set(gca, 'TickLength', [0 0])
ylim([0 120])
xlim([0.75 1])
[rho, pval] = corr([weight(:, i), dur(:, i)], 'Rows', 'complete');
txt = sprintf('m%d; r = %.2g; p = %.2g', i, rho(2, 1), pval(2, 1));
title(txt)
end
subplot(3, 2, 6)
for i = mice
    scatter(weight(:, i), dur(:, i), 'filled', c(i))
    hold on
end
box off
set(gca, 'TickLength', [0 0])
% ylim([0 120])
% xlim([0.75 1])
xlabel('Norm. Weight')
ylabel('Running Duration [s]')
[rho, pval] = corr([weight(:), dur(:)], 'Rows', 'complete');
txt = sprintf('all; r = %.2g; p = %.2g', rho(2, 1), pval(2, 1));
title(txt)

%%% correlation correct and weight %%%
figure
for i = mice
subplot(3, 2, i)
scatter(weight(:, i), correct(:, i), 'filled', c(i))
box off
set(gca, 'TickLength', [0 0])
% ylim([0 120])
% xlim([0.75 1])
[rho, pval] = corr([weight(:, i), correct(:, i)], 'Rows', 'complete');
txt = sprintf('m%d; r = %.2g; p = %.2g', i, rho(2, 1), pval(2, 1));
title(txt)
end
subplot(3, 2, 6)
for i = mice
    scatter(weight(:, i), correct(:, i), 'filled', c(i))
    hold on
end
box off
set(gca, 'TickLength', [0 0])
% ylim([0 120])
% xlim([0.75 1])
xlabel('Norm. Weight')
ylabel('Correct [%]')
[rho, pval] = corr([weight(:), correct(:)], 'Rows', 'complete');
txt = sprintf('all; r = %.2g; p = %.2g', rho(2, 1), pval(2, 1));
title(txt)

%%% weight across sessions %%%
clear weight
clear weightRaw
for i = mice
    weightRaw(:, i) = m{i}.weight;
    weight(:, i) = m{i}.weight / BLweight(i);
end
weightRaw = [BLweight; weightRaw];
weight = [1 1 1 1 1; weight];

sticks = [1, 2, 4, 6, 7, 9 : 2 : 24];
sidx = sticks(11);

f = figure;

subplot(2, 1, 1)
plot(weightRaw(idx, :))
hold on
errorbar(nanmean(weightRaw(idx, :), 2), std(weightRaw(idx, :), [], 2, 'omitnan'), 'k', 'LineWidth', 2)
% ylim([0.75 1])
xlim([1 size(weightRaw(idx, :), 1)])
ylim([20 32])
y = ylim;
x = xlim;
plot([sidx + 0.5, sidx + 0.5], y, '--k', 'LineWidth', 2)
yticks([20, 32])
xticks([])
box off
set(gca, 'TickLength', [0 0])
ylabel('Weight [g]')

subplot(2, 1, 2)
plot(weight(idx, :))
hold on
errorbar(nanmean(weight(idx, :), 2), std(weight(idx, :), [], 2, 'omitnan'), 'k', 'LineWidth', 2)
box off
set(gca, 'TickLength', [0 0])
xlabel('Day [#]')
% xticks([1 : size(correct, 1)])
axis tight
xlim([1 size(dur(idx, :), 1)])
y = ylim;
x = xlim;
plot([sidx + 0.5, sidx + 0.5], y, '--k', 'LineWidth', 2)
plot(x, [0.8 0.8], '--r')
plot(x, [0.85 0.85], '--b')
yticks([0.8 0.85 1])
xticks(sticks)
xticklabels(num2str([1 : 13]'))

%%% pre-post osmotic pump
clear x
post = 22 : 24;
pre = 19 : 21;

figure
subplot(3, 1, 1)
x(1, :) = nanmean(correct(pre, :));    
x(2, :) = nanmean(correct(post, :));   
plot(x)
hold on
errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 2)
box off
set(gca, 'TickLength', [0 0])
xlim([0.8 2.2])
ylabel('Success Rate [%]')
xticks([])
[~, p, ~, stats] = ttest(x(1, :), x(2, :));
txt = sprintf('t = %.2g; p = %.2g', stats.tstat, p);
title(txt)

subplot(3, 1, 2)
x(1, :) = nanmean(dur(pre, :));    % mean of last three days
x(2, :) = nanmean(dur(post, :));   % mean of last three nights
plot(x)
hold on
errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 2)
box off
set(gca, 'TickLength', [0 0])
xlim([0.8 2.2])
ylabel('Running Duration [s]')
xticks([])
[~, p, ~, stats] = ttest(x(1, :), x(2, :));
txt = sprintf('t = %.2g; p = %.2g', stats.tstat, p);
title(txt)

subplot(3, 1, 3)
x(1, :) = nanmean(weight(pre, :));    % mean of last three days
x(2, :) = nanmean(weight(post, :));   % mean of last three nights
plot(x)
hold on
errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 2)
box off
set(gca, 'TickLength', [0 0])
xlim([0.8 2.2])
ylabel('Norm. Weight')
xticks([1 2])
xticklabels({'Pre', 'Post'})
[~, p, ~, stats] = ttest(x(1, :), x(2, :));
txt = sprintf('t = %.2g; p = %.2g', stats.tstat, p);
title(txt)
sgtitle('Pre- vs. Post-Surgery')


%%% w/ w/o EtOH
clear x
post = 25;
pre = 20 : 1 : 24;

figure
subplot(2, 1, 1)
x(1, :) = nanmean(correct(pre, :));    
x(2, :) = correct(post, :);   
plot(x)
hold on
errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 2)
box off
set(gca, 'TickLength', [0 0])
xlim([0.8 2.2])
ylabel('Success Rate [%]')
xticks([])
[~, p, ~, stats] = ttest(x(1, :), x(2, :));
txt = sprintf('t = %.2g; p = %.2g', stats.tstat, p);
title(txt)

subplot(2, 1, 2)
x(1, :) = nanmean(dur(pre, :));    % mean of last three days
x(2, :) = dur(post, :);   % mean of last three nights
plot(x)
hold on
errorbar(nanmean(x, 2), std(x, [], 2, 'omitnan'), 'k', 'LineWidth', 2)
box off
set(gca, 'TickLength', [0 0])
xlim([0.8 2.2])
ylabel('Running Duration [s]')
xticks([])
[~, p, ~, stats] = ttest(x(1, :), x(2, :));
txt = sprintf('t = %.2g; p = %.2g', stats.tstat, p);
title(txt)
sgtitle('EtOH')
