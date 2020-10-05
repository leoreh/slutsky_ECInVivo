
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
    
    micenames{i} = ['m', num2str(mnames(i))];
    [num, ~, raw] = xlsread(filename, micenames{i});
    m{i}.name = micenames{i};
    m{i}.date = num(1, 2 : 2 : end);
    m{i}.time = num(2, 2 : 2 : end);
    m{i}.weight = num(3, 2 : 2 : end);
    m{i}.delay = num(4, 2 : 2 : end);
    m{i}.dur = num(5 : end, 2 : 2 : end);
    
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

%%% state metrics per session 
for i = 1 : nmice
    correct(:, i) = nanmean(m{i}.correct, 1);
    dur(:, i) = nanmean(m{i}.dur, 1);
    weight(:, i) = m{i}.weight / BLweight(i);
end

% select specific data
date = 10.9;
didx = m{i}.date == date;
idx = didx & sidx;

data = correct(didx, mice)'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mice = [1 : 6]; 
c = 'kkkkkr';
delay = 00;
clear data
for i = mice
    idx = find(m{i}.delay == delay);
    data(:, i) = correct(idx, i);
end
mean(data, 2)
mean(data(end : -1 : end - 2, :), 1)

figure
subplot(2, 1, 1)
hp = plot(data);
% for i = 1 : length(mice)
%     hp(i).Color = c(i);
% end
hold on
stdshade(data', 0.3, 'k');
ylim([0 1])
x = xlim;
xlim([1 x(2)])
plot(x, [0.5 0.5], '--k', 'LineWidth', 1)
yticks([0 0.5 0.75 1])
xticks([])
box off
set(gca, 'TickLength', [0 0])
ylabel('Success Rate [%]')

data = dur(sidx, mice); 
subplot(2, 1, 2)
hp = plot(data);
for i = 1 : length(mice)
    hp(i).Color = c(i);
end
hold on
stdshade(data', 0.3, 'k');
y = ylim;
ylim([0 y(2)]);
box off
set(gca, 'TickLength', [0 0])
ylabel('Trial Duration [s]')
xlabel('Time [d]')

mean(data(end - 1, :))





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

