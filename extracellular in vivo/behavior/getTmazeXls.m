
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = 'E:\Leore';
filename = 'b4.xlsx';
filename = fullfile(basepath, filename);
mnames = [74 : 80];

forceL = false;
saveFig = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmice = length(mnames);
if forceL || ~exist('micenames', 'var')
    for i = 1 : nmice
        
        micenames{i} = ['lh', num2str(mnames(i))];
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
                end
            end
        end
    end
    
    % read weight
    [num, ~, raw] = xlsread(filename, 'Weight');
    BLweight = num(1, 2 : end);
    
    % read water intake
    %     [num, ~, raw] = xlsread(filename, 'WaterIntake');
    %     water = num(:, 2 : nmice + 1);
    %     waterDate = num(:, 1);
end

% find max number of sessions
maxsessions = 1;
for i = 1 : nmice
    maxsessions = max([maxsessions, length(m{i}.date)]);
    if length(m{i}.date) >= maxsessions
        dates = m{i}.date;
    end
end
maxdays = 1;
for i = 1 : nmice
    maxdays = max([maxdays, length(unique(m{i}.date))]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select specific data
mice = [1 : 7];
delay = 30;

% initialize mats (session x mouse)
correct = nan(maxsessions, nmice);
dur = nan(maxsessions, nmice);
weight = nan(maxsessions, nmice);
for i = mice
    nsessions = length(m{i}.date);
    delayIdx = find(m{i}.delay == delay);
    for ii = 1 : length(delayIdx)
        correct(find(dates == m{i}.date(delayIdx(ii)), 1), i) =...
            nanmean(m{i}.correct(:, delayIdx(ii)), 1);
        dur(find(dates == m{i}.date(delayIdx(ii)), 1), i) =...
            nanmean(m{i}.dur(:, delayIdx(ii)), 1);
        weight(find(dates == m{i}.date(delayIdx(ii)), 1), i) =...
            nanmean(m{i}.weight(:, delayIdx(ii)), 1);
    end
end

% remove nan
rmidx = find(all(isnan(correct), 2));
correct(rmidx, :) = [];
dur(rmidx, :) = [];
weight(rmidx, :) = [];
% waterIdx = zeros(length(waterDate), 1);
% for i = 1 : length(dates)
%     waterIdx = waterIdx | dates(i) == waterDate;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - params vs. session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shadeIdx = [5 10];
shadeIdx = [];
mClr = 'nnn';
pflag = [1 1 0 0];      % (1) success, (2) duration, (3) weight, (4) water
avgFlag = [1];
cump = cumsum(pflag);

close all
% -------------------------------------------------------------------------
% success
if pflag(1)
    fh = figure;
    yLimit = [0 1];
    subplot(sum(pflag), 1, cump(1))
    hax(cump(1)) = gca;
    paramsVsSession(correct, yLimit, mClr, shadeIdx, avgFlag, dates, mice)
    ylabel('Success Rate [%]')
    yticks([0 0.5 0.75 1])
end

% -------------------------------------------------------------------------
% running duration
if pflag(2)
    yLimit = [0 ceil(max(max(dur)))];
    subplot(sum(pflag), 1, cump(2))
    hax(cump(2)) = gca;
    paramsVsSession(dur, yLimit, mClr, shadeIdx, avgFlag, dates, mice)
    ylabel('Trial Duration [s]')
end

% -------------------------------------------------------------------------
% weight
if pflag(3)
    yLimit = [50 100];
    subplot(sum(pflag), 1, cump(3))
    hax(cump(3)) = gca;
    paramsVsSession(weight, yLimit, mClr, shadeIdx, avgFlag, dates, mice)
    ylabel('Norm. Weight [%]')
end

% -------------------------------------------------------------------------
% water
if pflag(4)
    yLimit = [0 ceil(max(max(data)))];
    subplot(sum(pflag), 1, cump(4))
    paramsVsSession(water, yLimit, mClr, shadeIdx, avgFlag, dates, mice)
    ylabel('Water Intake [ml]')
end

xlabel('Time [date]')
legend('best', micenames(mice))
xticks([1 : length(dates)])
xticklabels(split(num2str(dates)))

if saveFig
    figname = fullfile(basepath, 't-maze');
    print(fh, figname, '-dpdf', '-bestfit', '-painters');
    export_fig(figname, '-tif', '-transparent', '-r300')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% success vs. delay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select specific data
pflag = 0;
if pflag
    clear data
    mice = [1 : 7];
    rmDates = [3.9, 4.9, 5.9, 6.9];      % dates to exclude from calculation
    
    delays = unique(m{i}.delay);
    k = 1;
    for i = mice
        nsessions = length(m{i}.date);
        delayIdx = zeros(1, nsessions);
        dateIdx = zeros(1, nsessions);
        for ii = 1 : rmDates
            dateIdx = dateIdx | m{i}.date == rmDates(ii);
        end
        for ii = 1 : length(delays)
            delayIdx = m{i}.delay == delays(ii);
            idx = delayIdx & ~dateIdx;
            data(k, ii) = mean(correct(idx, i));
        end
        k = k + 1;
    end
    
    figure
    bar(mean(data, 1))
    hold on
    errorbar(mean(data, 1), std(data, [], 1), 'k', 'LineWidth', 2)
    xticklabels(split(num2str(delays)))
    xlabel('Delay [s]')
    ylabel('Success Rate [%]')
    set(gca, 'TickLength', [0 0])
    title('Mnemonic Difficulty')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% no. trials on variation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mice = [1 : 7];
delay = 00;
trials = 1 : 10;
% trials = [];

correct = nan(maxsessions, nmice);
dur = nan(maxsessions, nmice);
weight = nan(maxsessions, nmice);
if isempty(trials)
    trials = ':';
end
for i = mice
    nsessions = length(m{i}.date);
    delayIdx = find(m{i}.delay == delay);
    for ii = 1 : length(delayIdx)
        correct(find(dates == m{i}.date(delayIdx(ii)), 1), i) =...
            nanmean(m{i}.correct(trials, delayIdx(ii)), 1);
    end
end
rmidx = find(all(isnan(correct), 2));
correct(rmidx, :) = [];
dur(rmidx, :) = [];
weight(rmidx, :) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data vs session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function paramsVsSession(data, yLimit, mClr, shadeIdx, avgFlag,...
    ~, mice)
if avgFlag
    if length(avgFlag) > 1
        [~, ia] = intersect(mice, avgFlag);
        stdshade(data(:, ia)', 0.1, 'k');
    else
        stdshade(data', 0.1, 'k');
    end
end
hold on
if ~isempty(shadeIdx)
    fill([shadeIdx fliplr(shadeIdx)]', sort([yLimit yLimit]), 'r',...
        'FaceAlpha', 0.1, 'EdgeAlpha', 0)
end
hp = plot(data);
if length(mice) == length(mClr)
    for i = 1 : length(mice)
        if mClr(i) == 'n'
            hp(i).Color = 'none';
        else
            hp(i).Color = mClr(i);
            hp(i).Color(4) = 1;
        end
    end
end

% plot(x, [0.5 0.5], '--k', 'LineWidth', 1)
ylim(yLimit)
ylabel('Success Rate [%]')
x = xlim;
xlim([1 x(2)])
xticks([])
box off
set(gca, 'TickLength', [0 0])
end


