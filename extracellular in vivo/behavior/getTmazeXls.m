
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepath = 'D:\Google Drive\PhD\Slutsky\Data Summaries';
filename = 'b3.xlsx';
filename = fullfile(basepath, filename);
mnames = [60 61 62 63 66 67];

forceL = false;
saveFig = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and arrange data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nmice = length(mnames);
if forceL || ~exist('micenames', 'var')
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
                end
            end
        end
    end
    
    % read weight
    [num, ~, raw] = xlsread(filename, 'Weight');
    BLweight = num(1, 2 : end);
    
    % read water intake
    [num, ~, raw] = xlsread(filename, 'WaterIntake');
    water = num(:, 2 : nmice + 1);
    waterDate = num(:, 1);
end

% find max number of sessions
maxsessions = 1;
for i = 1 : nmice
    maxsessions = max([maxsessions, length(m{i}.date)]);
end

% arrange mats (session x mouse)
correct = nan(maxsessions, nmice);
dur = nan(maxsessions, nmice);
weight = nan(maxsessions, nmice);
for i = 1 : nmice
    nsessions = length(m{i}.date);
    correct(1 : nsessions, i) = nanmean(m{i}.correct, 1);
    dur(1 : nsessions, i) = nanmean(m{i}.dur, 1);
    weight(1 : nsessions, i) = m{i}.weight / BLweight(i) * 100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - params vs. session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mice = [1 : 6];
% mice = [2, 6];
delay = 10;
dates = [];
% dates = [2.9 : 9.9];
dates = [2.1 : 14.1];
shadeIdx = [5 10];
shaedIdx = [];
% shadeIdx = [8.5 9.5];
mClr = 'nnn';
pflag = [1 1 1 1];      % (1) success, (2) duration, (3) weight, (4) water
avgFlag = [1];
% avgFlag = 1;
cump = cumsum(pflag);

% select specific data
clear data
for i = mice
    nsessions = length(m{i}.date);
    delayIdx = m{i}.delay == delay;
    dateIdx = zeros(1, nsessions);
    if isempty(dates)
        dates = unique(m{1}.date);
    end
    for ii = 1 : length(dates)
        dateIdx = dateIdx | m{i}.date == dates(ii);
    end
    idx = delayIdx & dateIdx;
end
waterIdx = zeros(length(waterDate), 1);
for i = 1 : length(dates)
    waterIdx = waterIdx | dates(i) == waterDate;
end

close all
% -------------------------------------------------------------------------
% success
if pflag(1)
    fh = figure;
    data = correct(idx, mice);
    yLimit = [0 1];
    subplot(sum(pflag), 1, cump(1))
    hax(cump(1)) = gca;
    paramsVsSession(data, yLimit, mClr, shadeIdx, avgFlag, dates, mice)
    ylabel('Success Rate [%]')
    yticks([0 0.5 0.75 1])
end

% -------------------------------------------------------------------------
% running duration
if pflag(2)
    data = dur(idx, mice);
    yLimit = [0 ceil(max(max(data)))];
    subplot(sum(pflag), 1, cump(2))
    hax(cump(2)) = gca;
    paramsVsSession(data, yLimit, mClr, shadeIdx, avgFlag, dates, mice)
    ylabel('Trial Duration [s]')
end

% -------------------------------------------------------------------------
% weight
if pflag(3)
    data = weight(idx, mice);
    yLimit = [50 100];
    subplot(sum(pflag), 1, cump(3))
    hax(cump(3)) = gca;
    paramsVsSession(data, yLimit, mClr, shadeIdx, avgFlag, dates, mice)
    ylabel('Norm. Weight [%]')
end

% -------------------------------------------------------------------------
% water
if pflag(4)
    data = water(waterIdx, mice);
    yLimit = [0 ceil(max(max(data)))];
    subplot(sum(pflag), 1, cump(4))
    paramsVsSession(data, yLimit, mClr, shadeIdx, avgFlag, dates, mice)
    ylabel('Water Intake [ml]')
end

xlabel('Time [date]')
% legend('best', micenames(mice))
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
    mice = [1, 3, 4, 6];
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


